# AMSS-NCKU DevPod 调试日志

> 本文件记录每次优化在 **DevPod** 上的调试运行:验证程序能编译、能跑完、正确性 PASS。
> **不用于计时**(devpod 共享 64 核,时间不准),只验证正确性。
> 集群正式计时见 [OPTIMIZATION_LOG.md](OPTIMIZATION_LOG.md)。
>
> **perf 采样过程也记录在本文件**(采样命令、参数、产出文件、job ID)。

## DevPod 调试参数(固定)
- 环境:鲲鹏 920B devpod, aarch64, TaiShan-v120, 64 核共享, 4 NUMA
- 演化:`Final_Evolution_Time = 4.0`(4 步快速验证)
- MPI:`MPI_processes = 30`(与生产一致)
- OMP:与本次优化设置一致
- 用 `./run.sh --twop-cache`(初值不变时跳过 TwoPuncture,加速调试)
- ⚠️ 若本次优化**改了 TwoPuncture 相关代码**,不能用 cache

## 提交坑
`hpc submit` 用提交时 cwd。devpod 上直接 `cd /home/h3250104945/HPC101/src/lab4 && ./compile.sh && ./run.sh --twop-cache && ./check.sh`。

## perf 采样限制(重要,避免重复踩坑)
- **计算节点**:`perf stat` 可用;`perf record` 受 `perf_event_mlock_kb=516`(RO)限制。
  - `-m` 太大(1024)→ "Permission error mapping pages"
  - `-m` 太小(1)→ 100% 丢样本
  - **可行**:MPI=1 + `-m 64`(job 104167 成功,81638 样本)
- **计算节点 perf report 无法解析符号**(刷屏 `unwind: get_proc_name unsupported`)→ 原始 perf4.data 要回 devpod 重解析:`perf report -i perf4.data --stdio ... | grep -vE "unwind|get_proc_name"` 落盘干净版。
- **devpod**:perf record 可用,但 30 rank 4 步要 >40min(共享机太慢)。
- **多 rank record 已可行**：30 rank 使用 `-m 4` 已成功（见 perf 采样 #4/#5）；更大的 mmap 缓冲仍会超过 516KB 的全局 mlock 预算。

## 记录模板

### DevPod 调试记录(每次优化复制)
```markdown
## DevPod 调试 #N: <简短标题>  (YYYY-MM-DD)
- 优化内容: <改了什么, 对应 OPTIMIZATION_LOG.md 优化 #N>
- 编译: 成功/失败(错误信息)
- 运行: 能跑完4步? <输出摘要>
- 正确性(check.sh): RMS=?%, constraint max=?  (PASS/FAIL)
- 备注: <调试中发现的问题>
```

### perf 采样记录(每次采样复制)
```markdown
## perf 采样 #N: <目的>  (YYYY-MM-DD)
- 目标: <看什么, 如 "30rank通信占比" / "MPI=1函数热点">
- 采样位置: 计算节点/devpod
- 参数: MPI=?, 步数=?, `perf record -F 99 -m ? --call-graph dwarf`
- job ID / 产出: perfX.data, perfX_*.txt
- 结果摘要: <样本数, 是否丢失, 关键发现>
- 备注: <成功/失败原因>
```

---

## perf 采样 #1: baseline 函数热点  (2026-08-16)
- 目标: 函数级热点(MPI=1,无通信干扰,看纯计算分布)
- 采样位置: 计算节点(job 104167, profile4.sh)
- 参数: MPI=1, 步数=4, `perf record -F 99 -m 64 --call-graph dwarf`,无 cache(完整链路)
- 产出: perf_profiles/perf4.data (674M), perf4_self/symbols/total.txt(已回 devpod 重解析干净版)
- 结果摘要: 81638 样本,0 丢失。Top: compute_rhs_bssn 16.7%, polint_ 9.9%, LineRelax_be 9.6%, __cos 8.4%, kodis_ 8.0%
- 备注: 计算节点生成 report 不解析符号,回 devpod 重解析才得到干净版。MPI=1 单步 ~131s(慢于 30rank 的 43s,因为没并行)。

## perf 采样 #2: baseline 硬件计数器  (2026-08-16)
- 目标: 30 rank 宏观硬件行为(IPC/cache/利用率)
- 采样位置: 计算节点(job 103849, profile.sh, perf stat)
- 参数: MPI=30, 步数=4, `perf stat -d -d -d`,无 cache(完整链路)
- 产出: perf_profiles/perf_stat.txt
- 结果摘要: wall=469s。IPC=1.93, LLC-miss=39.41%, dTLB-miss=6.87%, page-faults 4.4M。CPUs utilized = 11.7(⚠️ 全程平均,不代表演化阶段,见备注)。
- 备注: stat 用计数器不需 mmap,所以 30rank 随便跑。⚠️ **曾误读为"有效CPU 11.7/30 = 60% 资源空等",已修正**:11.7 是全程平均(perf stat attach 整个 run.sh,含 ~287s 单进程 matplotlib 绘图;单进程 CPU 有硬上限 ≈ 其墙钟)。排除法:总 CPU 5511s − 绘图 ≤287s ≈ 5224s 只可能来自 30 个 ABE 进程,5224/173.6s(演化墙钟)≈ 30.1 核 → **演化阶段实际满核,不是 60% 空等**。真正的问题是满核时 ~76% CPU 周期在 opal_progress 忙等轮询(见采样 #4)。

## perf 采样 #3: 尝试 30rank record(失败)  (2026-08-16)
- 目标: 想拿 30rank 函数级热点(看 MPI_Waitall 占比)
- 采样位置: 计算节点(job 104118 profile3.sh: -m 1; 多次尝试)
- 参数: 30rank + 大 `-m` 或 `-m 1`
- 结果: 全失败。大 `-m` 报 Permission error;`-m 1` 100% 丢样本。
- 备注: 待用折中方案(4-8 rank + 中等 -m)重试。

---

## perf 采样 #4: 30rank record 成功(perf5,决定性证据)  (2026-08-16)
- 目标: 拿到 30 rank 的函数级热点,确认 MPI 忙等占比
- 采样位置: 计算节点(job 107382, profile5.sh)
- 参数: MPI=30, 步数=4(--twop-cache), `perf record -F 99 -m 4 --call-graph dwarf`
- 产出: perf_profiles/perf5_m4.data (4.2G), perf5_m4.log, perf5_symbols.txt, perf5_self.txt
- 结果摘要: 520013 样本(丢 760 个)。**~76% 样本在 MPI 库**(57.71% libopen-pal 轮询循环 + 5.53%/5.23% libmpi + 4.68% 其他),ABE 计算符号合计 <10%(compute_rhs_bssn 6.92% 最大)。
- 备注: `-m 4`(30核×16KB=480KB < 516KB mlock 预算)是唯一可行配置。计算节点 report 不解析符号,回 devpod 重解析。**忙等是用户态执行,占 CPU → 演化阶段满核(采样#2 已证)时,76% ≈ 每步 43s 里 ~33s 在 opal_progress 轮询**。
- ⚠️ **本次运行输出数据/绘图是坏的**:perf record 写 4.2GB 致 /workspace IO 过载(lost 253 chunks),ABE 输出文件可能写失败,plot 画空图仅 ~1s(正常绘图 ~110s,见采样#1 的 perf4 占比拆分)。**该次结果不能用于正确性检查;perf 演化采样本身可信**。

---

---

## DevPod 调试 #2: OpenMP 并行正确性验证  (2026-08-17)
- 优化内容: 热点 Fortran stencil 循环加 `!$omp parallel do`(diff_new.f90 11 处 + kodiss 1 处 + lopsidediff 1 处,commit 1b9c77d);`AMSS_ENABLE_OPENMP=ON` 编译
- 编译: 成功(devpod,g++/gfortran + OpenMP)
- 运行: MPI=9 4 步,OMP=1 → 527.2s(单步 131.8s);OMP=4 → 926.5s(单步 231.6s,**devpod 上并行反而慢 75%**——共享机噪声 + 可能并行开销,集群为准)
- 正确性(check.sh): **OMP=1 和 OMP=4 都 FINAL: PASS,约束值逐位相同**(Ham=0.22588831, Px=0.02448088, Py=0.0091883125, Pz=0.014353656)→ 并行不改变数值结果 ✓
- 备注: Input.py 已恢复(MPI=30, OMP=1)。集群 OMP 组合实验: job 107588(30:1 / 15:2 / 10:3 / 6:5)

---

## DevPod 调试 #1: rank_sweep.sh 脚本验证(MPI=9)  (2026-08-17)
- 优化内容: 验证 rank 扫描脚本 `rank_sweep.sh`(备份/恢复 Input.py、4 步 + twop-cache、提取演化时间)——对应集群优化 #1(MPI rank 数探索)
- 编译: 成功(devpod 256 核 4 NUMA,compile.sh 全项目)
- 运行: MPI=9 跑完 4 步,Total Evolve Time = 541.5s(devpod 共享机,仅验证流程,时间不作数)
- 正确性: 未单独 check(集群正式计时时统一 check)
- 备注: **Input.py 通过 trap EXIT 成功恢复为 MPI_processes=30**;脚本可安全用于集群。顺带采集 devpod 拓扑:256 核 4 NUMA,node 0-1 近(10-12)、0↔2/3 远(35-40),mpiexec 默认无绑定

<!-- DevPod 调试 + perf 采样记录追加在此线之下 -->

## DevPod 调试 #3: 忙等根因诊断插桩  (2026-08-17)
- 优化内容: Parallel::transfer 累计 Waitall 时间 + bssn_class::Step 按步打印 wall/tfer/ana 分解(commit 47989cd + 7d1f777)——对应集群优化 #3(球面插值 OpenMP)的诊断阶段
- 编译: 成功
- 运行(集群 job 108473,30 rank 4 步):每粗步 42.8s 里 `STEP lev=0 YN=1 wall=32.9 tfer=0.12 ana=32.89` → **transfer 仅 0.12s,ana(分析)占 32.89s**
- 根因定位: perf5 按 pid 分析,30 rank CPU 时间"均衡"(3.32-3.33%)是假象——29 个等待 rank 在 opal_progress 忙等轮询烧 CPU,掩盖了同步点到达时间差。Allreduce 等待即负载不均衡:1 个 rank 串行做球面插值(polint 1.49% 全在 pid 377=rank 4),29 个到同步点干等 ~33s → 这就是 perf5 的 76% opal_progress
- 量化: 40 步里分析 ≈ 40×32.9 = 1316s,占演化时间 77%,比计算(10.3s/步)大 3 倍
- 结论: 瓶颈不是 transfer、不是 rank 间计算不均,是分析阶段单 rank 串行插值 + Allreduce 同步。→ 优化 #3(rank 内 OpenMP)先缓解,优化 #4(跨 rank 分散)根治

---

## DevPod 调试 #4: 分布式分散正确性验证  (2026-08-17)
- 优化内容: Interp_Points CPU 路径(NN≥256)改为广播中心块数据 + 按 rank 切片插值(对应集群优化 #4)
- 编译: 成功(devpod,g++/gfortran)
- 运行: MPI=30 OMP=1 4 步 + twop-cache,Total Evolve Time = 449.3s(单步 ~13.9s,devpod 共享机,时间不作数)
- STEP ana 分解: `wall=14.0 tfer=2.7 ana=13.7`(devpod 共享核争抢,ana 未到理论 ~1s;集群独占核后 ana 1.4s,见优化 #4)
- 正确性(check.sh): **FINAL: PASS**,约束值 Ham=0.22588831/Px=0.02448088/Py=0.0091883125/Pz=0.014353656 与 30:4(#3)逐位完全一致 → 切片+SUM 还原等价于原单 rank 全算 ✓
- 踩坑: 首版用"找含原点块"匹配球壳,但 PatList_Interp_Points(单点 BH 探针,NN=1)也走该路径,点不在原点 → abort。修复:加 `NN≥256` 门控,小查询走原路径;"找含原点"改"找含第一个插值点的块"

---

## perf 采样 #5: 分布式分散后热点复查  (2026-08-17)
- 目标: 验证分散后忙等是否消失、定位新热点
- 采样位置: 计算节点(job 110600, profile5.sh)
- 参数: MPI=30, OMP=1, 4步+cache, `perf record -F 99 -m 4 --call-graph dwarf`(分布式分散代码,Input 临时 4.0)
- 产出: perf_profiles/perf6_m4.data (1.2G, 从 perf5_m4 重命名避免覆盖 baseline)、perf6_m4.log、perf6_symbols.txt、perf6_addr.txt、perf6_dso.txt(后三者回 devpod 重解析,`grep -vE 'unwind|get_proc_name'`)
- 样本: 141K,丢 743
- **DSO 归类**:
  | DSO | perf6 分散后 | perf5 baseline |
  |---|---|---|
  | ABE | 53.92% | <2% |
  | libopen-pal | 27.64% | ~76% |
  | libc | 12.17% | ~3% |
  | libmpi | 5.13% | ~5% |
- **函数符号**: compute_rhs_bssn 25.48%、0x00f13c4(未解析,在 libopen-pal)22.96%、polint 5.40%、__memcpy_sve 4.81%、kodis 4.80%、fdderivs 4.57%、lopsided 3.98%、__memset_sve_zva64 3.37%、prolong3 2.34%、opal_progress 0.32%
- **关键结论**:
  - ✅ 分散成功:ABE 计算从 <2% 回升到 53.92%,忙等腾出的 CPU 给了有效计算
  - ⚠️ 忙等没消失,从 ~76% 降到 27.64%。**分析阶段(ana)忙等消失了**(ana 23.9→1.4s),但**计算阶段同步等待仍在**:RecursiveStep 各层 RestrictProlong→transfer 集体同步、RK4 子步间 Sync、新加 16 次 Bcast,30 rank 通信点快等慢仍产生 opal_progress 轮询
  - 下一瓶颈:计算阶段 AMR 层间通信的集体同步(27.64% libopen-pal),要降需动 transfer/RestrictProlong 通信结构(方向 5)
- 备注: ⚠️ profile5.sh 输出文件名写死 perf5_m4.*,本次覆盖了 baseline 的 perf5_m4.data 原始数据(baseline 只剩 perf5_symbols.txt/perf5_self.txt 干净版文本)。已重命名为 perf6_*。以后 perf 采样脚本应按次编号输出,避免覆盖

---

## DevPod 调试 #5: POINTWISE 模板阶数提升验证  (2026-08-18)
- 优化内容: 每个 bulk 点先计算一次 `pw_order`,内联一阶和二阶模板复用这个结果,对应集群优化 #11。
- 编译: POINTWISE + OpenMP 编译成功。
- 运行: MPI=30,4 步,`--twop-cache`；输出 `pointwise_runs/pointwise_order2_smoke_20260818_a`。
- 正确性: **FINAL: PASS**,Trajectory RMS=0,约束值与 golden 逐位一致。

## DevPod 调试 #6: 二阶导点内标量化验证  (2026-08-18)
- 优化内容: POINTWISE 的 66 个二阶导结果不再写入三维数组,改为点内标量,移位项和 Ricci 相关计算直接消费,对应集群优化 #12。
- 编译: POINTWISE 和 HALO 均编译成功；LEGACY 预处理也通过。
- 运行: MPI=30,2 步,`--twop-cache`；输出 `pointwise_runs/pointwise_d2scalar_smoke_20260818_b`。
- 正确性: **FINAL: PASS**,Trajectory RMS=0；集群 4 步 perf 输出再次检查为 **FINAL: PASS**。
- 踩坑: 首版漏改移位向量二阶导的公共消费点,出现 NaN；补齐 18 个标量消费后恢复正确。
