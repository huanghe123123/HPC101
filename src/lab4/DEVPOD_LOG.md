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
- **多 rank record 尚未做**(30 rank 缓冲超限;4/8 rank 折中方案待用)。

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

<!-- DevPod 调试 + perf 采样记录追加在此线之下 -->
