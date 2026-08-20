# Lab4 工作交接文档

> 最后更新：2026-08-20
> 项目：HPC101 实验四 AMSS-NCKU 数值相对论程序优化（CPU 任务一）
> 当前有效代码基线：优化 #42（BSSN #26 + #37 + fast-math #41 + global_interp uniform6 #42 + TwoPuncture #30–#35）

这份文件是继续工作前的唯一状态入口。它区分“当前有效代码”“已完成但被拒绝的实验”和“仍未测试的方向”，不能把历史实验配置当成当前可执行配置。

## 1. 当前结论

- 当前有效最优是 #42（BSSN #26 + #37 + fast-math #41 + global_interp uniform6 #42 + TwoPuncture #30–#35）：40 步 `This Program Cost=342.549s`，`Total Evolve=329.018s`，`Total Running=332.715s`，`FINAL: PASS`，匹配 golden 40/100 组，Trajectory RMS=4.05e-6%，Constraint 全部 PASS。相对 #37，端到端下降约15.0%、演化下降约15.8%。
- 阶段一（#22–#25）40 步检查点 job 115838 为 `681.446s`、Evolve `381.533s`、`FINAL: PASS`；阶段二加入 #26 后的正式结果以 job 115969 为准。
- #22/#23 消除发布 NaN 全扫描并让细层 predictor 走演化路径；#24 pooled 111 个 legacy workspace 数组；#25 编译掉真空物质数组加载；#26 融合 B=2 lopsided+KO 尾部扫描，全部通过短测并累计接受。
- #14–#17 都已在集群完成 paired `perf stat`、`perf record` 和 4 步正确性检查；四项均通过数值检查，但都没有达到稳定的性能接受标准，因此全部回退。
- #19（RHS 运算顺序重排）仍被拒绝；#27 简单 SIMD directive、#28 cache coloring、#29 j×k tile sweep 也已完成 paired 4 步实验但没有稳定收益，全部回退。#28 的首个 slab 实现曾 SIGSEGV，后改为独立 colored allocation 后正确性通过但性能仍回退。
- #38（fused RHS 主循环 hybrid OpenMP，`AMSS_RHS_OMP_ASSEMBLY`）正确性通过（private 列表完整：30×1 与 15×2 开关 ON vs OFF 均 bit-identical，无竞争），但性能负向被拒绝：30×1 ON 回退 +5%（超 ±0.5% 门），15×2 端到端慢 2.8×。阶段C SMT sweep（30×2/15×4/10×6/6×10）亦负向：30×2 比 #37 OFF 慢 3.3%（SMT 第二线程在 compute-bound 双精度 stencil 上无贡献，sibling 共享 L1/L2 无法掩盖 LLC miss），减 rank 配置慢 2.2–2.8×。hybrid MPI+OpenMP/SMT 方向整体放弃，开关默认 OFF，#37 不变。
- #39（THP slab，`AMSS_FGFS_HUGEPAGE_SLAB`）正确性通过（RMS=0、constraint 逐位一致），但性能严重负向被拒绝：Total Evolve +50% 回退（38.54→57.83s）。根因是每字段对齐 2MiB 把工作集从 40–141MiB 撑到 334MiB/Block，L3 装不下；字段间隔 2MiB 空洞破坏原有空间局部性，预取失效。dTLB miss 未降（3.54%→3.50%），LLC miss rate 降是慢跑的假象（绝对数 +10%）。大页对工作集远大于 L3 的流式场景是负优化。开关默认 OFF，#37 不变。
- #41（CPU ABE `-ffast-math`）通过正式 40 步：`Total Evolve=379.041s`、`Total Running=381.984s`、`This Program Cost=390.893s`、`FINAL: PASS`，相对 #37 端到端约快3.0%。它主要减少指令和 cycles，LLC miss 绝对数基本不变，因此不把它误判成访存优化；`AMSS_FAST_MATH` 已默认开启，TwoPuncture 和 GPU target 不受影响。
- #42（`AMSS_POLINT_UNIFORM6`）将 `global_interp→polin3` 的固定 `ordn=6` 路径改为专用均匀网格 tensor-product Lagrange 核；AMR Restrict/Prolong、`interp_2` 和 Sommerfeld 仍保留通用 `polint/polin3`。4-step paired 中 `Total Evolve` 为 39.41/39.42→33.14/33.52s，两个版本均 `FINAL: PASS`；正式 40 步为 `Total Evolve=329.018s`、`This Program Cost=342.549s`、`FINAL: PASS`。`polint_` self 从9.05%降至0.14%，开关已默认开启。
- 阶段一 miss 归因（诊断，不入成绩）：`perf record -e LLC-load-misses` 直接采样 miss 事件，推翻了 LLC/dTLB 计划阶段三的前提——LLC miss 最大来源是 `rungekutta4_rout_` 占 31.2%（whole-array triad，算法固有），不是 RHS（12.6%）+lopsided（10.9%）；dTLB 真实访存主源是 compute_rhs 占 28%（51% 是 opal 忙等假象）。RHS/lopsided 不是 LLC 主源，阶段三 derivative tile 暂缓。
- 当前源码正式路径为 #42；TwoPuncture 使用连续 JFD、预计算谱/几何表、NRELAX=10 和30线程 OpenMP。OpenMP 初版曾因共享临时量 SIGSEGV，已修复并重新验证；未修复的旧二进制不得使用。
- TwoPuncture 固定开销已不再是主导；`polint_` 已被 global_interp 专用核压低，当前主要剩余成本是 `Total Evolve=329.018s` 中的 RHS、batch stencil 和 AMR 层级等待。

## 2. 当前有效配置

### 2.1 源码与编译

当前有效路径由以下已有优化组成：

1. #8 `compute_rhs_bssn_fused_`：组装部分显式单点循环，点内中间量使用标量，去掉大量组装临时数组。
2. #11 `pw_order`：每个点只计算一次模板阶数，所有一阶/二阶模板复用结果。
3. #12：POINTWISE 二阶导结果使用点内标量。
4. #13：POINTWISE 一阶导结果使用点内标量。
5. #18 `batch_stencils`：kodis/lopsided 批处理（B=2），把 B 个场打包进一次扫描，一次 halo setup，每点一次模板阶数判断，lopsided 每点只读一次 shift。
6. #22/#23：关闭发布版 RHS NaN 全扫描；细层 predictor 使用 POINTWISE 演化路径，level 0 predictor 保留约束模式。
7. #24：legacy/constraint 分支的111个导数数组使用按 Block shape 精确复用的 workspace pool，POINTWISE 热路径不逐次创建这些数组。
8. #25：真空 BSSN 编译掉 `rho/S_i/S_ij` 恒零加载，保留原始 ABI 和字段布局。
9. #26：B=2 路径把 lopsided 与 KO 合并为一次 RHS 扫描。
10. #30：TwoPuncture JFD 连续布局、line workspace 复用、原地 Thomas 和 JFD 点级栈临时量。
11. #31：Chebyshev/Fourier 谱表、A/B/phi 几何表缓存，并跳过已收敛质量搜索的冗余 final Newton。
12. #32：TwoPuncture `NRELAX=10`，保留失败时回退到200的通用路径。
13. #33：TwoPuncture 独立30线程 OpenMP；谱导数、J-times 和红黑 line relaxation 使用线程私有 workspace。

正式有效配置使用：

```text
MPI=30
OMP=1
GPU=no
AMSS_FAST_MATH=ON
AMSS_POLINT_UNIFORM6=ON
AMSS_TWOP_OPENMP=ON
AMSS_TWOP_THREADS=30
AMSS_TWOP_NRELAX=10
AMSS_TWOP_PARALLEL_DERIVS=ON
AMSS_TWOP_PARALLEL_J=ON
AMSS_TWOP_PARALLEL_RELAX=ON
AMSS_ENABLE_OPENMP=ON
AMSS_RHS_LOCALITY=POINTWISE
AMSS_BATCH_STENCIL=ON
AMSS_BATCH_STENCIL_B=2
AMSS_OPT=-O3
Final_Evolution_Time=40.0
TwoPuncture cache=OFF
AMSS_RHS_NAN_CHECK=OFF
AMSS_FINE_PREDICTOR_EVOLVE=ON
AMSS_RHS_WORKSPACE_POOL=ON
AMSS_VACUUM_BSSN=ON
AMSS_FUSED_RHS_TAIL=ON
AMSS_RHS_BULK_SIMD=OFF
AMSS_RHS_SKIP_CONSTRAINT_STORES=OFF
AMSS_RHS_SKIP_CHIN1_SCAN=ON
AMSS_RHS_RAW_DIAG_LOPSIDED=ON
AMSS_FGFS_COLORED_SLAB=OFF
AMSS_RHS_TILE_J=0
AMSS_RHS_TILE_K=0
```

注意：仓库中的 `AMSS_NCKU_Input.py` 当前是 4.0 的 smoke 配置，具体为 `MPI_processes=30`、`OMP_threads=1`、`GPU_Calculation="no"`、`Final_Evolution_Time=4.0`。正式 40 步脚本会临时改写该值并在退出时恢复，不能把当前输入文件的 4.0 当作正式成绩。

`CMakeLists.txt` 的默认值现在就是 CPU 正式路径（POINTWISE、B=2 batch、BSSN/TwoP OpenMP、#22–#26、#37、#41 和 #42 开关均开启；#36 保持关闭；GPU 仍默认关闭）。需要测试旧路径或负向候选时再显式传入覆盖开关：

```bash
./compile.sh -DAMSS_OPT="-O3" \
  -DAMSS_TWOP_OPENMP=ON \
  -DAMSS_RHS_LOCALITY=POINTWISE \
  -DAMSS_BATCH_STENCIL=ON \
  -DAMSS_RHS_NAN_CHECK=OFF \
  -DAMSS_FINE_PREDICTOR_EVOLVE=ON \
  -DAMSS_RHS_WORKSPACE_POOL=ON \
  -DAMSS_VACUUM_BSSN=ON \
  -DAMSS_FUSED_RHS_TAIL=ON
```

并在运行前设置 `AMSS_ENABLE_OPENMP=ON`、`OMP_NUM_THREADS=1`、`AMSS_BATCH_STENCIL_B=2`、`AMSS_TWOP_THREADS=30`、`AMSS_TWOP_NRELAX=10`。`AMSS_LEVEL_DIAG`、`AMSS_CACHE_COMM_PLANS` 和 `AMSS_ASYNC_LEVEL_SYNC` 当前均默认关闭；后两项不是当前正式基线的一部分。

### 2.2 当前正式结果

| 编号 | 方向 | 4 步/短测 | 40 步正式结果 | 状态 |
|---:|---|---:|---:|---|
| #8 | fused RHS 单点化 | — | 721.7s | 接受 |
| #9 | HALO 边界视图 | 约 47.00s | — | 无稳定收益，不作为独立结论 |
| #10 | POINTWISE 逐点导数 | 48.19s | — | LLC 改善但净变慢 |
| #11 | 模板阶数每点一次 | 46.20s | — | 接受为 #13 的组成部分 |
| #12 | 二阶导点内标量化 | 43.82s | — | 接受为 #13 的组成部分 |
| #13 | 一阶导点内标量化 | — | 698.59s，演化 411.75s | 被 #18 取代 |
| #14 | Step ERROR Allreduce 合并 | 55.27→57.42s | — | 拒绝 |
| #15 | 插值 Allgatherv | 42.64→42.97s | — | 拒绝 |
| #16 | 活动 rank communicator | 45.30→44.20s（stat 单次） | — | 采样反向，拒绝 |
| #17 | 通信计划/缓冲区缓存 | 43.03→43.00s | — | 仅 -0.08%，拒绝 |
| #18 | kodis/lopsided 批处理 | 44.13→41.15s（B=2） | 688.18s，演化 397.92s | 被 #26 取代，PASS |
| #19 | RHS 运算顺序重排 | 42.56→42.72s | — | bit-for-bit IDENTICAL 但 perf 持平，拒绝 |
| #22/#23 | NaN 扫描关闭 + 细层 predictor 演化路径 | 43.66→40.56s | 681.446s，演化 381.533s（阶段一） | 接受，PASS |
| #24 | legacy RHS workspace pool | 40.54→40.11s | 纳入阶段一 | 修复 shape 后接受 |
| #25 | 真空 BSSN 专用化 | 40.86→39.58s | 包含于阶段一 | 接受，PASS |
| #26 | lopsided+KO B=2 尾部融合 | 41.13→39.49s | 668.825s，演化 377.374s | 最终接受，PASS |
| #27 | fused RHS SIMD directive | 40.79→41.50s | — | 拒绝 |
| #28 | fgfs cache coloring | 38.31→41.79s（B=2） | — | 正确但明显变慢，拒绝 |
| #29 | j×k tile sweep | 38.80→38.83–39.41s | — | 无稳定 ≥1% 收益，拒绝 |
| #30 | TwoPuncture 连续 JFD/workspace/原地 Thomas | 286.79→277.46s | — | solver-only +3.25%，PASS |
| #31 | TwoPuncture 谱表/几何缓存 + final Newton 短路 | 275.70→160.97s | — | solver-only +41.6%，PASS |
| #32 | TwoPuncture NRELAX=10 | 160.02→69.84s | — | paired +56.4%，PASS |
| #33 | TwoPuncture 30线程 OpenMP | 69.88→8.47s | — | paired +87.9%，PASS |
| #34 | TwoPuncture full 4步链路 | 318.001→51.786s | — | `FINAL: PASS`，RMS=1.03e-5% |
| #35 | TwoPuncture正式40步 | — | 411.257s，演化397.011s | 历史基线，`FINAL: PASS` |
| #36 | `co=1` 跳过约束专用中间数组写回 | 57.840→57.435s（-0.70%） | job 116964，4步 `FINAL: PASS` | miss指标改善但总时间未达1%，保留实验开关，拒绝默认开启 |
| #37 | `chin1`/对角度规辅助扫描消除 | Total Evolve 两次约降1.32%/1.14% | job 120789，4步 `FINAL: PASS`；job 120900正式40步 `FINAL: PASS` | 接受，已被 #41/#42 叠加 |
| #41 | CPU `-ffast-math` | 40.168→38.544s（-4.05%，paired） | job 123942：390.893s，演化379.041s，`FINAL: PASS` | 接受，默认开启 |
| #42 | `global_interp` 固定 uniform6 `polin3` 核 | 39.41/39.42→33.14/33.52s（约-15%） | job 124023：342.549s，演化329.018s，`FINAL: PASS` | 接受，默认开启 |

## 3. #14–#17、#19 实验档案

所有候选都使用同一计算节点上的 baseline/candidate 配对运行，MPI=30、OMP=1、4 步、TwoPuncture cache。两个版本都执行 `check.sh`，结果均为 `FINAL: PASS`、Trajectory RMS=0。

### #14：Step 内 ERROR Allreduce 合并

- 任务：114175。
- 候选思路：把 predictor 和三个 corrector 的本地错误累积到 Step 末尾，每个 Step 只做一次规约。
- `perf stat`：`Total Evolve 55.2737s → 57.4241s`；task-clock +4.4%，cycles +3.2%，instructions +3.1%；LLC miss rate 57.06%→55.04%。
- total 调用树中 `PMPI_Allreduce` 热点消失，但 paired stat 没有墙钟收益。`perf record` 的第二次运行方向相反，不能用采样运行墙钟推翻 stat。
- 结论：拒绝并回退。

### #15：分布式插值改用 Allgatherv

- 任务：114272。
- 候选思路：将排序后的连续点切片紧凑收集，避免对完整零填充 `NN*num_var` 数组做 Allreduce。
- 候选实际执行 64 次大规模聚合，每次 `NN=36864`。
- `perf stat`：`Total Evolve 42.6377s → 42.9712s`；cycles +0.77%，instructions +0.07%；LLC miss rate 55.74%→56.99%。
- total 调用树中 `PMPI_Allreduce` 15.43%→14.90%，`PMPI_Allgatherv` 只有 0.49%；插值 collective 不是总关键路径。
- 结论：拒绝并回退。

### #16：每层活动 rank communicator

- 任务：114350。
- 候选思路：根据每层 Block owner 创建活动 rank communicator，让不持有该层数据的 rank 不参加该层 ERROR collective。
- 成员检查：level 0 只有 9 个活动 rank（mask `0x1ff`），level 1–8 均为 30 个。
- `perf stat` 单次应用计时 45.2985s→44.1994s，cycles -2.1%；但 `PMPI_Allreduce` 17.74%→17.75%，`opal_progress` 42.51%→43.75%，目标 collective 没有稳定下降。
- `perf record` 第二次运行反向（44.5332s→44.7841s），所以 stat 的单次下降视为噪声。
- 结论：拒绝并回退。

### #17：通信计划与工作缓冲区缓存

- 任务：114423。
- 候选思路：在已有 segment 计划基础上缓存 peer 长度、send/receive 工作区和 request/status 存储。
- 命中统计：`plans_hit=3848`、`plans_miss=135`、`buffers_hit=3561`、`buffers_miss=422`。
- `perf stat`：`Total Evolve 43.0325s→42.9994s`，仅 -0.08%；task-clock +0.1%，cycles -1.1%，instructions -0.69%，LLC miss rate 58.13%→58.50%。
- malloc self 1.93%→1.78%、cfree 1.73%→1.63%，`Parallel::transfer` total 27.98%→17.48%，但 `MPI_Waitall`、progress 和总墙钟没有同步下降。
- `perf record` 第二次运行反向（42.9962s→43.1940s）。
- 结论：拒绝并回退。

详细实现、原始计数和接受/拒绝理由见 [OPTIMIZATION_LOG.md](OPTIMIZATION_LOG.md)；探索过程见 `tuningrecord`。tuningrecord 只保留理由和结果，不作为实现规范。

### #18：kodis/lopsided 批处理（正式 40 步 PASS，已接受）

- 正式任务：115371（B=2，40 步，FINAL: PASS）。
- 短测任务：115236（首次 B=2 失败，z-parity 符号错）；115311（修复后 B=2/4/8 全部通过）。
- 候选思路：compute_rhs 末尾 24 lopsided + 24 kodis 每次都整块扫描、各自建 halo、每点判模板阶数；24 次共享同一组 shift，阶数只与 (i,j,k) 有关。把 B 个场打包进一次扫描，摊销 halo setup 和 per-point guard。
- 正式结果（job 115371）：This Program Cost=688.18s（精确 688.179s），Total Evolve=397.92s，Total Running=401.03s，40 步全部完成。
- 正式正确性：FINAL: PASS；Trajectory RMS=0.000000%，matched 40/40 golden timesteps（完整，非 prefix）；Constraint maxima level 0：Ham=0.27739667，Px=0.028132512，Py=0.031488238，Pz=0.026503396，全部 PASS（≤2.0）。
- 对比基线（#13，698.59s）：端到端 -10.4s（-1.5%），Total Evolve -13.8s（-3.4%）。
- `perf stat` 短测（baseline #13 vs B）：B=2 Total Evolve 44.13→41.15s（-6.8%），cycles -7.3%，instructions -3.3%，LLC miss rate 56.51%→54.66%。B=4/8 收益回吐（B=8 近基线），batch 过宽工作集挤出 L1/TLB，与 #10 POINTWISE 同一陷阱。
- 短测正确性：B=2/4/8 均 `FINAL: PASS`，RMS=0，constraint/bssn_BH 与基线逐位一致（snapshot 仅行首时间戳注释差异，剥除后 md5 全匹配；STATE=Failed 为 `set -e`+`pipefail` 下 `diff|head` 的脚本假象）。
- 结论：正式 40 步通过，#18（B=2）接受为新的当前有效最优，取代 #13。

### #19：RHS 运算顺序重排（perf 中性，拒绝）

- 任务：115428。
- 候选思路：compute_rhs_bssn_fused_ 循环体约 70 个标量同时 live，假设随意运算顺序导致栈溢出。做纯 lifetime-scheduling 重排——Phase 1 把浅 gauge RHS（Lap/beta）提到 chi_rhs/gij_rhs 后立即 store，Phase 4 把 Gamma_rhs/dtSf_rhs 提到组装完成后立即 store，让专属临时量尽早死亡。不改公式、不动模块内部 FMA 顺序。开关 `AMSS_RHS_REORDER`（默认 OFF）。
- 正确性：`FINAL: PASS`，RMS=0；constraint.dat 与 bssn_BH.dat 相对 baseline **逐位 IDENTICAL**。重排数值安全，依赖分析（Gamxa_s 跨 Ricci、Christoffel conformal/physical 双态、BAM 原子块）得到验证。
- `perf stat`：Total Evolve 42.56→42.72s（+0.38%），task-clock +0.41%，cycles +0.26%，instructions -0.07%，LLC miss rate 56.75%→56.56%（-0.19pp）。`perf record` 42.69→43.15s（+1.1%），方向一致（略升），属噪声。`compute_rhs_bssn_fused_` self% 18.58%→18.03%（分母效应）。
- 结论：bit-for-bit 正确但 perf 持平（全在 ±0.4%，<1% 门槛），**拒绝**。根因：gfortran -O3 指令调度阶段已做良好寄存器分配，源码语句重排无额外信息；能早死的浅 RHS 标量少，gup*/Christoffel/g*_t 真正长 live 节点跨 Ricci 不可早死，压不动峰值。代码保留在 `AMSS_RHS_REORDER` 宏后作为实验档案。

## 4. 当前 perf 事实

### 4.1 主热点

- `compute_rhs_bssn_fused_` 仍是最大的应用计算内核（fused 循环约200s/40步，占演化过半）；TwoPuncture 已不再是当前主热点。
- `kodis`、`lopsided` 已由 #18 批处理和 #26 尾部融合压缩；#27–#29 已证明简单 SIMD、字段 coloring、普通 j×k 重排没有稳定收益。
- `polint_` 已由 #42 的 `global_interp` uniform6 专用核压低；`polin3_uniform6_` 只占约0.14% self，剩余主要热点回到 `compute_rhs_bssn_fused_`、batch stencil 和 MPI 层级等待。
- `libopen-pal` 的 `mca_btl_sm_poll_handle_frag`/`opal_progress` 是 MPI 共享内存轮询，不应直接等同于 BSSN 数组 LLC miss。

### 4.2 通信路径

- **MPI 通信延迟只占演化的 ~4.4%**（lev8 ghost exchange ~14s）。perf % 里 42% opal_progress 是忙等 CPU 采样占比，不是墙钟占比。
- per-rank 负载插桩（MPI_Reduce min/max/sum wall per Step）确认：lev8 的 30 rank 各 1 block，但 block 大小差 60%（两个 BH moving patch 不对称），max/min wall = 1.41，浪费 ~34s（8.7%）。
- **AMR 负载不均衡是真实的但不能通过 rank 重分配修复**：greedy 按大小重分配破坏 patch 局部性，ghost exchange 从 on-rank memcpy 变成跨 rank MPI，tfer 翻倍，总时间 -24%。
- Berger-Oliger 层级串行中"粗层等细层"的部分（~60s）有约一半是 lev8 这个负载不均衡导致的——是计算不均衡伪装成通信等待。
- #14–#17、#19、#21、#27–#29 已实测证明：降低小 collective、缩小活动 rank、减少插值全数组归约、缓存传输对象、RHS 语句重排、greedy rank 重分配、简单 SIMD、cache coloring、普通 tile 重排，都没有形成稳定端到端收益。不要仅因为调用树某一行下降就接受优化。

### 4.3 演化时间构成（#18 reference breakdown；#35 正式 Evolve=397.011s）

| 组成 | 时间 | 占比 | 性质 |
|---|---:|---:|---|
| compute_rhs（lev8+lev7 为主） | ~200s | 50% | 计算 |
| Berger-Oliger 层级串行等待 | ~60s | 15% | 算法（含 lev8 负载不均衡） |
| MPI_Waitall 纯通信等待 | ~14s | 4% | 通信延迟 |
| data_packer pack+unpack | ~21s | 5% | 内存搬运 |
| analysis（插值，lev0 only） | ~54s | 14% | #4 已优化 |
| RestrictProlong + 其他 | ~49s | 12% | AMR 层间 |

## 5. 必须遵循的实验准则

### 5.1 Perf 是接受优化的硬条件

1. 静态分析只能提出假设，不能单独决定保留代码。
2. 每个候选至少做同节点 paired baseline/candidate 的 4 步短测。
3. `perf stat` 是墙钟和硬件计数器的主要判据；`perf record` 用于 self/total/DSO 调用树和热点归因，不把单次 record 墙钟当作稳定成绩。
4. 候选必须同时满足：
   - `FINAL: PASS`、Trajectory RMS=0；
   - paired `perf stat` 的目标指标确实改善；
   - `Total Evolve` 稳定下降，通常至少 1%；
   - cycles/instructions/cache miss 的恶化必须有明确解释，不能用一个局部热点下降掩盖总时间上升。
5. 如果 stat 与 record 方向相反，或总时间差小于约 1%，必须视为噪声/未证明收益；除非重复实验仍然一致，否则拒绝。
6. 只有通过短测 perf gate 的候选才提交 40 步正式任务。#14–#17 均未达到该条件，所以没有正式 40 步任务。

### 5.2 正确性顺序

- 短测：4 步、MPI=30、OMP=1、`--twop-cache`，运行 `check.sh`。
- 正式：40 步、MPI/OMP 与基线一致、TwoPuncture cache 关闭，记录完整 `Total Evolve` 和 `This Program Cost`。
- 修改通信顺序、collective 或数组布局后，必须检查 trajectory、constraint 和每个输出时间组；不能只看程序是否正常退出。
- 失败候选不得进入下一次累计基线；下一项必须从最近一次被接受的版本开始。

### 5.3 Perf 文件和报告命名

- 原始文件必须同时包含优化编号和代码版本，例如：
  `opt17_baseline_opt13.stat.txt`、`opt17_candidate_comm_cache.data`。
- 本地 canonical 报告只保留：
  - `{name}_self.txt`
  - `{name}_total.txt`
  - `{name}_dso.txt`
- 不再生成 `comm` 报告。
- 集群只编译、运行、采集 `perf stat` 和 `perf record` 原始数据；`perf report` 必须在本地执行。
- 计算节点 `perf record` 使用 `-F 99 -m 4 --call-graph fp`；`-m 4` 受节点 516KB mlock 限制，是当前 30 rank 可用配置。
- 版本映射和 #14–#17 报告清单见 [canonical README](perf_profiles/canonical_20260818_a/README.md)。

### 5.4 DevPod 只做编译检查

- DevPod 是共享 64 核机器，30 rank 4 步需 ~13min/步（集群独占节点 ~1s/步），4 步也跑不完，时间不可用且噪声极大。
- **不在 DevPod 上跑性能或正确性测试。** 正确性和 perf 一律上集群 paired 跑。
- DevPod 的唯一用途：验证代码能编译通过。

### 5.5 文档准则

- `HANDOVER.md`：只维护当前有效基线、实验结论、原始数据路径、下一步和不可违反的准则；失败方向要明确写“拒绝/回退”。
- `OPTIMIZATION_LOG.md`：详细记录假设、调用路径、编译开关、测试实现、集群任务、原始指标、报告路径和接受依据。
- `tuningrecord`：使用既定参考文风和第一人称，只写“为什么尝试、perf/计时结果、接受或回退”；不要写“用户告诉我”，不要堆砌类名、循环细节或通信 buffer 实现。不要写详细 perf 数据表或具体硬件计数器数值（perf stat 表、self% 明细、counter delta 等）——这些放进 OPTIMIZATION_LOG.md；tuningrecord 只留尝试的理由和时间收益结论（接受/回退加粗略时间差）。
- 不要把未执行的理论分析写成已验证结论；理论推断必须标明“未实测”。

## 6. 集群与本地工作流

### 6.1 集群提交

集群使用 `hpc`，不是本地 `sbatch`。提交前必须位于仓库目录 `/home/h3250104945/HPC101/src/lab4`，因为任务按提交时工作目录运行：

```bash
hpc submit -p lab4 -c 60 -m 100Gi -t 30m -d \
  -n <明确的任务名> -o optimization_logs/<明确的输出名> <script.sh>
```

当前 `lab4` 限制为最多一个活动任务；不要同时提交多个正式任务。使用持久输出文件和 `hpc info`/`hpc ls` 做必要检查，不要高频轮询。#14–#17 的有效任务号是 114175、114272、114350、114423；114160 和 114266 是脚本入口错误导致的失败提交，不包含实验数据。

### 6.2 本地报告

数据取回后，在本地执行：

```bash
./generate_perf_reports_local_20260818_a.sh <raw-label>
```

报告脚本会生成 self、total、DSO 三份文件；不要在集群上运行 `perf report`，不要把 stderr 中的 unwind 噪声当作应用结果。

## 7. 下一步建议

#42 已通过正式40步计时（job 124023，This Program Cost 342.549s，Evolve 329.018s，FINAL: PASS），是当前有效最优。TwoPuncture solver-only paired 由约277s降至8.47s；#27–#29、THP slab 和 hybrid MPI+OpenMP 仍然拒绝。

当前主要剩余热点是 `compute_rhs_bssn_fused_`、batch stencil 和 AMR 层级等待。后续候选仍需4步 paired perf/correctness gate，之后再做正式40步验证；不要重复已拒绝的 rank/SMT、THP 或普通 tile 方案。

## 8. 代码结构与已知约束

```text
AMSS_NCKU_Input.py        输入参数（正式脚本临时改写 Final_Evolution_Time）
AMSS_NCKU_Program.py      生成参数、运行 TwoPuncture/ABE、绘图
compile.sh/run.sh/check.sh
CMakeLists.txt            编译选项
src/bssn_class.C          RecursiveStep、Step、分析和层间控制
src/bssn_rhs.f90/.h       BSSN RHS 与 fused/POINTWISE 路径
src/MPatch.C               Patch、分布式插值
src/Parallel.C/.h          Sync、transfer、Restrict/Prolong、计时
src/cgh.C/.h               AMR Patch/Block 组织和 regrid
perf_profiles/             原始 perf 数据与本地 canonical 报告
optimization_logs/         集群运行日志
HANDOVER.md                本文件
OPTIMIZATION_LOG.md        详细实验主日志
tuningrecord               精简探索记录
```

- 网格层次：`cgh(levels=9) → Patch → Block → fgfs[i][idx]`。
- `fgfs` 是 SoA 的 167 个 double 数组；线性索引为 `i + j*nx + k*nx*ny`。
- `ghost_width=3`、`buffer_width=6`；AMR 层间数据经过 Sync、Restrict/Prolong、OutBd 等路径。
- 3D stencil 的 x/y/z 访问步长不同，没有一个简单循环顺序能同时让三个方向都连续向量化；任何 stencil 重排都必须先看汇编和 perf，再做数值验证。
- 必须保持 `AMSS_NCKU_Input.py` 的物理参数不被实验脚本永久改写；MPI/OMP/GPU 参数和临时演化时间要在日志中明确记录。

## 9. 仍未开始的任务

- GPU 任务尚未开始；`ABEGPU`、CUDA-aware MPI 和 device-resident 通信不能与当前 CPU #26 结果混淆。
- TwoPuncture 初值求解已完成通用 solver 优化并通过正式40步：当前配置为 `AMSS_TWOP_THREADS=30`、`AMSS_TWOP_NRELAX=10`，solver-only 约8.5s。后续若改变 TwoPuncture，仍需单独 solver-only paired gate，不能与 CPU RHS 候选混测。
