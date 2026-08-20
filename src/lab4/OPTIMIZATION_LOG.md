# AMSS-NCKU 优化日志（集群实验与正式计时）

> 本文件记录集群计算节点上的实验结果。记录分为 `探索`、`诊断` 和 `正式计时` 三种状态；只有明确标为 `正式计时` 的记录用于最终评分对比。
> 正式计时默认使用 MPI=30、Final_Evolution_Time=40.0、完整链路（无 `twop-cache`）；探索/诊断记录会明确标出临时参数，不能与正式成绩混用。
> 每次优化只追加一条记录，并保留作业 ID、原始运行结果、正确性结果和与前一基线的比较。
> DevPod 调试记录见 [DEVPOD_LOG.md](DEVPOD_LOG.md)。

## 工作范式(两阶段)

每次优化遵循以下两阶段流程:

### 阶段一:DevPod 调试(验证程序无误,不计时)→ 记入 DEVPOD_LOG.md
- 演化:`Final_Evolution_Time = 4.0`(临时,4 步快速验证)
- MPI:`MPI_processes = 30`(与生产一致)
- 用 `./run.sh --twop-cache`(初值不变时跳过 TwoPuncture)
- 验证:能编译 + 能跑完 + `./check.sh` 返回 **PASS**。短调试运行允许只匹配 golden 的前缀，但必须记录 `matched N/M`；只有 `N=M` 才能作为正式完整正确性证明。
- ⚠️ 若本次优化**改了 TwoPuncture 相关代码**,不能用 cache(初值会变)
- 详见 [DEVPOD_LOG.md](DEVPOD_LOG.md)

### 阶段二:集群正式计时(拿正式时间)→ 记入本文件
- 改回正式参数:`MPI_processes = 30`,`Final_Evolution_Time = 40.0`
- **不用 cache**(完整链路含 TwoPuncture)
- 提交命令(固定模板):
  ```bash
  cd /home/h3250104945/HPC101/src/lab4
  hpc submit -p lab4 -c 60 -t 30m -d -n <优化名> --chdir /home/h3250104945/HPC101/src/lab4 \
    bash -c 'bash compile.sh && bash run.sh && echo "=== CHECK ===" && bash check.sh'
  ```
- 跑完读 `This Program Cost` 和 check 结果,记入下方日志

> 提交坑:`hpc submit` 用提交时 cwd,务必先 `cd` 到 lab4 或用 `--chdir`。

### 步数不够时的外推准则
30min walltime 可能跑不完 40 步(~33 分钟)。当出现 Timeout/步数不足时:
1. **同时记录原始值和外推值**,缺一不可:
   - 原始值:实际跑到的步数 N、实际每步时间(列出或给平均)、实际状态(Timeout/Succeeded)
   - 外推值:用 N 步的实测平均单步时间 × 40,加固定开销,得到 40 步外推总时间
2. 外推前提:**同一次运行内每步演化时间稳定**(本程序已验证:34 步内单步 43.19±0.3s)。若优化后单步时间波动大,不可外推,必须想办法跑完(如降 rank 提速、或分两次跑)。
3. 固定开销(TwoPuncture + 绘图 + setup)经验值约 **287–297s**。后续完整运行（#3–#8）可直接拆出 **TwoPuncture 约 286s**、绘图约 1s，其余为 setup/运行波动。早期依据 perf4 CPU 样本占比推得的“TwoPuncture 176s + 绘图 110s”不能等同于墙钟拆分，已废弃。若优化 TwoPuncture、切换环境或冷缓存行为改变，必须重测，不能沿用该常数。
4. **评分/对比以外推值为准**(因为正式评测就是 40 步),但原始值保留以证明外推可靠。

## 集群运行记录模板（每次优化复制使用）

```markdown
## 优化 #N：<简短标题>（YYYY-MM-DD）
- 记录状态：`探索` / `诊断` / `正式计时`
- 优化内容：<针对什么瓶颈，改了哪些文件/代码，做了什么>
- 运行参数：MPI=?，OMP=?，Final_Evolution_Time=?，是否 `twop-cache`，编译器/编译选项，绑核/NUMA
- 作业与产物：job ID、脚本和关键报告路径
- 运行结果：This Program Cost = ?；Total Evolve Time = ?；实际步数/状态
- 正确性：`check.sh`、RMS/constraint、是否完整 PASS
- 对比基线：时间变化 ?（?%），正确性是否保持
- profiling/分析：<优化前后热点变化>
- 结论：<是否接受该方向，是否计入正式成绩>
- 备注：<失败/异常/外推依据/限制>
- 下一步：<后续方向>
- DevPod 调试记录：见 `DEVPOD_LOG.md` 对应条目
```

---

## Baseline（参考点）

## 优化 #0：Baseline（2026-08-16）
- 记录状态：`基准/外推`（正式 40 步未完成）
- 优化内容：无；原始代码 `-O3 -fno-strict-aliasing`
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0；编译器=g++/gfortran 14.2（OpenMPI 5）；编译选项=-O3；绑核=默认（未设置）
- 作业与产物：job 107261（正式链路，超时）；baseline20 job 102714（20 步正确性参考）；perf_stat job 103849；perf4 job 104167
- 运行结果：
  - **原始值(job 107261, Timeout)**: 实际跑到 **step 34/40**,状态 Timeout(30min walltime 用尽)
  - 实测每步时间(step 1–34, 单位秒):
    43.16, 43.23, 43.42, 43.38, 43.58, 43.59, 43.78, 43.83, 43.66, 43.46,
    43.13, 42.96, 43.05, 42.87, 43.03, 43.33, 43.35, 43.24, 43.01, 43.00,
    43.36, 42.92, 43.15, 42.88, 42.90, 43.06, 43.06, 43.32, 43.60, 43.24,
    43.06, 42.65, 42.60, 42.91
  - 平均单步 = **43.19s**(σ ≈ 0.30s,极稳定 → 外推可靠)
  - **外推值(40 步)**: 40 × 43.19 = **1727.6s** 演化
    + 固定开销初始估计 ≈ 287s（该次记录中的 perf4 占比拆分见“备注”，后续已修正）
    = **外推总时间 ≈ 2015s ≈ 33.6 分钟**(与 Timeout 吻合)
- 正确性：正式 job 未跑完 40 步，无 check 结果；baseline20（job 102714，20 步）的 RMS=0.000000%，constraint max(Ham=0.269，Px=0.028，Py=0.021，Pz=0.022)，**PASS**（证明代码正确，但正式链路 walltime 不足）
- 对比基线：—（本条即基准）
- profiling/分析（perf_stat job 103849，4 步，MPI=30）：
  - IPC=1.93, LLC-miss=39.41%, dTLB-miss=6.87%, page-faults 4.4M, CPUs utilized=11.7(⚠️ 全程平均,含 ~287s 单进程绘图;演化阶段实际满核,详见 DEVPOD_LOG 采样#2)
  - perf4 函数热点(job 104167, 4步 MPI=1, 81638样本):
    compute_rhs_bssn 16.7%, polint_ 9.9%, LineRelax_be 9.6%, __cos 8.4%, kodis_ 8.0%, LineRelax_al 6.7%, fdderivs_ 5.8%, Thomas 5.1%, lopsided_ 4.2%, malloc/cfree ~7%
  - 模块: BSSN演化~40%, TwoPuncture~23%, 内存管理~11%, 分析~10%, 数学库~10%
  - 负载不均衡: level0 只用 9/30 rank(40点切 min_width=12 只能切9块)
- 结论：作为后续正式计时的基准；40 步正式运行需要依靠优化压缩演化时间。
- 备注：
  - ⚠️ perf4 是 MPI=1,演化单步 ~131s(比 30 rank 的 43s 慢 3 倍,因为没并行);函数占比在"计算部分"可信,但"通信/等待占比"看不到。30 rank 的 MPI 等待占比已由 **perf5(job 107382, 30 rank)确认:~76% CPU 周期在 MPI 库忙等轮询**,即每步 ~33s 在 opal_progress(见 DEVPOD_LOG 采样#4)。
  - **固定开销拆分的版本说明**：本次早期依据 perf4 CPU 样本占比的估计为约 287s = **TwoPuncture ~176s**（单进程串行）+ **绘图 ~110s** + setup ~2s；这是历史估计，不能直接等同于墙钟拆分。后续 #3–#8 的完整运行直接测得固定开销约 287–297s，其中 TwoPuncture 约 286s、绘图约 1s。曾误写“绘图 277 + TwoPuncture 10”，已废弃。
  - 正式评测若也是 30min walltime,40 步外推 ~2015s 会 Timeout → **这是优化必须解决的核心矛盾**:要么提速到 30min 内,要么评测可能有不同的步数/walltime 配置(待评测说明确认)。
- 下一步：优先降低演化阶段时间，并保留完整 40 步正确性验证。
- DevPod 调试记录：无（baseline 不需要调试）

---

## 优化 #1：MPI rank 数探索（2026-08-17）

- 记录状态：`探索`（4 步 + `twop-cache`，非正式 40 步计时）；正式 40 步计时在选定 rank 后补。

- 优化内容：固定资源（60 核节点）下，探索 MPI rank 数对单步演化时间的影响。
- 运行参数：OMP=1，Final_Evolution_Time=4.0（临时），4 步 + `twop-cache`，无绑核。
- 作业与产物：job 107520，脚本 `rank_sweep.sh`。
- 运行结果（单步演化时间，`Total Evolve Time ÷ 4`）：
  | MPI rank | 模式 | 单步(s) | 每步稳定度 | 备注 |
  |---|---|---|---|---|
  | 9 | 物理核 | **56.7** | 55.8-57.8 | 绘图 287s(matplotlib 首启) |
  | 15 | 物理核 | **48.3** | 47.6-49.0 | 绘图 1s(缓存复用) |
  | 30 | 物理核 | **43.1** | 42.9-43.4 | 绘图 1s(缓存复用) |
  | 40 | `--use-hwthread-cpus` | **58.1** | 57.5-58.5 | 超线程争抢,更慢 |
  | 60 | `--use-hwthread-cpus` | **56.4** | 55.7-57.0 | 超线程争抢,更慢 |
  | 36/40/45/60 | 物理核 | ✗ 失败 | — | **cpuset=0-59(60 逻辑核)= 30 物理核,OpenMPI 默认每物理核 1 slot → 上限 30** |
- 正确性：本次为性能探索，未执行正式完整 `check.sh`；基线正确性由 baseline20（job 102714）提供。
- 对比基线：MPI=30 的 4 步结果为 43.1s/步，作为本次 rank 探索的内部基线；未产生正式 40 步成绩。
- 结论：**方向一（rank 数）：30 rank 是平台上限且实测最优，已用满。** 9→30 单调加速；“减 rank 减忙等”假设被推翻：
  - 76% 忙等是 CPU 周期**比例**,不是绝对时间;30 rank 每步 = 计算 ~10.3s + 忙等 ~33s;9 rank = 计算 ~34.4s + 忙等 ~22s
  - 减 rank 后每 rank 计算量 ×3.3,计算时间上升 > 忙等下降 → 总时间变长
  - 超线程(--use-hwthread-cpus)让 60 逻辑核都可用,但 40/60 rank 共享物理核争抢执行单元,单步反而升到 56-58s
- profiling/分析：
  - **附加发现（修正此前认知）**：绘图 287s 不是渲染时间，而是 **matplotlib 首启一次性成本**（同 job 内第一次绘图 287s，后续 1s；各配置画的图完全相同）。真正渲染 35 张图 ≈1s。此前“绘图 ~110s/~277s 正常值”均受此影响（107382 的 1s 也非数据损坏，是缓存复用）；评测若为全新环境，首启 287s 计入评分，是公共税。
  - **结论②（方向五：通信瓶颈）**：rank 增大到 30，单步仍持续下降（通信/忙等开销 < 并行收益），通信在 30 rank 内不是主导瓶颈；忙等是同步/负载不均衡，不是带宽。
- 备注：本条只用于选择后续正式运行的 MPI rank，不纳入正式成绩。
- 下一步：OpenMP 方向（rank×OMP 组合：计算 ~10.3s/步可被 OpenMP 加速；30 rank×OMP>1 会超线程争抢，需测 15×2、10×3 等组合）；选定组合后进行 DevPod 正确性验证和正式 40 步计时。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 调试 #1（rank_sweep 脚本验证）。

---

## 优化 #2：OpenMP 探索（2026-08-17）

- 记录状态：`探索`（4 步 + `twop-cache`，非正式 40 步计时）。代码改动已 commit（1b9c77d，13 处 `!$omp parallel do`）；默认编译不启用 OpenMP（零影响）。

- 优化内容：给热点 Fortran stencil 循环（fderivs/fdderivs/fdd*/kodis/lopsided）加 OpenMP，并测试 rank×OMP 组合是否优于 30:1。
- 运行参数：集群 4 步 + `twop-cache`，`AMSS_ENABLE_OPENMP=ON`；DevPod OMP=1/4；未进行正式 40 步运行。
- 作业与产物：job 107588（集群 rank×OMP 组合测试）；DevPod 验证记录见 `DEVPOD_LOG.md` 调试 #2。
- 运行结果（集群单步时间）：
  | 组合 | 单步(s) | vs 30:1 |
  |---|---|---|
  | **30:1** | **43.6** | — |
  | 15:2 | 53.1 | +22% |
  | 10:3 | 58.3 | +34% |
  | 6:5 | 51.7 | +19% |
- 正确性：DevPod OMP=1/OMP=4 均 **PASS 且约束值逐位相同**（Ham=0.22588831 等），并行不改变数值。
- 对比基线：30:1 = 43.6s/步，相对本次组合测试基线；与无 OpenMP 编译的 43.1s/步接近；未产生正式 40 步成绩。
- 结论：**方向二：加 OpenMP 不值得。** 网格块太小（40×40×20），13 处 parallel do 的线程同步/调度开销 > 并行收益；计算仅占单步 24%（10.3s/43.1s），收益上限 ~18% 被开销吃掉。DevPod 上 OMP=4 慢 75% 与此一致。
- profiling/分析：OpenMP 编译 + OMP=1（43.6s）≈ 无 OpenMP 编译（43.1s）；`!$omp` 指令在未启用时零开销，因此代码保留但正式计时不启用。
- 备注：本条为负向探索，不纳入正式成绩。
- 下一步：方向三/四（绑核 + NUMA）——比较默认无绑核与 `--bind-to core` 等配置。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 调试 #2。

---

## 优化 #3：球面插值 OpenMP（rank 内并行）（2026-08-17）

> 根因诊断见 [DEVPOD_LOG.md](DEVPOD_LOG.md) 调试 #3(transfer Waitall 计时 + Step 级 ana 分解,commit 47989cd/7d1f777):每粗步 32.9s 里 AnalysisStuff 占 32.89s,transfer 仅 0.12s → 瓶颈是分析阶段的单 rank 串行球面插值,29 rank 在 Allreduce 忙等。

- 记录状态：`正式计时`（附带 4 步 + `twop-cache` 快速迭代）。
- 优化内容：Interp_Points CPU 点循环加 `#pragma omp parallel for`（commit 9d37530，每点独立；llb/uub/varl 循环内私有），并按 3D 位置排序插值点以改善缓存局部性（commit 7512252）。
- 运行参数：正式运行 MPI=30，OMP=4，Final_Evolution_Time=40.0，完整链路（无 cache），`AMSS_ENABLE_OPENMP=ON`，无绑核。
- 作业与产物：快速迭代 jobs 108578/108636（4 步 + cache）；正式计时 job 108721；诊断 job 108677。
- 快速迭代结果（4 步 + `twop-cache`）：
  | 组合 | 单步(s) | ana | 说明 |
  |---|---|---|---|
  | 30:1 | 43.0 | 32.9s | 基线 |
  | 30:2 | 36.5 | 23.9s | 插值并行,线程被忙等 rank 挤占,仅 1.38× |
  | 30:4 | 35.4 | 23.9s | 同上(线程没核可用,4 线程≈2 线程) |
- 运行结果（job 108721）：**This Program Cost = 1650.9s ≈ 27.5 分钟，40 步全部跑完（≤30min walltime）**。
  - Total Evolve Time = 1363.0s,平均 34.1s/步
  - 固定开销 = 287.9s(TwoPuncture 真算 ~286s + 绘图 1s + setup 0.4s)
- 正确性：**FINAL: PASS**（40 步完整验证）。
- 对比基线：单步 43.19 → 34.1s（**-21%**）；总时间 2015 → 1650.9s（**-18%**，从超时到按时完成）。
- profiling/分析：`ana=23.9s` 仍占单步约 70%；插值线程被同核上忙等的 rank 挤占（30:2≈30:4），`OMPI_MCA_mpi_yield_when_idle=1` 测试无效（job 108677）。
- 结论：rank 内 OpenMP 能够提速，但只有一个 rank 实际执行球面插值，受其他 rank 忙等限制；要继续降低 `ana`，需要把球面插值分散到多个 rank。
- 备注：正式运行已完整跑完 40 步，结果可用于正式比较。
- 下一步：实现球面插值跨 rank 分散（见优化 #4）。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 调试 #3。

---

## 优化 #4：球面插值分布式分散（跨 rank 并行）（2026-08-17）

- 记录状态：`正式计时`。
- 优化内容：Interp_Points CPU 路径（`#else` 分支，NN≥256 时）改为：① 找含球壳的中心块（owner=rank 4）；② owner 打包 fgfs+X，通过 `MPI_Bcast` 广播到全 30 rank；③ 每 rank 只插值 NN/30 个点切片（去掉 `myrank==BP->rank` 门控）；④ Allreduce 不变。NN<256（单点 BH 探针）走原路径。
- 原理：原来只有中心块 owner（rank 4）串行插值 36864 点，29 个 rank 干等 Allreduce；现在 30 个 rank 各算约 1229 点（每点独立，Allreduce SUM 还原完整结果，数值等价）。数据所有权机制（fgfs 只在 owner 分配）决定必须先广播中心块数据到全 rank。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，完整链路（无 cache），`AMSS_ENABLE_OPENMP=ON`（开关开但 OMP=1），无绑核。
- 作业与产物：正式计时 job 110492；相关诊断/采样见 `DEVPOD_LOG.md` 和 perf6 记录。
- 运行结果（job 110492）：**This Program Cost = 764.2s ≈ 12.7 分钟，40 步全部跑完**。
  - Total Evolve Time = 467.1s,平均 11.7s/步(σ 小,极稳定)
  - lev0 STEP ana 分解:`wall=1.4s tfer=0.2s ana=1.4s` ← ana 从 #3 的 23.9s 降到 ~1.4s(**17×**)
  - 固定开销 = 297.1s(TwoPuncture ~286s + 绘图 ~1s + setup;含无 cache 真算)
- 正确性：**FINAL: PASS**（Trajectory RMS=0，Constraint max Ham=0.27739667，与 #3 逐位一致）。
- 对比基线：单步 43.19 → 11.7s（**-73%**）；总时间 2015 → 764.2s（**-62%**）。
- 对比 #3：单步 34.1 → 11.7s（**-66%**）；ana 23.9 → 1.4s。
- profiling/分析（perf6，30 rank、OMP=1、4 步 + cache、`-m 4`、141K 样本；采样记录见 [DEVPOD_LOG.md](DEVPOD_LOG.md) perf 采样 #5）：
  - DSO 归类:ABE 计算 53.92%(baseline <2%,大幅回升)、libopen-pal 27.64%(baseline ~76%,大幅下降但未消失)、libc 12.17%、libmpi 5.13%
  - 函数符号:compute_rhs_bssn 25.48%、polint 5.40%、kodis 4.80%、fdderivs 4.57%、lopsided 3.98%、prolong3 2.34%、memcpy/memset/malloc/cfree 合计 ~14.5%、opal_progress 0.32%
  - **修正此前"忙等基本消失"的判断**:忙等没消失,而是从 ~76% 降到 27.64%。**分析阶段(ana)的忙等消失了**(对应 ana 23.9→1.4s),但**计算阶段的同步等待仍在**——来自 RecursiveStep 各层 RestrictProlong→transfer 的集体同步、RK4 子步间 Sync、新加的 16 次 Bcast。30 rank 在通信点同步时快等慢,仍产生 opal_progress 轮询
  - perf6 证明:ABE 计算占比从 <2% 回升到 53.92%,说明忙等腾出的 CPU 给了有效计算(分散成功);但残留 27.64% libopen-pal 是下一瓶颈——计算阶段 AMR 层间通信的集体同步等待,要降它得动 transfer/RestrictProlong 的通信结构(方向 5)
- 结论：收益来自“把 1 个 rank 的串行工作摊到 30 个 rank”；分析阶段瓶颈显著下降，但计算阶段同步等待仍是主要剩余瓶颈。
- 备注：广播开销（中心块约 534KB/次）被摊到 30 rank 后可忽略。OpenMP（#3 的 rank 内并行）在分散后收益有限，延续 #2 的“块太小，OpenMP 无效”结论，因此 OMP 用 1。lev0 中心块 owner 稳定（static_grid_level，AMR 不 regrid 静态层）。
- 下一步：① 764s 中 TwoPuncture 286s 占 37%（单进程串行，LineRelax/Thomas 未并行）；② 计算阶段 MPI 同步等待 27.64%（perf6 libopen-pal），需要更深 profile 定位具体同步点。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 对应调试/采样条目。

---

## 优化 #5：-march=native -mtune=native 编译实验（无效）（2026-08-17）

- 记录状态：`正式计时`（负向结果）。
- 优化内容：在 #4 分布式分散基础上，将编译开关设为 `AMSS_ARCH_FLAGS="-march=native -mtune=native"`，尝试启用集群鲲鹏 920B 的 SVE。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，完整链路（无 cache），`AMSS_ENABLE_OPENMP=ON`，`AMSS_ARCH_FLAGS="-march=native -mtune=native"`。
- 作业与产物：正式计时 job 110876。
- 运行结果：This Program Cost = 776.2s（Total Evolve 未单独记录）。
- 正确性：**FINAL: PASS**；约束值 Ham=0.27739667，与 #4 逐位一致。
- 对比基线（#4，764.2s）：+12s（+1.6%），**反而略慢**。
- profiling/分析：双精度 SVE lane 数与 NEON 相同（均 2-wide、128 位/8 字节），无直接吞吐提升；SVE 谓词/gather 在双精度 AXPY 段收益微小，可能因指令序列变长导致 icache 压力。此前汇编分析显示 compute_rhs 已有约 1100 条 NEON 向量指令，编译器向量化收益空间有限；手写 intrinsic 未实测，以上部分为分析推断。
- 结论：放弃该编译开关；本次未获得性能收益。
- 备注：此前还分析了 SIMD intrinsic/OpenMP/block-patch 并行/level 并行四条方向，其中只有 #2 stencil OpenMP 和本次 native 编译实测，其余（SIMD 手写、RK4 OpenMP、block 负载重分配、level 并行）仅为理论分析，尚未验证。详见 tuningrecord。
- 下一步：转向 TwoPuncture 并行或 rank9–29 计算阶段工作重分配。
- DevPod 调试记录：未在本条中单列；正确性来自正式 job 的 FINAL PASS。

---

## 优化 #6：stencil 子程序 fh 缓冲区复用（2026-08-18）

> #4 分布式分散后,我盯着 perf6 的热点看,发现 malloc 加 cfree 占了 3.3%,想着这块能不能省点。仔细看调用栈,发现是 fderivs、fdderivs、kodis、lopsided 这几个 stencil 子程序在反复 malloc 一个临时数组 fh。

- 记录状态：`正式计时`。
- 优化内容：将 stencil 子程序里反复 malloc/free 的边界填充缓冲区 fh 改为 module 级 save 数组，只 allocate 一次后复用。
  - **问题根因**: fh 是 dimension(-1:ex(1),-1:ex(2),-1:ex(3)) 的 automatic array,gfortran 把大的 automatic array 放堆上,每次调用 fderivs/fdderivs 就 malloc 一次,用完 free。每次 compute_rhs 要调 fderivs 约 21 次加 fdderivs 约 11 次,每次各 malloc 一个 fh,lev8 时约 784KB。perf6 显示 malloc 加 cfree 占 3.3%
  - **fh 作用**: 边界镜像填充副本。stencil 计算在边界附近要读正负 2 个邻居格点,block 边界外没真实数据,用对称条件把边界外虚构出来。symmetry_bd 把输入数组 f 复制到 fh,然后把下标 0、-1 的位置按对称条件填上镜像值,这样 stencil 循环能统一访问 fh(i-2:i+2) 而不需要额外边界判断
  - **改动**: 新建 diff_buffers module,内含 fh_buf 全局一份,懒初始化,不够大就 realloc。每个子程序把 fh 改成 pointer,用 fh_ensure(ex) 保证够大,再 fh(-1:ex(1),...) => fh_buf 把 pointer 关联过去,带负下标
  - 改动文件: src/diff_new.f90 新增 diff_buffers module,fderivs/fdderivs 内 fh 改 pointer; src/kodiss.f90 / src/lopsidediff.f90 同样改动
- DevPod 调试记录：4 步、`twop-cache`，513.6s；与 #4 的 514s 接近。步数少，malloc 开销占比小；正确性 PASS。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，完整链路，无绑核。
- 作业与产物：正式计时 job 111015。
- 运行结果（job 111015）：**This Program Cost = 749.0s，40 步全部跑完**。
  - Total Evolve Time = 459.8s,平均单步 11.5s
  - 固定开销 ≈ 289s(TwoPuncture 286s + 绘图 + setup)
- 正确性：**FINAL: PASS**，约束值与 #4 逐位一致。
- 对比基线（#4，764.2s）：-15.2s（-2.0%）；malloc 加 cfree 从 3.3% 降到约 1.8%。
- profiling/分析：malloc 开销减半，但这 2% 也接近该方向上限，继续在 stencil 细节上抠收益递减。
- 结论：复用 fh 缓冲区减少了资源管理开销，保留该改动；不改变算法和数值结果。
- 备注：每次 compute_rhs 约省掉 32 次 malloc，在长时间演化中累积收益。fh 缓冲区复用是内存分配优化的典型案例。
- 下一步：764s 中 TwoPuncture 286s 占 37%，是最大单块；perf6 中 libopen-pal 残留 27.64% 同步等待也需要继续定位。
- DevPod 调试记录：本条已在上方列出；详细输出见 `DEVPOD_LOG.md` 对应条目。

---

## 优化 #7：sync 与 compute 重叠探索（实测否定）（2026-08-18）

> perf6 残留 27.64% libopen-pal 同步等待,我原本想通过"计算-通信重叠"消掉一部分——把 Sync 的 Waitall 和 compute 重叠。做了精细插桩 + 依赖链分析 + 跳过实验,结论是这个程序的重叠**没有收益**。

- 记录状态：`诊断/探索`（4 步，非正式计时；跳过 Sync 的对照运行故意产生数值错误，仅用于墙钟测量）。
- 优化内容：评估“Sync 的通信等待能否被 compute 隐藏”（计算-通信重叠）。
- 运行参数：4 步；正常对照与跳过 lev8 Sync 的配置均使用同一 MPI/OMP 环境；采用 `twop-cache`，不用于正式成绩。
- 作业与产物：插桩 job 111199；跳过 Sync 对照 job 111220；Parallel.C 中保留 pack/wait/unpack 计时插桩。
- 运行结果：正常对照 Total Evolve=47.07s；跳过 lev8 Sync Total Evolve=46.63s，差 0.44s。
- 正确性：正常对照链路保持正确；跳过 lev8 Sync **数值错误**，不作为正确性结果。
- 对比基线：跳过 lev8 Sync 相对正常对照仅快 0.44s，不能作为正式性能收益。
- 分析过程:
  1. **依赖链确认**: compute(f0^{k+1}) 依赖 Sync(f1^k) 完成(RK4 更新 f1 全数组含 ghost,Sync 修正 ghost,swap 后 f0=f1)→ 跨子步 pipeline 不可行,唯一出路是拆 compute_rhs 内部段+边界段
  2. **精细插桩**(Parallel.C pack/wait/unpack 分段计时,job 111199 4 步): wait 占 tfer 的 85%,pack/unpack 只 15%;粗层(lev4/5/6)wait 超过自身 wall 的 1.2-3 倍 → 粗层 wait 等的是深层 RestrictProlong 反馈(层级串行),不是通信
  3. **跳过实验**(job 111220): 临时跳过 lev8 的 Sync(数值错,只测墙钟),4 步 Total Evolve 47.07→46.63s,只省 0.44s → lev8 的 Sync 等待不在关键路径上
- 关键数据:
  | 配置 | 4 步 Total Evolve | lev8 YN=1 wall |
  |---|---|---|
  | 正常(job 111199) | 47.07s | 165ms |
  | 跳过 lev8 Sync(job 111220) | 46.63s | 166ms |
- profiling/分析结论：**实测否定**。重叠无收益——通信本身已是异步（Isend/Irecv），Waitall 等数据时别的层还在算，等待被 Berger-Oliger 层级计算自然隐藏；关键路径是层级串行，不是通信。即使完美重叠，收益 ≈ 0。
- 备注：
  - 依赖链分析确认：compute 读 f0.ghost = 上轮 Sync 的结果；swap 不破坏 ghost；RK4 用 f0.ghost 复制 f1.ghost（错）再被 Sync 修正 → compute 必须等上轮 Sync。
  - 唯一理论出路（拆 compute_rhs 内部/边界，内部 57% 不读 ghost）因实测收益 ≈ 0 而放弃。
  - 插桩代码保留（Parallel.C pack_time/wait_time/unpack_time），未来分析可用。
- 结论：不接受该方向，不纳入正式成绩；层间并行/通信重叠需要先解决 Berger-Oliger 层级依赖。
- 下一步：关键路径是 Berger-Oliger 层级串行（RecursiveStep），不是通信；要提速需改层级并行（深水区）或 TwoPuncture（286s，最大单块）。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 对应插桩诊断条目。

---

## 优化 #8：compute_rhs 单点化重写（寄存器融合）（2026-08-18）

> compute_rhs_bssn 是最大热点(perf6 25.48% + 子调用 15%),43 个临时 3D 数组反复读写(LLC-miss 39%)。单点化把组装部分(134-817 行)重写为显式 do k,j,i 循环,循环内中间量用标量(寄存器),消除临时数组。

- 记录状态：`正式计时`。
- 优化内容：新增 compute_rhs_bssn_fused 函数（与原函数共存，宏切换）。
  - **结构**: 32 个 stencil 调用(fderivs×21 + fdderivs×11)保留在循环前(整数组);组装 175 条公式在 do k,j,i 循环内(标量中间量 + (i,j,k) 下标);24 lopsided + 24 kodis 在循环后
  - **两义变量按代拆分**: gxxx-gzzz 两代(fderivs 输出 vs 组装覆盖),fxx-fzz 被 8 个 fdderivs 复用 → 每代独立数组(bxxx/fxx_dxx 等 66 个新数组)
  - **验证器**: 175 条组装公式规范形式逐条比对(0 不匹配),修复 2 处转写错误(gupyz 误写 gupyy、96 处未声明变量)
- DevPod 调试记录：4 步，`twop-cache`；FINAL PASS，约束值 Ham=0.22588831 等与原版**逐位一致**。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，完整链路，`AMSS_ENABLE_OPENMP=ON`。
- 作业与产物：正式计时 job 111499。
- 运行结果（job 111499）：**This Program Cost = 721.7s，40 步全部跑完**。
  - Total Evolve Time = 433.7s,平均单步 10.84s(前 15 步实测 10.5-11.0s)
  - 固定开销 ≈ 288s(TwoPuncture 286s + 绘图 + setup)
- 正确性：**FINAL: PASS**；约束值 Ham=0.27739667，与 #6 逐位一致。
- 对比基线（#6，749.0s）：**-27.3s（-3.6%）**。
- profiling/分析：
  - 编译器自动融合(-floop-nest-optimize)实测无效(函数体 29485 行/1111 NEON fmla 完全一样)→ 必须源码级单点化
  - 独立基准: 单点循环 vs whole-array,求逆段实测快 36%(0.220s→0.140s)
  - 实际收益 -3.6% 小于预估(-10%),因为导数整数组(40个)仍是访存大头,单点化只消除了组装临时数组
- 结论：源码级单点化带来 3.6% 端到端收益；保留 fused 路径，后续继续围绕剩余访存热点优化。
- 备注：宏切换（bssn_rhs.h fortran3 分支 → compute_rhs_bssn_fused_）；回退原版改回 compute_rhs_bssn_ 即可。fused 函数保留在 bssn_rhs.f90 末尾。
- 下一步：① 分块（cache blocking）进一步压 lev8 访存；② TwoPuncture 并行（286s，端到端大块）；③ 目标端到端 340s（100 分），当前 721.7s。
- DevPod 调试记录：本条已在上方列出；详细输出见 `DEVPOD_LOG.md` 对应条目。

---

## 优化 #9：HALO 边界视图（2026-08-18）

- 记录状态：`4 步 perf 探索`；未作为正式 40 步成绩。
- 优化内容：在等价对称的 z 边界只 materialize 需要的两到三层 halo，内部点直接读取原数组；保留对称阶数大于一时的旧路径作为正确性回退。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=4.0，`--twop-cache`，`-O3 -g -fno-omit-frame-pointer`；集群只采集 `perf stat` 和 `perf record` 原始数据。
- 作业与产物：job 112350；原始数据为 `perf_profiles/perf9_halo_20260818_a.stat.txt`、`.data`；规范报告为 `perf_profiles/canonical_20260818_a/opt09_halo_20260818_a_{self,total,dso}.txt`。
- 正确性：**FINAL: PASS**；Trajectory RMS=0，约束值与 #8 控制逐位一致。
- 结果：相对 #8 LEGACY，LLC miss rate 从 64.84% 降至 63.75%，但 task-clock 从 1.488e6 ms 升至 1.517e6 ms，指令数从 9.242T 升至 9.523T；4 步演化约 47.60s 降至 47.00s，处于运行波动范围。
- 结论：halo 视图确实减少了少量 LLC 流量，但没有形成稳定的端到端收益；不把它单独作为性能结论，保留为 POINTWISE 的边界基础设施。
- 下一步：继续观察 POINTWISE 的逐点导数路径。

---

## 优化 #10：POINTWISE 逐点导数路径（2026-08-18）

- 记录状态：`4 步 perf 探索`；未作为正式 40 步成绩。
- 优化内容：对等价对称的 bulk 点直接展开一阶和二阶差分模板，只在反射 z 边界调用通用 helper；跳过完整块级导数数组的预处理。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=4.0，`--twop-cache`，同 #9 的编译和 perf 参数。
- 作业与产物：job 112363；原始数据为 `perf_profiles/perf10_pointwise_20260818_a.stat.txt`、`.data`；规范报告为 `perf_profiles/canonical_20260818_a/opt10_pointwise_20260818_a_{self,total,dso}.txt`。
- 正确性：**FINAL: PASS**；Trajectory RMS=0，约束值与 #8 控制逐位一致。
- 结果：相对 #9，LLC loads 下降 15.5%，LLC misses 下降 23.9%，miss rate 从 63.75% 降至 57.40%；但指令数上升 10.8%，L1D miss 上升 45.2%，dTLB miss 上升 57.1%，4 步演化变为 48.19s。
- 结论：POINTWISE 明显降低了 LLC 流量，却把太多数组同时带入 L1 和 TLB，净性能反而变差。保留该路径作为后续标量化的基线，不接受本版本作为独立性能优化。
- 下一步：先把每个点的模板阶数判断从 26 个导数组中提出来。

---

## 优化 #11：每点只判断一次模板阶数（2026-08-18）

- 记录状态：`4 步 perf 探索`；接受为 POINTWISE 的有效增量优化。
- 优化内容：在每个 bulk 点先计算一次 `pw_order`，取四阶、二阶或无模板；一阶和二阶内联模板只读取这个结果，不再为每个导数组重复比较 i、j、k 边界。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=4.0，`--twop-cache`，同 #10 的 perf 参数。
- 作业与产物：job 112589；原始数据为 `perf_profiles/perf11_pointwise_order2_20260818_a.stat.txt`、`.data`；规范报告为 `perf_profiles/canonical_20260818_a/opt11_pointwise_order_hoisted_20260818_a_{self,total,dso}.txt`。
- 正确性：DevPod 4 步 smoke 和集群 4 步输出均 **FINAL: PASS**；Trajectory RMS=0，约束值逐位一致。
- 结果：相对 #10，Total Evolve 48.19s 降至 46.20s，约快 4.1%；指令数 10.55T 降至 9.39T，下降 11.1%；cycles 下降 8.9%，dTLB miss 下降 5.8%，LLC miss 基本不变。
- 结论：重复的模板边界判断和随之产生的下标计算确实是实质开销；`pw_order` hoist 保留。
- 下一步：继续去掉 POINTWISE 二阶导结果的逐点三维数组写回。

---

## 优化 #12：POINTWISE 二阶导结果点内标量化（2026-08-18）

- 记录状态：`4 步 perf 探索`；接受为有效访存优化。
- 优化内容：POINTWISE 路径的 66 个二阶导结果改存点内标量，直接供移位项、Ricci 张量、chi 和 Lap 的协变二阶导计算使用；LEGACY/HALO 路径仍保留原有 fdderivs 三维数组。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=4.0，`--twop-cache`，同 #11 的 perf 参数。
- 作业与产物：job 112752；原始数据为 `perf_profiles/perf12_pointwise_d2scalar_20260818_a.stat.txt`、`.data`；规范报告为 `perf_profiles/canonical_20260818_a/opt12_pointwise_d2_scalar_20260818_a_{self,total,dso}.txt`。
- 正确性：DevPod 2 步和集群 4 步均 **FINAL: PASS**；Trajectory RMS=0，约束值与 #8 控制逐位一致。
- 结果：相对 #11，Total Evolve 46.20s 降至 43.82s，约快 5.1%；L1D miss 40.82B 降至 33.97B，下降 16.8%；dTLB miss 224.08B 降至 196.81B，下降 12.2%；LLC miss 仅下降 2.5%，基本没有变化。指令数升至 9.68T，但 cycles 基本不升。
- 结论：二阶导数组写回确实是 POINTWISE 工作集过宽的来源之一；标量化降低了 L1/TLB 压力，保留该版本。
- 下一步：测量并评估 45 个一阶导结果的点内标量化；同时把 `kodis`、`lopsided` 的重复整块扫描列为独立候选，不与当前 RHS 改动混测。

---

## 优化 #13：POINTWISE 一阶导结果点内标量化（2026-08-18）

- 记录状态：`正式 40 步运行`。
- 优化内容：POINTWISE 路径的 45 个一阶导结果改为当前网格点的局部变量；公共组装段通过点内别名直接消费这些变量，LEGACY/HALO 路径继续从原有三维数组读取。
- 运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，无 TwoPuncture cache，`AMSS_ENABLE_OPENMP=ON`，`-O3`。
- 作业与产物：job 113370；脚本 `run_pointwise_d1scalar_full_20260818_a.sh`；完整日志为 `optimization_logs/opt13_pointwise_d1scalar_full_20260818_a.run.log`。
- 运行结果：**This Program Cost=698.59s**，Total Evolve=411.75s，40 步全部完成。
- 正确性：**FINAL: PASS**；Trajectory RMS=0，约束值与 #8 基线逐位一致，level 0 Ham=0.27739667。
- 对比基线（#8，721.7s）：下降 23.1s，约 3.2%。本次是完整任务结果，不是 4 步缓存外推。
- 结论：一阶导中间数组写回也确实带来额外访存；标量化版本保留，下一步单独验证 `kodis/lopsided` 点内批处理。
- 下一步：暂停 `kodis/lopsided` 流水线融合，先用现有 perf 和通信计时定位 MPI 开销来源。

---

## 优化 #14：每个 Step 合并 ERROR Allreduce（perf 否定，2026-08-18）

- 记录状态：`4 步 perf 探索`；正确性通过，性能不接受，代码回退到 #13 基线。
- perf 假设：40 步正式运行约 2640 次 `Step`，当前 predictor 加三个 corrector 每次都对单个 ERROR 整数执行 world `MPI_Allreduce`，累计约 10560 次；将本地错误累积到 Step 末尾并只规约一次，应该减少 `opal_progress` 和 collective 同步。
- 实现开关：`AMSS_STEP_ERROR_REDUCE=ON`。候选在 `bssn_class::Step` 中延迟 predictor/corrector 的错误归约，Step 末尾使用一次 `MPI_MAX`；默认关闭时保留原路径。该候选不改变正常无错误的数值公式。
- 测试脚本：`profile_mpi_opt14_17_cluster_20260818_a.sh`，脚本在同一集群任务中先后运行 #13 基线和 #14 候选；集群只生成 stat/data/log，本地用 `generate_perf_reports_local_20260818_a.sh` 生成报告。
- 作业与产物：job 114175；原始文件为 `perf_profiles/opt14_baseline_opt13.{stat.txt,data,run.log}`、`perf_profiles/opt14_candidate_step_error.{stat.txt,data,run.log}`；报告为 `perf_profiles/canonical_20260818_a/opt14_{baseline_opt13,candidate_step_error}_{self,total,dso}.txt`。
- 正确性：两个 4 步结果均 `FINAL: PASS`，Trajectory RMS=0，候选约束值与基线一致。
- `perf stat` 配对结果：

  | 指标 | #13 baseline | #14 candidate | 变化 |
  |---|---:|---:|---:|
  | Total Evolve | 55.2737s | 57.4241s | +3.9% |
  | task-clock | 1,692,398.55ms | 1,766,255.55ms | +4.4% |
  | cycles | 4.671e12 | 4.820e12 | +3.2% |
  | instructions | 1.1097e13 | 1.1446e13 | +3.1% |
  | LLC miss rate | 57.06% | 55.04% | -2.02pp |

- `perf record` 的候选运行 Total Evolve 为 44.5974s，而基线为 53.4093s；该方向与 stat 相反，说明采样运行存在较大噪声，不能以单次 record 墙钟接受优化。
- total 调用树中候选不再出现 `PMPI_Allreduce`，说明改动确实减少了该 collective；但 paired stat 没有显示端到端收益。结论：**拒绝 #14**，下一项 #15 直接基于 #13，不带 `AMSS_STEP_ERROR_REDUCE`。

## 优化 #15：分布式插值改用 MPI_Allgatherv（perf 否定，2026-08-18）

- 记录状态：`4 步 perf 探索`；正确性通过，性能不接受，代码回退到 #13 基线。
- perf 假设：分布式 `Patch::Interp_Points` 已把点划成各 rank 的连续切片，但仍对完整零填充数组做 Allreduce；改为紧凑切片 Allgatherv 可以减少无效零数据和 weight 归约。
- 实现开关：`AMSS_INTERP_ALLGATHERV=ON`。候选只替换 `NN >= 256` 的 CPU 分布式插值聚合，保留字段/坐标 Bcast、排序和原始索引逆置换；其余路径仍使用 legacy Allreduce。
- 测试脚本：`profile_mpi_opt14_17_cluster_20260818_a.sh`；候选日志额外输出 Allgatherv 次数和累计 gather 时间。集群只生成 stat/data/log，本地生成 self/total/DSO 三份报告。
- 作业与产物：job 114272；原始文件为 `perf_profiles/opt15_baseline_opt13.{stat.txt,data,run.log}`、`perf_profiles/opt15_candidate_interp_allgatherv.{stat.txt,data,run.log}`；报告为 `perf_profiles/canonical_20260818_a/opt15_{baseline_opt13,candidate_interp_allgatherv}_{self,total,dso}.txt`。
- 正确性：两个版本均 `FINAL: PASS`，Trajectory RMS=0；候选在两个 perf 运行中各执行 32 次大规模插值聚合（`NN=36864`）。
- `perf stat` 配对结果：

  | 指标 | #13 baseline | #15 candidate | 变化 |
  |---|---:|---:|---:|
  | Total Evolve（stat 运行） | 42.6377s | 42.9712s | +0.8% |
  | task-clock | 1,384,688.00ms | 1,396,238.72ms | +0.8% |
  | cycles | 3.907e12 | 3.937e12 | +0.77% |
  | instructions | 9.716e12 | 9.723e12 | +0.07% |
  | LLC miss rate | 55.74% | 56.99% | +1.25pp |

- total 调用树中 `PMPI_Allreduce` 仅由 15.43% 降至 14.90%，`PMPI_Allgatherv` 占 0.49%；该 collective 不是总通信瓶颈。
- `perf record` 第二次运行的 Total Evolve 为 42.8949s → 42.7036s，差异很小且与 stat 方向相反；以 paired stat 为接受依据，**拒绝 #15**。下一项 #16 直接基于 #13。

## 优化 #16：每层活动 rank communicator（目标 MPI 无收益，2026-08-18）

- 记录状态：`4 步 perf 探索`；正确性通过，目标通信没有改善，性能结论不接受，代码回退到 #13 基线。
- perf 假设：粗层只由少数 rank 持有 Block，却仍让全 world 参与 ERROR Allreduce；为每层建立活动 rank communicator，既可缩小 collective，也可能允许无任务 rank 提前进入后续层。
- 实现开关：`AMSS_LEVEL_ACTIVE_COMM=ON`。初始化和 regrid 后根据 Block owner 位图创建/刷新 `MPI_Comm_split`；ERROR collective 只由该层活动 rank 调用，点对点 transfer 仍保留 world communicator。
- 测试脚本：`profile_mpi_opt14_17_cluster_20260818_a.sh`；候选输出每层活动 rank 位图，集群只生成 stat/data/log，本地生成 self/total/DSO 报告。
- 作业与产物：job 114350；原始文件为 `perf_profiles/opt16_baseline_opt13.{stat.txt,data,run.log}`、`perf_profiles/opt16_candidate_level_active_comm.{stat.txt,data,run.log}`；报告为 `perf_profiles/canonical_20260818_a/opt16_{baseline_opt13,candidate_level_active_comm}_{self,total,dso}.txt`。
- 正确性：两个版本均 `FINAL: PASS`，Trajectory RMS=0；候选 level 0 活动 rank=9，level 1–8 活动 rank=30。
- `perf stat` 配对结果：

  | 指标 | #13 baseline | #16 candidate | 变化 |
  |---|---:|---:|---:|
  | Total Evolve（stat 运行） | 45.2985s | 44.1994s | -2.4% |
  | task-clock | 1,490,261.94ms | 1,457,313.11ms | -2.2% |
  | cycles | 4.177e12 | 4.091e12 | -2.1% |
  | instructions | 1.0223e13 | 1.0126e13 | -0.95% |
  | LLC miss rate | 57.07% | 57.09% | +0.02pp |

- total 调用树中 `PMPI_Allreduce` 17.74% → 17.75%，`PMPI_Waitall` 19.34% → 18.61%，`opal_progress` 42.51% → 43.75%；没有目标 collective 的稳定下降。
- `perf record` 第二次运行 Total Evolve 44.5332s → 44.7841s，方向与 stat 相反；单次 stat 的下降不能作为接受依据。结论：**拒绝 #16**，下一项 #17 直接基于 #13。

## 优化 #17：通信计划与工作缓冲区缓存（perf 无稳定收益，2026-08-18）

- 记录状态：`4 步 perf 探索`；正确性通过，目标分配/传输路径有所改善但端到端没有稳定收益，代码回退到 #13 基线。
- perf 假设：现有通信计划只缓存 segment 几何，`Parallel::transfer` 每次仍重新分配 send/receive 数组、request/status 数组并重新探测长度；缓存这些运行期对象可以降低 malloc/cfree 和 pack 前的管理开销。
- 实现开关：`AMSS_CACHE_COMM_BUFFERS=ON`（自动启用 `AMSS_CACHE_COMM_PLANS`）。按 `CommPlan + VarList source/target + symmetry` 缓存 peer 长度、vector 工作缓冲区和 request/status 存储；regrid 时由 `clear_plan_cache()` 一并失效。程序结束打印命中统计。
- 测试脚本：`profile_mpi_opt14_17_cluster_20260818_a.sh`；集群只生成 stat/data/log，本地生成 self/total/DSO 报告。
- 作业与产物：job 114423；原始文件为 `perf_profiles/opt17_baseline_opt13.{stat.txt,data,run.log}`、`perf_profiles/opt17_candidate_comm_cache.{stat.txt,data,run.log}`；报告为 `perf_profiles/canonical_20260818_a/opt17_{baseline_opt13,candidate_comm_cache}_{self,total,dso}.txt`。
- 正确性：两个版本均 `FINAL: PASS`，Trajectory RMS=0；候选统计 `plans_hit=3848`、`plans_miss=135`、`buffers_hit=3561`、`buffers_miss=422`。
- `perf stat` 配对结果：

  | 指标 | #13 baseline | #17 candidate | 变化 |
  |---|---:|---:|---:|
  | Total Evolve（stat 运行） | 43.0325s | 42.9994s | -0.08% |
  | task-clock | 1,324,202.64ms | 1,325,868.00ms | +0.1% |
  | cycles | 3.713e12 | 3.674e12 | -1.1% |
  | instructions | 9.324e12 | 9.259e12 | -0.69% |
  | LLC miss rate | 58.13% | 58.50% | +0.37pp |

- 采样中 malloc 1.93%→1.78%、cfree 1.73%→1.63%，`Parallel::transfer` total 27.98%→17.48%，说明缓存路径确实命中并减少了部分管理开销；但 `MPI_Waitall`、progress 和总墙钟没有稳定改善。
- `perf record` 第二次运行 Total Evolve 42.9962s → 43.1940s，方向相反；按 paired stat 和 1% 接受门槛，**拒绝 #17**。#14–#17 全部保留为负结果，当前正式基线仍为 #13。

## 优化 #18：kodis/lopsided 批处理（2026-08-18）

- 记录状态：`正式计时`。短测 perf gate 通过后提交 40 步正式任务，FINAL: PASS，已接受为新的当前有效最优。
- 优化内容：compute_rhs 末尾的 24 次 lopsided + 24 次 kodis 调用每次都做整块扫描、各自建 halo 副本、每点单独判模板阶数。批处理把 B 个场打包进一次扫描：一次 halo setup，每点一次模板阶数判断，lopsided 每点只读一次 shift。
  - **问题根因**: 24 次 lopsided 调用共享同一组位移数组 betax/betay/betaz；24 次调用的模板阶数只和 (i,j,k) 有关，与场无关。每次扫描重复建立 halo、重复判断、重复读 shift。
  - **改动**: 新建 `src/batch_stencils.f90`，提供 `kodis_batch` / `lopsided_batch` 两个入口。编译开关 `AMSS_BATCH_STENCIL=ON`，运行时 env `AMSS_BATCH_STENCIL_B=2/4/8` 控制 batch 宽度。门控 `Symmetry<=1`（赤道对称），八分限走原 per-variable 路径。每点算术与 `kodis_loop.fh` / `lopsided_loop.fh` 活跃 `#else` 分支逐位一致；赤道 z 反射由 `zref()` 内联，匹配 `symmetry_bd`。
- 正式运行参数：MPI=30，OMP=1，Final_Evolution_Time=40.0，无 TwoPuncture cache（完整链路），`AMSS_ENABLE_OPENMP=ON`，`AMSS_RHS_LOCALITY=POINTWISE`，`AMSS_BATCH_STENCIL=ON`，`AMSS_BATCH_STENCIL_B=2`，`-O3`。
- 短测运行参数：MPI=30，OMP=1，Final_Evolution_Time=4.0，`--twop-cache`，`AMSS_RHS_LOCALITY=POINTWISE`，`AMSS_BATCH_STENCIL=ON`，B=2/4/8，`-O3`。
- 作业与产物：正式 job 115371（B=2，40 步，FINAL: PASS）；脚本 `run_batch_stencil_full_cluster_20260818_a.sh`；完整日志 `optimization_logs/opt18_batch_b2_full_20260818_a.run.log`。短测 job 115236（首次提交，B=2，trajectory RMS=35% 失败）；job 115311（修复后，B=2/4/8 全部通过）。
- 正式运行结果（job 115371，40 步）：**This Program Cost=688.18s**（精确 688.179s），Total Evolve=397.92s，Total Running Time=401.03s，40 步全部完成。
- 正式正确性（job 115371）：**FINAL: PASS**；Trajectory RMS=0.000000%，matched 40/40 golden timesteps（完整，非 prefix）；Constraint maxima level 0：Ham=0.27739667，Px=0.028132512，Py=0.031488238，Pz=0.026503396，全部 PASS（≤2.0）。
- 对比基线（#13，698.59s）：端到端 **-10.4s（-1.5%）**，Total Evolve **-13.8s（-3.4%）**。短测 4 步 B=2 的 44.13→41.15s（-6.8%）是 perf-gate 证据，正式 40 步为权威结果。
- 短测正确性：B=2/4/8 均 `FINAL: PASS`，Trajectory RMS=0，Ham=0.22588831，constraint.dat 与 bssn_BH.dat 与基线逐位一致（snapshot 仅行首时间戳注释差异，剥除后 md5 全匹配）。作业 STATE=Failed 是 `set -e`+`pipefail` 下 `diff | head` 管道的脚本假象。
- `perf stat` 配对结果（baseline #13 vs B，4 步短测）：

  | 指标 | #13 baseline | B=2 | B=4 | B=8 |
  |---|---:|---:|---:|---:|
  | Total Evolve (ABE) | 44.13s | 41.15s | 41.54s | 43.41s |
  | perf wall | 59.57s | 54.12s | 54.70s | 56.11s |
  | task-clock (msec) | 1,462,259 | 1,345,671 | 1,364,117 | 1,405,194 |
  | cycles | 4.078e12 | 3.779e12 | 3.826e12 | 3.943e12 |
  | instructions | 1.005e13 | 9.715e12 | 9.608e12 | 9.539e12 |
  | IPC | 2.46 | 2.57 | 2.51 | 2.42 |
  | LLC-load-misses | 1.100e10 (56.51%) | 1.070e10 (54.66%) | 1.080e10 (55.72%) | 1.115e10 (57.37%) |

- profiling/分析（self% 热点，baseline → B=2）：
  - `kodis_._omp_fn.{1,2,3}`：5.45% → 合并为 `__batch_stencils_MOD_kodis_batch._omp_fn.0` = 5.01%
  - `lopsided_._omp_fn.{1,2,3}`：5.93% → 合并为 `__batch_stencils_MOD_lopsided_batch._omp_fn.0` = 7.31%（self% 升因为总 cycles 分母降得更快）
  - `__rhs_halo_MOD_halo_get`：1.16% → 0.18%（24 次 halo rebind 折成约 12 次）
  - `compute_rhs_bssn_fused_`：16.60% → 18.33%（dominator 效应，绝对 cycles 下降）
- 结论：**接受 #18**。短测 perf gate 通过（B=2 Total Evolve -6.8%，cycles -7.3%，instructions -3.3%，LLC miss rate -1.85pp）；正式 40 步 FINAL: PASS，This Program Cost 688.18s（vs #13 698.59s，-1.5%），Total Evolve 397.92s（-3.4%），40/40 golden timesteps 完整匹配。B=4/8 短测收益回吐（B=8 近基线），batch 过宽把工作集挤出 L1/TLB，与 #10 POINTWISE 同一陷阱，故 B=2 为正式采纳值。#18 batch stencil B=2 成为新的当前有效最优。
- 备注：
  - 首次提交 job 115236 B=2 trajectory RMS=35% 失败，根因是 24 个场中 10 个（gxy/gyz/Axy/Ayz/Gamx/Gamz/betax/betaz/dtSfx/dtSfz）的 z-parity 符号取反。逐条解码原代码 SSS/AAS/ASA/SAA/ASS/SAS/SSA 三元组取 z 分量，0/24 mismatch 后修复。
  - 作业 STATE=Failed 为脚本假象：`set -e`+`pipefail` 下 `diff | head` 在 snapshot 行首时间戳不同时返回非零。剥除注释头后 md5sum 全匹配。
- 下一步：#18 已接受为新的当前有效最优。RHS 操作序重排候选（#19）是下一项，已提交集群正式 4 步 paired perf/正确性任务（job 115428，脚本 `profile_rhs_reorder_cluster_20260818_a.sh`）。
- DevPod 调试记录：见 `DEVPOD_LOG.md` 对应条目。

---

## 优化 #19：RHS 运算顺序重排（perf 中性，拒绝，2026-08-18）

> #18 之后我盯上 compute_rhs_bssn_fused_ 的栈压力：循环体里 gup*、18 个 Christoffel、18 个 first-kind connection、6 个 Ricci、6 个 fxx..fzz 等约 70 个标量同时 live，超出 AV 寄存器后栈溢出。假设随意运算顺序导致中间量长 live，重排让“一个 RHS 完成后尽早 store，专属临时量尽快死亡”能压低同时 live 的向量数。

- 记录状态：`4 步 perf 探索`；正确性通过（bit-for-bit IDENTICAL），性能持平（<1%，噪声），代码回退到 #18 基线。
- 优化内容：在 `compute_rhs_bssn_fused_` 的 fused do k,j,i 循环内做纯 lifetime-scheduling 重排，不改任何公式、不动模块内部 FMA 顺序。新增编译开关 `AMSS_RHS_REORDER`（默认 OFF）。
  - **Phase 1**：浅 gauge RHS（`Lap_rhs = -2*alpha*K`，`beta_rhs = FF*dtSf`）从循环末尾提到 chi_rhs/gij_rhs 之后立即 store。这些只依赖 loop-invariant 输入（trK、dtSfx/y/z），无下游消费，早存可让 alpn1_s/dtSf 读尽早死亡。
  - **Phase 4**：Gamma_rhs 组装完成（含 shift 项）后立即 store `Gamx_rhs/Gamy_rhs/Gamz_rhs` 并完成 `dtSf_rhs = Gamma_rhs - eta*dtSf`，让 Gamx_rhs_s/Gamy_rhs_s/Gamz_rhs_s 标量早死。注意 Gamxa_s/Gamya_s/Gamza_s 被 Ricci 的 Γ^i·Γ_{ijk} 项复用，不能在此处死，保留。
  - **Phase 5/6 顺序保持**：trK_rhs 必须在 Ricci 之后——Christoffel 在代码里是同一组标量被 conformal（Ricci 用）和 physical（lapse Hessian→trK_rhs 用）复用，physical 修正会覆盖 conformal，trK_rhs 前置会破坏 Ricci。BAM 段（trK_rhs + Aij_rhs）视为原子块整体移动，内部浮点顺序不动。
- 实现开关：`AMSS_RHS_REORDER=ON`。`CMakeLists.txt` 新增 `option AMSS_RHS_REORDER`；`bssn_rhs.f90` 用 `#ifdef AMSS_RHS_REORDER` 包裹 Phase 1/Phase 4 的 store 前置和循环末尾对应 store 的跳过。默认 OFF 时代码路径与 #18 完全一致。
- 测试脚本：`profile_rhs_reorder_cluster_20260818_a.sh`，baseline(#13 无重排) vs candidate(AMSS_RHS_REORDER=ON)，同节点 4 步 paired `perf stat`+`perf record`+`check.sh`+constraint/BH 逐位 diff。
- 作业与产物：job 115428；原始文件为 `perf_profiles/opt19_baseline_opt13.{stat.txt,data,run.log}`、`perf_profiles/opt19_candidate_rhs_reorder.{stat.txt,data,run.log}`；报告为 `perf_profiles/canonical_20260818_a/opt19_{baseline_opt13,candidate_rhs_reorder}_{self,total,dso}.txt`。
- 正确性：两个版本均 `FINAL: PASS`，Trajectory RMS=0。**constraint.dat 与 bssn_BH.dat 相对 baseline 逐位 IDENTICAL**（剥除行首时间戳注释后 md5sum 全匹配）。重排数值安全。
- `perf stat` 配对结果：

  | 指标 | #13 baseline | #19 candidate | 变化 |
  |---|---:|---:|---:|
  | Total Evolve（stat 运行） | 42.563s | 42.724s | +0.38% |
  | task-clock | 1,313,451.60ms | 1,318,825.91ms | +0.41% |
  | cycles | 3.670e12 | 3.680e12 | +0.26% |
  | instructions | 9.282e12 | 9.275e12 | -0.07% |
  | LLC miss rate | 56.75% | 56.56% | -0.19pp |

- `perf record` 的 Total Evolve baseline 42.69s → candidate 43.15s（+1.1%），与 stat 方向一致（均略升），属运行噪声。
- `self%`：`compute_rhs_bssn_fused_` 18.58% → 18.03%（-0.55pp，分母效应，绝对 cycles 略升），kodis/lopsided 几乎不动。
- 分析：重排数值安全（bit-for-bit IDENTICAL 证明 lifetime scheduling 未破坏依赖），但性能持平——wall-clock、cycles、task-clock 全在 ±0.4% 内。根因：gfortran -O3 自身在指令调度阶段已做良好寄存器分配，源码层语句重排未提供额外信息；且能早死的只有少数浅 RHS 标量，gup*(7)、Christoffel(18)、g*_t(18) 是真正长 live 的核心节点（跨 Ricci 不可早死），Phase 重排压不动峰值。#18 批处理已把 kodis/lopsided 访存压力压掉，compute_rhs 内部边际收益进一步缩小。
- 结论：**拒绝 #19**，不纳入正式成绩。代码保留在 `AMSS_RHS_REORDER` 宏后（默认 OFF）作为实验档案。当前正式基线仍为 #18。
- 备注：重排的 bit-for-bit 正确性验证了依赖分析（Gamxa_s 跨 Ricci、Christoffel conformal/physical 双态、BAM 原子块）的正确性，为后续若需在 #18 基线上做更深 cache-blocking 时提供安全的语句移动边界。
- 下一步：#18（688.18s）距 340s 目标仍差 ~348s，其中 TwoPuncture 286s 固定块占主导。演化 397.92s 若要继续压缩，需考虑 cache blocking（lev8 分块）或 TwoPuncture 并行化（单进程串行，LineRelax/Thomas 未并行）。
- DevPod 调试记录：无（devpod 共享机过慢，4 步需 ~13min/步，正确性验证直接上集群 job 115428）。

## 优化 #22/#23：关闭 RHS NaN 扫描并让细层 predictor 走演化路径（2026-08-19）
- 记录状态：`4 步 perf 探索`；正确性通过，作为阶段一临时累计基线。
- 改动：新增 `AMSS_RHS_NAN_CHECK` 和 `AMSS_FINE_PREDICTOR_EVOLVE`；level 0 predictor 保留 constraint 模式，level>0 predictor 改走 POINTWISE 演化模式。
- 作业：115792；baseline `opt22_baseline_opt18`，candidate `opt23_candidate_nan_off_finepred`。
- `perf stat`：Total Evolve `43.6576s → 40.5572s`（-7.1%）；cycles `3.985e12 → 3.718e12`；LLC miss rate `57.11% → 55.62%`。
- `perf record`：Total Evolve `43.2405s → 40.5427s`，方向一致。
- 正确性：两版本均 `FINAL: PASS`，Trajectory RMS=0，constraint maxima 逐位一致。
- 结论：接受为阶段一临时累计基线。

## 优化 #24：legacy fused RHS workspace pool（2026-08-19）
- 记录状态：`4 步 perf 探索`；首次实现因 workspace shape 复用错误拒绝，修复后接受。
- 首次作业：115801；candidate constraint Ham=802.23342，立即拒绝。
- 修复作业：115807；按当前 Block shape 精确重建 pooled workspace，避免旧大 shape 的 pointer extent 泄漏。
- 改动：`AMSS_RHS_WORKSPACE_POOL=ON` 时，legacy/constraint 分支的111个导数数组进入 module workspace；POINTWISE 演化调用不再逐次创建这些数组。
- `perf stat`：Total Evolve `40.5436s → 40.1077s`（-1.08%）；dTLB miss `174.99B → 162.84B`（-6.9%）。
- 正确性：修复后 `FINAL: PASS`，Trajectory RMS=0，constraint maxima 通过。
- 结论：接受修复后的 workspace pool。

## 优化 #25：vacuum BSSN 专用化（2026-08-19）
- 记录状态：`4 步 perf 探索`；正确性通过，阶段一检查点组成项。
- 作业：115817；baseline `opt25_baseline_opt24`，candidate `opt25_candidate_vacuum`。
- 改动：`AMSS_VACUUM_BSSN=ON` 编译掉恒为零的 `rho/S_i/S_ij` 元素加载，保留原始 ABI 和字段布局。
- `perf stat`：Total Evolve `40.8602s → 39.5823s`（-3.13%）；cycles `3.474e12 → 3.432e12`；LLC miss rate `55.46% → 54.81%`。
- 正确性：两版本均 `FINAL: PASS`，Trajectory RMS=0，constraint maxima 通过。
- 结论：接受真空专用化。

## 阶段一正式检查点（2026-08-19）
- 作业：115838；配置包含 #22/#23/#24/#25，B=2，FUSED tail 关闭。
- 结果：`This Program Cost=681.446s`，`Total Evolve=381.533s`，`Total Running=385.272s`。
- 正确性：40/40 time groups，Trajectory RMS=0，`FINAL: PASS`；level 0 Ham=0.27739667，Px=0.028132512，Py=0.031488238，Pz=0.026503396。
- 对比 #18：端到端 `688.18s → 681.446s`（-0.98%），演化 `397.92s → 381.533s`（-4.12%）。

## 优化 #26：lopsided+KO B=2 尾部融合（2026-08-19）
- 记录状态：`4 步 perf 探索`；正确性通过，接受为阶段二临时累计基线。
- 作业：115828；baseline `opt26_baseline_opt24`，candidate `opt26_candidate_fused_tail`。
- 改动：B=2 batch 路径中把 lopsided 和 KO 合并为一次 RHS 扫描；保留原 lopsided 累加顺序，KO 在同一 RHS 值上追加。
- `perf stat`：Total Evolve `41.125s → 39.485s`（-3.99%）；L1D miss `0.73% → 0.73%`，dTLB miss `3.67% → 3.66%`。
- 正确性：两版本均 `FINAL: PASS`，Trajectory RMS=0，constraint maxima 通过。
- 结论：接受 B=2 fused tail。

## 优化 #27：fused RHS i-loop SIMD directive（2026-08-19）
- 记录状态：`4 步 perf 探索`；正确性通过，性能拒绝。
- 作业：115855；baseline `opt27_baseline_opt26`，candidate `opt27_candidate_bulk_simd`。
- 编译器报告仍显示主循环因 unsupported control flow 未向量化。
- `perf stat`：Total Evolve `40.7869s → 41.5044s`（+1.76%）；task-clock +1.46%，cycles +1.24%，dTLB miss +3.2%。
- 结论：拒绝；简单 SIMD directive 不改变主循环向量化结果，且引入额外开销。

## 优化 #28：fgfs cache coloring（2026-08-19）
- 记录状态：`4 步 perf 探索`；首次 slab 实现 SIGSEGV，改为独立 allocation+cache-line offset 后完成隔离测试，性能拒绝。
- 作业：115862 首次 slab；115883 隔离；115912 修复后的独立 allocation。
- 修复后 baseline `opt28j_baseline_opt25`，colored B=2 `opt28j_candidate_colored_b2`，colored B=4 `opt28j_candidate_colored_b4` 均 `FINAL: PASS`。
- `perf stat` baseline B=2 Total Evolve `38.3081s`；colored B=2 `41.7931s`；colored B=4 `41.4852s`。
- colored B=2 L1D miss `0.76% → 0.94%`，dTLB `3.69% → 3.91%`；B=4 L1D `1.11%`。
- 结论：拒绝 coloring；未纳入最终配置，开关默认关闭。

## 优化 #29：j×k tile sweep（2026-08-19）
- 记录状态：`4 步 perf 探索`；四组正确性通过，性能均未达到1%门槛。
- 作业：115936；baseline `opt29_baseline_best`，测试 `(TJ,TK)=(2,2),(4,2),(4,4),(8,2)`。
- baseline Total Evolve `38.7964s`；四组分别为 `39.4128s`（+1.59%）、`38.8913s`（+0.24%）、`38.8336s`（+0.10%）、`38.9240s`（+0.33%）。
- 结论：拒绝所有 tile 配置；普通 j×k 重排未形成稳定收益，最终保持未分块遍历。

## 阶段二最终正式检查点（#22–#26，2026-08-19）
- 作业：115969；MPI=30、OMP=1、Final_Evolution_Time=40.0、无 TwoPuncture cache、POINTWISE、batch B=2、workspace pool、vacuum specialization、fused tail。
- 结果：`This Program Cost=668.825s`，`Total Evolve=377.374s`，`Total Running=381.304s`，40步全部完成。
- 正确性：matched 40/100 golden timesteps（本次正式目标为40步），Trajectory RMS=0，40 time groups，`FINAL: PASS`；level 0 Ham=0.27739667，Px=0.028132512，Py=0.031488238，Pz=0.026503396。
- 对比 #18：端到端 `688.18s → 668.825s`（-2.81%），Total Evolve `397.92s → 377.374s`（-5.16%）。
- 结论：#22/#23/#24/#25/#26 累计路径成为新的正式 CPU 基线；#27/#28/#29 不纳入最终配置。

## 优化 #30：TwoPuncture 连续 JFD、workspace 复用与原地 Thomas（2026-08-19）
- 记录状态：`solver-only paired perf`；通用输入，未使用 cache，NRELAX=200、OMP=1。
- 作业：116059；baseline `build_cache/TwoPunctureABE`，candidate 为连续矩阵和复用 line workspace 版本。
- 改动：JFD/列索引改为连续行存储；5个 line scratch 在一次 BiCG 求解中复用；Thomas 算法原地消元；JFD 点级 `derivs` 改为栈上 nvar=1 临时量。
- `perf stat`：`286.792s → 277.464s`（-3.25%）；instructions -7.1%，L1D miss -26.4%，dTLB miss -25.0%。
- 正确性：ADM/Newton 日志一致，去除时间头后的 `Ansorg.psid` 数值逐行一致。
- 结论：接受，作为后续 TwoPuncture 基线。

## 优化 #31：TwoPuncture 谱表/几何缓存与冗余 final Newton 短路（2026-08-19）
- 记录状态：`solver-only paired perf`；NRELAX=200、OMP=1。
- 作业：116083；phase1 版本 `275.7046s`，candidate `160.9721s`。
- 改动：缓存 Chebyshev/Fourier sin/cos 表、A/B/phi 几何表；删除恒零 source 数组；当 ADM 误差和 Newton 残差均达标时跳过冗余 final Newton。
- `perf stat`：`275.70s → 160.97s`（-41.6%）；instructions -55.7%。
- 正确性：ADM 误差通过，Newton 残差 `7.53e-13`；`Ansorg.psid` 差异约 `2e-11` 绝对量级，保留端到端数值门槛。
- 结论：接受。

## 优化 #32：TwoPuncture NRELAX sweep（2026-08-19）
- 记录状态：`solver-only sweep + serial paired perf`；job 116097 sweep，job 116121 paired。
- sweep 结果（秒）：`NRELAX=1:96.47`，`2:92.04`，`5:78.28`，`10:74.60`，`25:85.87`，`50:103.18`，`100:124.11`，`200:167.02`。
- paired：`160.020s → 69.844s`（NRELAX=10，-56.4%）。
- 改动：新增 `AMSS_TWOP_NRELAX`；低 sweep 不收敛或达到迭代上限时回退到历史 `NRELAX=200`。
- 正确性：最终 ADM/Newton 通过，`Ansorg.psid` 数值差异仍为浮点级。
- 结论：正式候选选择 NRELAX=10。

## 优化 #33：TwoPuncture 独立 OpenMP（2026-08-19）
- 记录状态：`solver-only paired perf`；30个物理核，BSSN 仍保持 MPI=30/OMP=1。
- 作业：116184；baseline NRELAX=10、OMP=1，candidate NRELAX=10、TwoPuncture OMP=30。
- 改动：谱导数独立线、`J_times_dv` 的 i 线和红黑 line relaxation 并行；每线程使用独立 scratch；新增 `AMSS_TWOP_OPENMP` 和 `AMSS_TWOP_THREADS`。
- 初版曾因共享 `J_times_dv` 临时量 SIGSEGV，修复 private 作用域后重新作业通过。
- `perf stat`：`69.880s → 8.472s`（-87.9%，8.25×）；candidate task-clock 约25核。
- 正确性：ADM/Newton 通过，solver 输出有限，端到端 gate 后确认。
- 结论：接受30线程版本；旧 SIGSEGV 二进制拒绝。

## 优化 #34：TwoPuncture full-chain 4步 gate（2026-08-19）
- 作业：116194；无 cache，baseline 为 #30 版本，candidate 为 #30–#33 累计版本，TwoPuncture=30线程/NRELAX=10。
- 结果：`This Program Cost 318.001s → 51.786s`；baseline/candidate `Total Evolve=38.435s/38.845s`，演化差异约1%。
- 正确性：两侧 `FINAL: PASS`；candidate trajectory RMS `1.03e-5%`，constraints PASS。
- 结论：允许进入正式40步验证。

## 优化 #35：TwoPuncture 正式40步（2026-08-19）
- 作业：116215；MPI=30、BSSN OMP=1、TwoPuncture OMP=30、NRELAX=10、无 cache、Final_Evolution_Time=40.0。
- 结果：`This Program Cost=411.257s`，`Total Evolve=397.011s`，`Total Running=401.019s`；TwoPuncture 初值约8.5s。
- 正确性：matched 40/100 golden groups；Trajectory RMS `4.05e-6%`；level 0 Ham=0.27739668、Px=0.028132512、Py=0.031488237、Pz=0.026503396，`FINAL: PASS`。
- 对比 #26：端到端 `668.825s → 411.257s`（约-38.5%）。本次 ABE 演化为397.011s，较 #26 正式记录有约5%运行波动；full 4步 paired 的演化差异仅约1%，不把该波动归因于TwoPuncture。
- 结论：#35 成为当前正式 CPU 基线；下一步回到 BSSN/AMR/analysis 热点。

<!-- 后续集群正式运行记录追加在此线之下 -->

## 优化 #36：`co=1` 跳过约束专用中间数组写回（2026-08-19）
- 作业：116964；MPI=30、OMP=1、4步、TwoPuncture cache、B=2；baseline/candidate 使用完全相同编译参数。
- 改动：在 `compute_rhs_bssn_fused` 中，`gup*`、`R*`、`Gam***` 共30个只供 `co=0` 约束块读取的数组，仅在 `co==0` 时写回；`co=1` 仍执行标量计算，不改变演化公式和ABI。
- 正确性：baseline/candidate 均 `FINAL: PASS`，Trajectory RMS=0，匹配4/100组；level-0 constraint maxima 均为 Ham=0.22588831、Px=0.02448088、Py=0.0091883125、Pz=0.014353656。
- `perf stat`：Total elapsed `57.8399s → 57.4351s`（-0.70%）；cycles `3.9832e12 → 3.9754e12`（-0.20%）；L1D miss `29.45e9 → 25.94e9`、LLC miss `6.769e9 → 6.535e9`、dTLB miss `180.35e9 → 167.05e9`。
- perf report 中 `compute_rhs_bssn_fused_` self share `14.55% → 13.45%`，但总演化下降未达到1%门槛。
- 结论：访存计数改善，墙钟收益不足，**不默认开启**；保留 `AMSS_RHS_SKIP_CONSTRAINT_STORES` 作为后续组合实验开关。

## 优化 #37：`chin1` 与对角度规辅助扫描消除（2026-08-19）
- 作业：120789；MPI=30、OMP=1、4步、TwoPuncture cache、B=2；baseline/candidate 使用完全相同编译参数。
- 改动：`co=1` 跳过 `chin1=chi+1` 整数组扫描；batch lopsided 对对角度规直接使用 `dxx/dyy/dzz`，利用 lopsided 系数和为0消除 `gxx/gyy/gzz` 三个辅助数组的生成与读取。#36 保持关闭。
- 正确性：两侧 `FINAL: PASS`，Trajectory RMS=0，匹配4/100组；level-0 constraint maxima 两侧均为 Ham=0.22588831、Px=0.02448088、Py=0.0091883125、Pz=0.014353656。
- `perf stat`：Total elapsed `54.1887s → 53.6783s`（-0.942%）；cycles `3.8242e12 → 3.7804e12`（-1.146%）；instructions -0.628%，L1D miss -2.629%，dTLB miss -1.881%，LLC miss绝对数 +0.891%。
- `Total Evolve` 两次 paired 运行分别为 `41.453s → 40.9057s`（-1.32%）和 `40.696s → 40.2313s`（-1.14%）。
- 结论：演化时间稳定达到1%门槛，**接受并默认开启**；下一步提交正式40步验证。

## 优化 #37 正式40步验证（2026-08-19）
- 作业：120900；默认开启 #37，MPI=30、OMP=1、TwoPuncture OMP=30、NRELAX=10、无 cache、Final_Evolution_Time=40.0。
- 结果：`This Program Cost=402.9969s`，`Total Evolve=390.744s`，`Total Running=393.718s`。
- 正确性：匹配40/100 golden groups，Trajectory RMS `4.05058975e-6%`，40 time groups、9 levels，`FINAL: PASS`；约束全部通过。
- 对比 #35：端到端 `411.257s → 402.997s`（-2.01%），演化 `397.011s → 390.744s`（-1.58%）。
- 结论：#37 进入当前正式 CPU 基线。

## 优化 #38：BSSN 演化 hybrid MPI+OpenMP（fused RHS 主循环，2026-08-19）

> 依据 [BSSN_HYBRID_OMP_PLAN.md](../BSSN_HYBRID_OMP_PLAN.md) Switch 1。历史 #2 给旧版 stencil 子程序加过 OpenMP（30×1 最快），但 #2 没覆盖 #8 才出现的 fused RHS 主循环；本轮专门并行该循环。

- 记录状态：`4 步 perf 探索`；正确性通过（private 列表完整、无竞争），性能负向，代码保留在开关后默认关闭，不进正式路径。
- 优化内容：在 `compute_rhs_bssn_fused_`（src/bssn_rhs.f90:1056–3325）的非 tiled k/j 循环（`:1492` `#else` 分支）外加 `!$omp parallel do collapse(2) schedule(static) if(omp_get_max_threads()>1 .and. pointwise_mode) private(<194 个点级标量>)`。`pointwise_mode=(co==1 .and. Symmetry<=1)`（`:1319`），故 co==0 约束路径和 octant 串行、无竞争。新增编译开关 `AMSS_RHS_OMP_ASSEMBLY`（CMakeLists.txt，默认 OFF）；`use omp_lib`（`if()` 子句在 `implicit none` 下需显式声明 `omp_get_max_threads`）。循环体逐字节未动。
- 实现细节：
  - private 列表由脚本从声明块 `bssn_rhs.f90:1209–1269` 生成（194 个），含所有 `*_s` 标量、`gxxx_t..gzzz_t`、`gxxx2_s` 组、`pw_*` 导数结果标量、`pw_order`、`SSS..SSA`、`i,j,k`。共享量（留 function scope）：`PI`（循环后复用）、`dX/dY/dZ`、`pw_d12x…pw_fyzc`、`pw_imin/jmin/kmin`、所有输入/输出/存储数组。
  - `if()` 子句双重门控：OMP=1 时为假→循环串行（无团队开销），co==0 时为假→约束路径串行。意图是开关 OFF 时 #37 逐位不变、开关 ON 且 OMP=1 时也不回退。
- 运行参数：阶段A paired，4 步，`--twop-cache`，`AMSS_RHS_LOCALITY=POINTWISE`、`AMSS_BATCH_STENCIL=ON` B=2、TwoPuncture OMP=30/NRELAX=10。
- 作业与产物：阶段A job 121399（baseline 30×1 OFF、cand 30×1 ON、cand 15×2 ON）；阶段A2 job 121433（baseline 15×2 OFF，用于同 rank 数 race-check）；脚本 `profilescripts/profile_hybrid_stageA_cluster.sh`、`profile_hybrid_stageA2_cluster.sh`、`hybrid_stageA2_diff.sh`；日志 `optimization_logs/hybrid_stageA_*.run.log`、`hybrid_stageA2_*.run.log`；运行输出 `perf_runs/{baseline_30x1_off,cand_30x1_on,cand_15x2_on,baseline_15x2_off}/`。
- 运行结果（4 步，Total Evolve）：

  | config | switch | MPI×OMP | bind | Total Evolve | This Program Cost |
  |---|---|---|---|---:|---:|
  | baseline 30×1 | OFF | 30×1 | — | 38.7196s | 39.888s |
  | cand 30×1 | ON | 30×1 | — | 40.6403s | 42.317s |
  | cand 15×2 | ON | 15×2 | `--bind-to core --map-by slot:pe=2` | 108.839s | 110.378s |
  | baseline 15×2 | OFF | 15×2 | `--bind-to core --map-by slot:pe=2` | 134.291s | 135.778s |

- 正确性（race-check，关键）：
  - 30×1 ON vs 30×1 OFF：`FINAL: PASS`，Trajectory RMS=0，level-0 constraint maxima `Ham=0.22588831, Px=0.02448088, Py=0.0091883125, Pz=0.014353656` 与 OFF 逐位一致。
  - 15×2 ON vs 15×2 OFF（同 rank 数、仅开关不同）：剥除 `bssn_BH.dat` 时间戳注释头后数据行 **bit-identical**（`hybrid_stageA2_diff.sh` 输出 `BIT-IDENTICAL`）。→ private 列表完整，OMP=2 下无数据竞争。
  - 15×2 ON 对 30-rank golden：`FINAL: FAIL`，Trajectory RMS=12.443238%，但 constraint maxima 仍逐位为 `Ham=0.22588831` 等。→ 这是拿 15-rank 轨迹比 30-rank golden 的测试设计问题（rank 数变了→collective 模式变→轨迹合理发散），不是竞争、不是数值发散。
- 性能分析：
  - 30×1 ON vs OFF：演化 +5.0%（38.72→40.64s），超计划 §6 的 ±0.5% 无回归门。`if()` 为假时循环虽串行，但 `private(194 标量)` + `collapse(2)` 改了 gfortran 代码生成，串行路径也付开销（每 (k,j) 迭代的 private 栈设置/代码布局变化）。其中含节点负载成分（本次 baseline 38.72s 偏快于历史 #37 4 步的 ~40–41.5s），但 directive 静态开销确实存在。
  - 15×2 ON vs 15×2 OFF：演化 -19%（134.29→108.84s），OpenMP 在每 rank 的多块上确有并行收益（打破了 #2「块太小、同步开销吃掉收益」的旧结论——因为 #2 测的是旧 stencil 子程序，不是 fused 主循环）。
  - 15×2 ON vs 30×1：端到端慢 2.8×（110s vs 40s）。减半 rank 让每 rank 计算量翻倍，OMP=2 只能部分补偿，与 #1 rank 探索（9→30 单调加速）一致。
- 结论：**hybrid 物理核方向判定失败，拒绝**。硬门未过：30×1 ON 回退 +5% > ±0.5%；15×2 端到端慢 2.8×。private 列表正确性已验证（OMP=2 bit-identical），代码保留在 `AMSS_RHS_OMP_ASSEMBLY` 宏后默认 OFF，#37 正式基线不变。
- 备注：SMT 模式（30×2/15×4/10×6/6×10）随后单独测试（阶段C）；若 SMT 亦负向，hybrid OpenMP 方向整体放弃，回到 THP slab / LLC-dTLB 备选。
- 下一步：阶段C SMT sweep；无论结果如何不修改正式默认路径。

## 优化 #38 阶段C：SMT hybrid OpenMP sweep（2026-08-19）

> 阶段A 物理核 hybrid 失败后，按 BSSN_HYBRID_OMP_PLAN.md §3.2/§5 阶段C 补完 SMT 矩阵。SMT 保持 30 物理核、用满 60 逻辑 CPU：每 rank 的 OMP 线程落到其物理核的 SMT sibling。理论动机是 SMT 隐藏内存等待。

- 记录状态：`4 步 perf 探索`；正确性通过（constraint 逐位一致、无竞争），性能负向，hybrid OpenMP 方向整体放弃。
- 运行参数：4 步，`--twop-cache`，switch ON（`AMSS_RHS_OMP_ASSEMBLY=ON`），POINTWISE/B=2/TwoPuncture OMP=30/NRELAX=10。SMT 绑核：`--bind-to core --map-by slot:pe=<phys/rank>`，`OMP_PLACES=threads`、`OMP_PROC_BIND=spread`。编译一次共用 `build_hybrid_smt`。
- 作业与产物：job 121465；脚本 `profilescripts/profile_hybrid_stageC_smt_cluster.sh`；日志 `optimization_logs/hybrid_stageC_cand_*_smt_on.run.log` + `hybrid_stageC_compile.log`；运行输出 `perf_runs/cand_{30x2,15x4,10x6,6x10}_smt_on/`。
- 运行结果（4 步，Total Evolve；参考 30×1 OFF #37 = 38.72s/4step，阶段A job 121399）：

  | config | MPI×OMP | phys/rank (pe) | Total Evolve | This Program Cost | check.sh | Ham (level 0) |
  |---|---|---:|---:|---:|---|---|
  | cand 30×2 SMT | 30×2 | 1 | **39.9813s** | 41.336s | FINAL: PASS, RMS=0 | 0.22588831 |
  | cand 15×4 SMT | 15×4 | 2 | 84.7801s | 86.262s | FAIL(RMS 12.44%, rank≠30) | 0.22588831 |
  | cand 10×6 SMT | 10×6 | 3 | 86.4326s | 87.876s | PASS(matched prefix 4/100, RMS=0) | 0.22588831 |
  | cand 6×10 SMT | 6×10 | 5 | 107.379s | 108.949s | PASS(matched prefix, RMS=0) | 0.22588831 |

- 正确性：四个配置 level-0 constraint maxima 全部逐位 `Ham=0.22588831, Px=0.02448088, Py=0.0091883125, Pz=0.014353656`（与 #37 一致）→ SMT 下无竞争、无数值发散。30×2（rank=30）对 golden `FINAL: PASS, RMS=0` 是直接判定通过。15×4 的 `FAIL RMS=12.44%` 与阶段A 15×2 物理核完全相同（同一 15-rank 轨迹比 30-rank golden 的测试设计问题，rank 数变→collective 模式变→轨迹合理发散）；10×6/6×10 的 `PASS RMS=0` 是 matched-prefix 4/100 的巧合（前 4 步 BH 位置恰在 `--time-tolerance 1e-8` 内匹配 golden），不表示 10/6 rank 等价于 30 rank。
- 性能分析：
  - **30×2 SMT（39.98s）比 30×1 ON（40.64s）略快 1.6%，但比 30×1 OFF（38.72s）仍慢 3.3%**。SMT 第二线程在 compute-bound 双精度 stencil 上几乎无贡献（与 #5 SVE lane 不增、#1 `--use-hwthread-cpus` 争抢更慢一致），仅摊掉了 OMP=1 ON 那部分 directive 静态开销，未带来净正收益。
  - 减 rank 明显变慢：30×2(39.98) → 15×4(84.78) → 10×6(86.43) → 6×10(107.38)。与阶段A 物理核趋势一致（#1 rank 探索：9→30 单调加速，减 rank 必慢），OMP 线程无法补偿每 rank 计算量翻倍。
  - 15×4 SMT（84.78s）比 15×2 物理核 ON（108.84s）快 22%——SMT 在减 rank 场景下能部分补偿（每 rank 块多，OMP 多线程有更多并行面），但绝对值仍是 #37 的 2.2×。
- 结论：**SMT 方向亦失败，hybrid OpenMP 整体放弃**。所有配置相对 #37（38.72s/4step）无稳定收益：30×2 慢 3.3%、其余慢 2.2–2.8×。`AMSS_RHS_OMP_ASSEMBLY` 保持默认 OFF，#37 正式基线不变。SMT 未隐藏 compute-bound 双精度 stencil 的内存等待（该等待本就是 LLC-miss 39% 的访存瓶颈，SMT sibling 共享 L1/L2，无法并行掩盖同级 cache miss）。
- 下一步：hybrid MPI+OpenMP/SMT 方向全部否决，回到 BSSN_HYBRID_OMP_PLAN.md §6 末尾的备选（THP slab / LLC-dTLB miss 降低），或其它不改变 rank/线程模型的访存优化。


## 阶段一：LLC/dTLB miss 归因（2026-08-19）

> BSSN_LLC_DTLB_THP_PLAN.md §阶段一。先归因 miss 来源，再决定阶段二/三是否值得做。#37 基线，4 步，30×1，`--twop-cache`，`-O3 -g -fno-omit-frame-pointer -march=native`。

- 记录状态：`诊断`（miss 归因，不计时）。
- 作业与产物：job 121544（perf stat ×2 + cycle 采样 perf record）；job 121570（首次 miss 事件采样，`-e LLC-load-misses` 未加 `-m` 超 516KB mlock 预算失败）；job 121586（修 `-m 4` 后成功）。原始文件 `perf_profiles/llc_dtlb_opt37_{llc,dtlb}.stat.txt`、`perf_profiles/llc_dtlb_opt37.data`、`perf_profiles/missattr_opt37_{llcmiss,dtlbmiss}.data`；报告 `perf_profiles/canonical_20260818_a/{llc_dtlb_opt37,missattr_opt37_*}_{self,total,dso}.txt`。
- 聚合 miss（perf stat，证实计划前提）：

  | 指标 | 计数 | rate |
  |---|---:|---:|
  | LLC-load-misses | 6.636e9 | 50.37% |
  | dTLB-load-misses | 1.524e11 | 3.60% |
  | cache-misses | 2.849e10 | 0.73% |
  | task-clock | 1252s | 24.78 CPUs utilized |
  | IPC | 2.60 | — |

- 周期采样 self%（cycle 热点）：`libopen-pal 0xf13c4` 25.81%（opal 忙等，非访存）、`compute_rhs_bssn_fused_` 16.22%、`lopsided_batch` 15.06%、`polint_` 8.57%、`__memcpy_sve` 4.07%、`prolong3_` 2.84%、`point_d2` 2.54%、malloc/cfree 3.54%。total 显示 `GOMP_parallel` 15.07%（batch stencil 的 OMP 团队创建，OMP=1 时仍有 fork/join 开销，呼应 #38 阶段C）。
- **miss 事件采样 self%（直接归因，关键）**：

  | 函数 | LLC-load-misses self% | dTLB-load-misses self% |
  |---|---:|---:|
  | `rungekutta4_rout_` | **31.17%** | — |
  | `__memcpy_sve` | 23.09% | — |
  | `compute_rhs_bssn_fused_` | 12.64% | **28.37%** |
  | `lopsided_batch` | 10.86% | — |
  | `enforce_ga_` | 7.28% | — |
  | `average2_` | 4.80% | — |
  | libopen-pal（opal 忙等，非真实访存）| — | 51.50% |
  | libmpi 通信 | 0.76% | ~7% |

- 关键发现：
  1. **LLC miss 最大单一来源是 `rungekutta4_rout_`（31.2%），不是 RHS/lopsided**。RK4 是 whole-array triad（`f1=f0+a·dT·f_rhs`），40 个演化变量各一独立 malloc 数组，每子步流式扫 3 数组，lev8 工作集 ~20MiB 远超 per-core L3 → 首扫必 miss。计划 §阶段三前提「LLC miss 主要来自 RHS/lopsided」**被推翻**。
  2. **dTLB miss 的 51% 是 opal 忙等假象**（忙等轮询访问共享内存），真实访存里 compute_rhs 占 28.37% 是最大来源——与计划 §当前基线写的「167 个独立 malloc 不利于 THP 覆盖」一致，支持阶段二 THP slab 方向。
  3. memcpy（transfer pack/unpack）占 LLC miss 23%，是 AMR 通信的内存搬运。
- 结论：阶段二（THP slab）针对 dTLB 仍对路（dTLB 真实访存主源是 compute_rhs 扫 167 小数组）；阶段三（derivative tile）前提不成立（LLC 主源是 RK4 不是 RHS），暂缓。

## 优化 #39：THP slab（AMSS_FGFS_HUGEPAGE_SLAB，2026-08-19）

> BSSN_LLC_DTLB_THP_PLAN.md §阶段二。167 个独立 fgfs malloc 合并成一个 2MiB 对齐、MADV_HUGEPAGE 背书的 slab，降低 dTLB miss 和页表 walk。

- 记录状态：`4 步 perf 探索`；正确性通过（逐位一致），性能严重负向（+50% 回退），拒绝。开关默认 OFF，#37 不变。
- 优化内容：`src/Block.h` 新增 `double *fgfs_slab;`（switch 守护）；`src/Block.C` 构造时一次 `posix_memalign(2MiB)` 拿整块 slab，每字段 view 对齐到 2MiB 边界（`field_stride = round_up(nn*8, 2MiB)`），`madvise(MADV_HUGEPAGE)` 请求大页；析构 `free(fgfs_slab)` 一次；`swapList` 只换 view 指针，slab 所有权不动。`CMakeLists.txt` 新增 `option AMSS_FGFS_HUGEPAGE_SLAB OFF` + define 块 + 摘要打印。顺带修复 USE_GPU 分支悬空 `}` 的 `#endif` 错位。
- 运行参数：paired 4 步，30×1，`--twop-cache`，`-O3 -g -fno-omit-frame-pointer`。baseline OFF（#37）/ candidate ON。
- 作业与产物：job 121636；脚本 `profilescripts/profile_thp_slab_stage2_cluster.sh`；日志 `perf_profiles/thp_{baseline_off,candidate_on}.run.log` + stat/data；运行输出 `perf_runs/thp_{baseline_off,candidate_on}/`。
- 运行结果（4 步，Total Evolve 取自 `ABE_out.log`）：

  | 指标 | baseline OFF | candidate ON | 变化 |
  |---|---:|---:|---:|
  | Total Evolve | 38.5361s | 57.827s | **+50.0% 回退** |
  | task-clock (LLC run) | 1265s | 1936s | +53% |
  | cycles | 3.57e12 | 4.93e12 | +38% |
  | instructions | 9.62e12 | 1.06e13 | +10% |
  | IPC | 2.69 | 2.16 | -20% |
  | LLC-load-misses 绝对数 | 6.65e9 | 7.30e9 | +10% |
  | LLC-load-misses rate | 50.05% | 13.80% | -36pp（假象） |
  | dTLB-load-misses 绝对数 | 1.53e11 | 1.69e11 | +11% |
  | dTLB-load-misses rate | 3.54% | 3.50% | 基本不变 |

- 正确性：candidate `FINAL: PASS`，Trajectory RMS=0，level-0 constraint maxima `Ham=0.22588831, Px=0.02448088, Py=0.0091883125, Pz=0.014353656` 与 baseline **逐位一致**。slab 布局改对了，ABI 不变，无数据错误。
- 性能分析（拒绝依据）：
  1. **Total Evolve +50% 回退**，远超计划 §验收「≥1% 改善」门槛（方向反了 50 倍）。
  2. **LLC miss rate 从 50% 降到 13.8% 是假象**：candidate 跑得慢得多（57.8s vs 38.5s），分母（总访问）变小，rate 降低是数学假象；**绝对 miss 数反而 +10%**。task-clock/cycles 暴涨 +53%/+38%、IPC 从 2.69 跌到 2.16、instructions +10% 是真实信号：每条指令更慢、执行更多指令。
  3. **dTLB miss 完全没降**（3.54%→3.50%，绝对数 +11%）——THP 的核心目标未达成。
- 根因：`field_stride = round_up(nn*8, 2MiB)` 让每字段占满 2MiB（哪怕数据只有 0.24–0.84 MiB），167 字段 × 2MiB = 334MiB/Block，远超原来 40–141MiB。工作集暴涨，per-core L3（几 MiB）根本装不下；且字段间隔 2MiB 破坏了原有空间局部性（原来相邻字段地址连续，RK4/RHS 扫多字段时预取器能带动，现在每字段间隔 2MiB 空洞，预取失效）。THP 用大页减 TLB 项的收益被工作集膨胀 + cache 局部性破坏的代价完全压倒。dTLB miss 本就只有 3.5%（阶段一归因显示 51% 是 opal 忙等假象，真实 compute_rhs 只占 28%），优化空间有限。
- 结论：**THP slab 方向拒绝**。`AMSS_FGFS_HUGEPAGE_SLAB` 保持默认 OFF，#37 正式基线不变。代码保留在开关后作为实验档案。
- 备注：THP/大页对「工作集远大于 L3 + 字段间有空间局部性」的场景是负优化。若要降低 dTLB miss，应在不膨胀工作集的前提下做（如紧排 slab 而非 2MiB 对齐每字段，但那样 fgfs[i] 可能跨 2MiB 边界、大页背书失效——THP 与紧凑布局在 167 个小字段上本质冲突）。阶段三（derivative tile）前提在阶段一已被推翻（LLC 主源是 RK4 不是 RHS），亦暂缓。
- 下一步：LLC/dTLB 优化方向（THP slab、derivative tile）均不成立或负向，#37 维持。可考虑的剩余方向：RK4 整数组流式的 miss（31% LLC 来源，但属算法固有 whole-array triad，改写需动 RK4 数据流）、或不再追求访存优化转向其它。

## 优化 #40：BSSN 数值 kernel 候选（2026-08-20）

- 目的：验证 `batch stencil` 并行区、`enforce_ga` 点式化、RK4 下界融合、RK4 空间块化和 Sommerfeld 批处理；所有候选均由独立 CMake 开关控制，默认关闭，不改变 #37 正式路径。
- 开关：`AMSS_BATCH_STENCIL_SINGLE_REGION`、`AMSS_ENFORCE_GA_POINTWISE`、`AMSS_RK4_FUSE_LOWERBOUND`、`AMSS_RK4_BATCHED`、`AMSS_SOMMERFELD_BATCHED`。
- 正确性：五个单项候选均通过 tiny 单层/1-rank smoke，`bssn_BH.dat` 和 `bssn_constraint.dat` 去除注释时间头后与 OFF 逐项一致；batch/enforce 两个候选又完成 30-rank/1-step smoke，均 `FINAL: PASS`、Trajectory RMS=0、constraint maxima 一致。
- 30-rank/1-step 同节点配对（`TotalTime=1.0`，仅用于候选筛选）：

  | 配置 | Total Evolve | Total Running | 相对 OFF | 状态 |
  |---|---:|---:|---:|---|
  | OFF baseline | 558.307s | 575.405s | — | PASS |
  | batch single-region | 561.611s | 585.807s | Evolve +0.59% | 拒绝 |
  | enforce_ga pointwise | 522.406s | 542.587s | Evolve -6.43% | 保留候选 |

- tiny kernel 方向筛选：lowerbound 融合与 RK4 tiled 在 16×16×8/1-rank 上明显慢于 OFF；Sommerfeld batch 与 OFF 基本持平，因此这三项暂不安排长时间 30-rank 配对，开关保持关闭。
- 集群 4-step 配对：job `123751` 在 30 MPI ranks、1 OMP thread、TwoPuncture cache 开启的条件下，以 `AMSS_ENFORCE_GA_POINTWISE=ON` 完成运行；CMake 日志确认开关确实生效。候选结果为 `Total Evolve=38.6471s`、`Total Running=41.9494s`、`Program Cost=42.4853s`。作业内置的 `check.sh` 因相对路径重复拼接而误报失败，改用绝对路径复核后为 `FINAL: PASS`，Trajectory RMS=0、constraints PASS。
- 与 #38 的 30-rank/1-thread OFF 4-step 参考（`Total Evolve=38.7196s`）相比，Evolve 仅改善 **0.187%**，低于「≥1%」验收门槛，属于测量噪声范围；此前本地共享 DevPod 上看到的约 6.4% 不具备代表性。因此拒绝该优化，`AMSS_ENFORCE_GA_POINTWISE` 保持默认 OFF，#37 正式基线不变。

## 优化 #41：CPU fast-math（2026-08-20）

- 目的：验证在 CPU `ABE` 的 C++/Fortran 演化 target 上启用 `-ffast-math` 是否能减少浮点运算指令；不把该选项施加到 TwoPuncture 初值求解器或 GPU target。
- 实现：`CMakeLists.txt` 新增 `AMSS_FAST_MATH`，并由 `compile.sh` 转发；正式验证期间显式 `AMSS_ENABLE_GPU=OFF`、`AMSS_FAST_MATH=ON`。通过后将两个默认值设为 ON，仍可用 `-DAMSS_FAST_MATH=OFF` 回退。
- 短测 paired：job `123890`，MPI=30、OMP=1、4 step、TwoPuncture cache=ON；OFF/ON 均 `FINAL: PASS`、Trajectory RMS=0、constraint 逐位一致。

  | 配置 | Total Evolve | Total Running | Program Cost |
  |---|---:|---:|---:|
  | fast-math OFF | 40.1683s | 40.8018s | 41.5936s |
  | fast-math ON | 38.5437s | 39.1501s | 39.9417s |

  `Total Evolve` 下降约4.05%，`Program Cost` 下降约3.97%。
- 硬件计数器 paired：job `123914`，OFF/ON 的 Evolve 为 38.5693/37.9333s；cycles 下降3.33%，instructions下降2.31%，task-clock下降1.51%。LLC-load-misses 绝对数为 6.696e9→6.737e9、miss rate 50.35%→50.60%，没有访存收益；dTLB-load-misses 下降约2.97%。收益主要来自指令数和浮点调度，不应归入 LLC 优化。
- 正式验证：job `123942`，脚本 `profilescripts/run_opt41_fastmath_formal40_cluster.sh`，日志 `optimization_logs/opt41_fastmath_formal40`。CPU-only、MPI=30、OMP=1、B=2、TwoPuncture cache=OFF、40 step；`Total Evolve=379.041s`、`Total Running=381.984s`、`This Program Cost=390.893s`，Trajectory RMS=`4.05058975e-06%`，constraints PASS，`FINAL: PASS`。
- 结论：**接受 #41，`AMSS_FAST_MATH` 默认开启**。它不会改变默认的 GPU 编译路径（GPU 仍 OFF），TwoPuncture 仍使用不带该选项的 target。

## 优化 #42：global_interp 固定 uniform6 polint 核（2026-08-20）

- 目的：`perf record` 显示 `polint_` 是 analysis/interpolation 路径的高热点；检查是否能利用当前 `global_interp` 恒定 `ordn=2*ghost_width=6` 且每轴采样点为 `(0,1,...,5)` 的事实，去掉通用 Neville 插值的数组复制、循环边界和逐次调用开销。
- 实现：新增 `AMSS_POLINT_UNIFORM6`。在 `global_interp` 中仅当 `ORDN==6` 时调用 `polin3_uniform6`：先计算三轴 6 点 Lagrange 权重，再做固定大小的 6×6×6 tensor-product contraction。`polin3/polint` 原实现保留为其他 order 和其他调用者的 fallback；`interp_2`、AMR Restrict/Prolong、Sommerfeld 不走专用核。`dy` 在该调用路径未被消费者使用，因此专用核将其置零。
- 局部逐点检查：专用核与通用 `polin3` 对固定 6×6×6 数据、非整数坐标的绝对差为 `1.78e-15`；集群 4-step 输出中的 trajectory、constraint 和 golden 比较均通过。
- paired perf/correctness：job `123998`，脚本 `profilescripts/profile_opt42_polint_uniform6_cluster.sh`，结果目录 `perf_profiles/opt42_polint_uniform6/`。两套均 CPU-only、fast-math=ON、MPI=30、OMP=1、4 step、B=2、TwoPuncture cache=ON。

  | 配置 | 第一次 Total Evolve | 第二次 Total Evolve | check |
  |---|---:|---:|---|
  | uniform6 OFF | 39.4079s | 39.4224s | `FINAL: PASS` |
  | uniform6 ON | 33.1405s | 33.5164s | `FINAL: PASS` |

  candidate 两次均与 baseline 的 4-step trajectory RMS=0、level-0 constraints 逐位一致。
- `perf stat` 原始文件：`opt42_baseline_fastmath.stat.txt`、`opt42_candidate_uniform6.stat.txt`。基线/候选的 cycles 为 3.416e12/2.869e12（-16.03%），instructions 为 8.881e12/7.562e12（-14.85%），perf elapsed 为 49.063/42.945s（-12.47%）。LLC-load-misses 为 6.977e9/6.966e9（绝对数基本不变，rate 52.85%→52.70%）；dTLB-load-misses 为 1.461e11/1.379e11（绝对数约降5.61%，rate 3.59%→3.87%）。收益来自消除通用插值算术/临时量，不是依靠降低 LLC miss。
- 本地 `perf report`：`opt42_baseline_fastmath_self.txt` 中 `polint_` self=9.05%、`polin3_`=0.63%；`opt42_candidate_uniform6_self.txt` 中 `polin3_uniform6_`=0.14%，原 `polint_` 不再进入热点。调用树中 `global_interp_` 仍保留，说明只替换了内部数值核。
- 正式验证：job `124023`，脚本 `profilescripts/run_opt42_polint_formal40_cluster.sh`，日志 `optimization_logs/opt42_polint_formal40`。CPU-only、fast-math=ON、uniform6=ON、MPI=30、OMP=1、B=2、TwoPuncture cache=OFF、40 step；`Total Evolve=329.018s`、`Total Running=332.715s`、`This Program Cost=342.54874753952026s`，Trajectory RMS=`4.05058975e-06%`，constraints PASS，`FINAL: PASS`。相对 #37 的 390.744/402.997s，演化约快15.80%，端到端约快15.00%。
- 结论：**接受 #42，`AMSS_POLINT_UNIFORM6` 默认开启**。该开关只影响 CPU `ABE` 的 `global_interp` ordn=6 路径；其他插值和 GPU 代码保持原实现。
