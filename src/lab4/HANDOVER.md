# Lab4 工作交接文档

> 最后更新:2026-08-16
> 项目:HPC101 实验四 AMSS-NCKU 数值相对论程序优化(任务一 CPU)

---

## 一、当前状态与关键发现

### 最新进展:成功采到 30 rank 的 perf record

之前一直拿不到多 rank 的函数级热点(计算节点 mlock 预算限制),现用 **`-m 4`**(30核×16KB=480KB < 516KB 预算)成功采到(job 107382, profile5.sh, 520013 样本)。这是迄今最关键的一张 profile。

### 🔴 决定性发现:30 rank 下 76% 时间在 MPI 忙等

30 rank perf record 显示:

```
57.71%  libopen-pal.so.80.0.5  (OpenMPI OPAL 运行时, 进度轮询)
 5.53%  libmpi.so.40.40.7
 5.23%  libmpi.so.40.40.7
 3.94%  libopen-pal.so.80.0.5
 3.57%  libopen-pal.so.80.0.5
 0.86%  opal_progress
─────────────────
~76%   MPI 库(主要是忙等轮询, 非真实通信)
<2%    ABE 自己的计算代码
~1%    malloc/cfree
```

**结论:30 rank 不是最优配置。** 程序大部分时间不在算,而是在 MPI 库里反复轮询"消息来了没"(opal_progress 是 OpenMPI 的 progress engine)。真正的计算(compute_rhs_bssn 等)在 30 rank 下被淹没到 <2%。

### 之前的 MPI=1 profile(perf4)已删除

perf4 是 MPI=1 跑的,**看不到 MPI 等待、且 4 步放大了 TwoPuncture 占比**,有误导性。已删除相关文件,以 30 rank 的 perf5 为准。

### 文档确认:评测不固定 30 rank

查实验文档([index.md](docs/lab/Lab4-AMSS-NCKU/index.md)):`AMSS_NCKU_Input.py` 正式测试**只允许改 MPI/OpenMP/GPU 参数**,其他禁改。文档明确把 rank 数决策交给学生(715 行:"实际评测平台应当如何配置 MPI rank 和 OpenMP thread,需要结合节点核心数、NUMA 拓扑、绑定策略和程序负载来决定")。**30 只是 baseline 默认值,不是评测要求**。GPU 才有"建议设为 1",CPU 无此限制。

---

## 二、优化方向(基于新发现,重新排序)

| 优先级 | 方向 | 理由 |
|--------|------|------|
| 1️⃣ | **测试并减少 rank 数**(9/15) | 30 rank 下 76% 在 MPI 忙等。减 rank 直接减等待。1 行改动(`MPI_processes`),零成本验证。文档允许改 rank 数 |
| 2️⃣ | **加 OpenMP 补偿并行度** | 减 rank 会损失并行度,用 OpenMP 在 rank 内补。文档鼓励 MPI/OMP 权衡(思考题#2)。需 CMake 启用 `AMSS_ENABLE_OPENMP` |
| 3️⃣ | **BSSN 向量化**(SVE) | 计算虽在 30 rank 下占比小,但减 rank 后计算占比会回升,向量化仍有意义。SoA 布局利向量化 |
| 4️⃣ | 减 malloc/free | 预分配缓冲区,中收益 |
| 5️⃣ | MPI 通信优化(重叠/合并) | 难,且需更深 profile |

**建议第一项实际优化:测试 rank 数。** 改 `MPI_processes` 跑 9/15 对比单步时间。

---

## 三、Baseline 数据(参考点)

### 集群正式参数
- MPI=30, OMP=1, Final_Evolution_Time=40.0, 编译=g++/gfortran 14.2 + OpenMPI 5, -O3, 无绑核

### 结果(外推,详见 OPTIMIZATION_LOG.md 优化#0)
- **原始值**(job 107261, Timeout):跑到 step 34/40,平均单步 43.19s(σ=0.30s,极稳定)
- **外推值**(40步):40×43.19 + 固定开销~287s(TwoPuncture ~176s + 绘图 ~110s + setup ~2s)≈ **2015s ≈ 33.6 分钟**
- 正确性:baseline20(job 102714)RMS=0%,constraint max≤0.27,**PASS**
- ⚠️ 外推 ~2015s > 30min walltime → 若评测 walltime 也是 30min,40 步会 Timeout,优化需解决

### perf 数据
- **perf5(30 rank, -m 4, job 107382)**:520013 样本。76% 在 MPI 库忙等,计算 <2%。报告在 `perf_profiles/perf5_symbols.txt`、`perf5_self.txt`。⚠️ **该次运行输出数据/绘图是坏的**(perf 写 4.2GB 致 /workspace IO 过载,plot 画空图仅 1s),不能用于正确性检查;演化采样本身可信
- **perf_stat(30 rank, 4步)**:IPC=1.93, LLC-miss=39.41%, dTLB-miss=6.87%, page-faults 4.4M, CPUs utilized=11.7 是**全程平均**(含 ~287s 单进程绘图;演化阶段实际满核,详见 DEVPOD_LOG 采样#2)

---

## 四、项目代码框架

### 目录结构
```
src/lab4/
├── AMSS_NCKU_Input.py        # 输入文件(只能改 MPI/OMP/GPU 参数)
├── AMSS_NCKU_Program.py      # 运行 driver
├── compile.sh / run.sh / check.sh
├── CMakeLists.txt            # 有 AMSS_ENABLE_OPENMP/ARCH_FLAGS 等开关
├── src/                      # C++/Fortran/CUDA 源码
├── scripts/                  # Python 辅助脚本
├── golden/                   # 正确性真值(只有 bssn_BH.dat)
├── OPTIMIZATION_LOG.md       # 集群正式计时主日志(含范式+外推准则)
├── DEVPOD_LOG.md             # DevPod调试+perf采样副日志
├── HANDOVER.md               # 本文件
└── perf_profiles/            # perf 报告
    ├── perf5_m4.data         # 30 rank record 原始数据(4.2G)
    ├── perf5_symbols.txt     # 函数热点(干净版)
    ├── perf5_self.txt        # self time
    ├── perf5_m4.log          # 采样日志
    └── perf_stat.txt         # 30 rank perf stat 计数器
```

### 三个可执行文件
| 文件 | 作用 | 主要源文件 |
|------|------|-----------|
| TwoPunctureABE | 生成初值 | TwoPunctureABE.C, TwoPunctures.C |
| ABE | CPU 演化(任务一) | ABE.C, bssn_class.C, bssn_rhs.f90 |
| ABEGPU | GPU 演化(任务二) | bssn_gpu_class.C, *_gpu.cu |

### 数据流
```
AMSS_NCKU_Input.py → driver 生成 input.par → TwoPunctureABE 解初值(Ansorg.psid+裸质量)
→ 回填 input.par → mpirun -n N ./ABE 演化 → binary_output/*.dat → 复制+绘图
```

### 网格存储四级
```
cgh(levels=9, PatL[lev]) → Patch(逻辑盒子+ghost) → Block(切块给rank, fgfs) → fgfs[i][idx]
```
- `fgfs` = double**, 167 个变量各一个独立数组(SoA),idx=i+j*nx+k*nx*ny
- `ghost_width=3`, `buffer_width=6`, `min_width=12`(域分解最小切块)

---

## 五、工作范式与准则

### 两阶段流程(详见 OPTIMIZATION_LOG.md)
- **阶段一 DevPod 调试**:演化4.0、MPI=30、`--twop-cache`,只验证正确性,记 DEVPOD_LOG.md
- **阶段二 集群计时**:正式参数(40步,rank 数随优化配置)、无 cache,记 OPTIMIZATION_LOG.md

### 外推准则(步数不够时)
- **同时记录原始值和外推值**,缺一不可
- 外推前提:同一运行内每步稳定(baseline 34步 σ=0.30s 已验证)
- 固定开销≈287s = **TwoPuncture ~176s**(perf4 样本 21.4% × 824.6s,单进程,与 rank 无关)+ **绘图 ~110s** + setup ~2s;三次实测 286.9/287.9/287.2s。⚠️ 曾误写"绘图277+TwoPuncture10",已按 perf4 修正

### 提交坑
`hpc submit` 用提交时 cwd,必须先 `cd src/lab4` 或用 `--chdir`。

### perf 采样要点
- 计算节点 mlock 预算 516KB(全局总额,RO)
- 30 rank record:**只有 `-m 4` 成功**(480KB<516);默认 -m 和大 -m 都超;`-m 1` 太小丢样本
- 计算节点生成的 report 不解析符号,要回 devpod `grep -vE "unwind|get_proc_name"` 重解析

---

## 六、⚠️ 当前临时状态

`src/lab4/AMSS_NCKU_Input.py` 现在是 profiling 临时值,正式跑前要改回:
```python
# 第47行(当前临时,要改回):
Final_Evolution_Time = 100.0 if GPU_Calculation == "yes" else 4.0  ## TEMP profile5
# 改回:
Final_Evolution_Time = 100.0 if GPU_Calculation == "yes" else 40.0  ## FORMAL: 40 steps
```

---

## 七、关键记忆文件
- 项目记忆:`~/.claude/projects/-home-h3250104945-HPC101/memory/`(lab4-overview/run-environment/data-flow/perf-hotspots/output-files)
- 仓库日志:OPTIMIZATION_LOG.md(主)、DEVPOD_LOG.md(副)、HANDOVER.md(本文件)

## 八、任务二(GPU)尚未开始
- 环境 x86-5418Y devpod(家目录不跨地域,用 git 同步)
- GPU 用 host staging 通信(MPI_CUDA_AWARE=0)
- baseline 未跑
