# AMSS-NCKU 优化日志(集群正式计时)

> 本文件只记录**集群计算节点**上的正式运行:MPI=30、Final_Evolution_Time=40.0、完整链路(无 twop cache)。
> 用于评分对比和最终报告。每次优化在此追加一条。
> DevPod 调试记录见 [DEVPOD_LOG.md](DEVPOD_LOG.md)。

## 工作范式(两阶段)

每次优化遵循以下两阶段流程:

### 阶段一:DevPod 调试(验证程序无误,不计时)→ 记入 DEVPOD_LOG.md
- 演化:`Final_Evolution_Time = 4.0`(临时,4 步快速验证)
- MPI:`MPI_processes = 30`(与生产一致)
- 用 `./run.sh --twop-cache`(初值不变时跳过 TwoPuncture)
- 验证:能编译 + 能跑完 + `./check.sh` 返回 **PASS**
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
3. 固定开销(TwoPuncture + 绘图 + setup)≈ **287s** = **TwoPuncture ~176s**(perf4 样本:LineRelax+Thomas 21.4% × 总 CPU 824.6s;单进程,与 rank 数无关)+ **绘图 ~110s** + setup ~2s。三次独立实测:j102714 286.9s / j103849 287.9s / j107433 287.2s。⚠️ 曾误拆为"绘图277s+TwoPuncture10s"(把 103849 误当缓存命中),已按 perf4 修正。若优化改了 TwoPuncture 才需重测。
4. **评分/对比以外推值为准**(因为正式评测就是 40 步),但原始值保留以证明外推可靠。

## 集群运行记录模板(每次优化复制使用)

```markdown
## 优化 #N: <简短标题>  (YYYY-MM-DD)
- 优化内容: <针对什么瓶颈, 改了哪些文件:行, 做了什么>
- 集群运行参数: MPI=30, OMP=?, Final_Evolution_Time=40.0, 编译器=?, 编译选项=?, 绑核/NUMA=?
- 集群运行结果: This Program Cost = Xs
- 正确性: RMS=?%, constraint max(Ham/Px/Py/Pz)=?  (PASS/FAIL)
- 对比 baseline: 时间变化 ?s (?%), 正确性是否保持
- profiling 复查: <优化前后热点变化>
- 备注: <失败/异常/思考>
- DevPod 调试记录: 见 DEVPOD_LOG.md 优化 #N
```

---

## Baseline (参考点)

## 优化 #0: Baseline 记录  (2026-08-16)
- 优化内容: 无,原始代码 `-O3 -fno-strict-aliasing`
- 集群运行参数: MPI=30, OMP=1, Final_Evolution_Time=40.0, 编译器=g++/gfortran 14.2 (OpenMPI 5), 编译选项=-O3, 绑核=默认(未设)
- 集群运行结果:
  - **原始值(job 107261, Timeout)**: 实际跑到 **step 34/40**,状态 Timeout(30min walltime 用尽)
  - 实测每步时间(step 1–34, 单位秒):
    43.16, 43.23, 43.42, 43.38, 43.58, 43.59, 43.78, 43.83, 43.66, 43.46,
    43.13, 42.96, 43.05, 42.87, 43.03, 43.33, 43.35, 43.24, 43.01, 43.00,
    43.36, 42.92, 43.15, 42.88, 42.90, 43.06, 43.06, 43.32, 43.60, 43.24,
    43.06, 42.65, 42.60, 42.91
  - 平均单步 = **43.19s**(σ ≈ 0.30s,极稳定 → 外推可靠)
  - **外推值(40 步)**: 40 × 43.19 = **1727.6s** 演化
    + 固定开销 ≈ 287s(实测:TwoPuncture ~176s + 绘图 ~110s + setup)
    = **外推总时间 ≈ 2015s ≈ 33.6 分钟**(与 Timeout 吻合)
- 正确性: 未跑完 40 步,无 check 结果;baseline20(job 102714, 20步)的 RMS=0.000000%, constraint max(Ham=0.269,Px=0.028,Py=0.021,Pz=0.022) **PASS**(证明代码正确,只是 walltime 不够)
- 对比 baseline: —(本条即基准)
- profiling 复查(perf_stat job 103849, 4步 MPI=30):
  - IPC=1.93, LLC-miss=39.41%, dTLB-miss=6.87%, page-faults 4.4M, CPUs utilized=11.7(⚠️ 全程平均,含 ~287s 单进程绘图;演化阶段实际满核,详见 DEVPOD_LOG 采样#2)
  - perf4 函数热点(job 104167, 4步 MPI=1, 81638样本):
    compute_rhs_bssn 16.7%, polint_ 9.9%, LineRelax_be 9.6%, __cos 8.4%, kodis_ 8.0%, LineRelax_al 6.7%, fdderivs_ 5.8%, Thomas 5.1%, lopsided_ 4.2%, malloc/cfree ~7%
  - 模块: BSSN演化~40%, TwoPuncture~23%, 内存管理~11%, 分析~10%, 数学库~10%
  - 负载不均衡: level0 只用 9/30 rank(40点切 min_width=12 只能切9块)
- 备注:
  - ⚠️ perf4 是 MPI=1,演化单步 ~131s(比 30 rank 的 43s 慢 3 倍,因为没并行);函数占比在"计算部分"可信,但"通信/等待占比"看不到。30 rank 的 MPI 等待占比已由 **perf5(job 107382, 30 rank)确认:~76% CPU 周期在 MPI 库忙等轮询**,即每步 ~33s 在 opal_progress(见 DEVPOD_LOG 采样#4)。
  - 固定开销 ~287s 里 **TwoPuncture ~176s**(40 步时 ~9%,单进程串行,文档明确优化对象:OpenMP 并行 LineRelax/Thomas)+ **绘图 ~110s**(~5%,文档范围外)+ setup ~2s。⚠️ 曾误写"绘图277+TwoPuncture10",按 perf4 样本占比(21.4% × 824.6s)修正。
  - 正式评测若也是 30min walltime,40 步外推 ~2015s 会 Timeout → **这是优化必须解决的核心矛盾**:要么提速到 30min 内,要么评测可能有不同的步数/walltime 配置(待评测说明确认)。
- DevPod 调试记录: 无(baseline 不需要调试)

---

<!-- 后续集群正式运行记录追加在此线之下 -->
