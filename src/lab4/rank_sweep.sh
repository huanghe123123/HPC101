#!/bin/bash
# rank_sweep.sh — 对比不同 MPI rank 数的单步演化时间(快速探索用)
#
# 用法: bash rank_sweep.sh 9 15:2 30   (rank:omp, 默认 omp=1)
# 行为: 对每个 rank/omp 组合, 备份 → 改 AMSS_NCKU_Input.py 的 MPI_processes 和
#        OMP_threads → 跑 ./run.sh --twop-cache (4 步, 初值缓存命中) → 恢复原 Input.py。
# 指标: 演化时间取 run.sh 日志的 "Total Evolve Time" (不受绘图异常影响)。
# 环境: 输出节点 lscpu / numactl, 供 NUMA/绑核方向分析。
# 安全: Input.py 用 cp 备份, trap EXIT 恢复; 不改其他任何文件。
set -uo pipefail

LAB=/home/h3250104945/HPC101/src/lab4
[[ -d "$LAB" ]] || LAB="$(pwd)"
cd "$LAB" || { echo "FATAL: cannot cd to $LAB"; exit 1; }
export AMSS_CACHE_DIR="${AMSS_CACHE_DIR:-$LAB/twopuncture_cache}"

INPUT=AMSS_NCKU_Input.py
[[ -f "$INPUT" ]] || { echo "FATAL: $INPUT not found"; exit 1; }
cp "$INPUT" "$INPUT.bak"
restore() { cp "$INPUT.bak" "$INPUT"; rm -f "$INPUT.bak"; }
trap restore EXIT

LOG=rank_sweep.log
echo "################ rank_sweep $(date -Is) ################" | tee -a "$LOG"

echo "==> 节点环境 (NUMA/绑核方向)" | tee -a "$LOG"
lscpu | grep -E "^CPU\(s\)|^Socket|^NUMA|^Core|^Thread|Model name" | tee -a "$LOG"
numactl -H 2>/dev/null | head -25 | tee -a "$LOG" || echo "numactl 不可用" | tee -a "$LOG"
echo "==> 当前 mpiexec: $(grep -E 'AMSS_MPIEXEC|mpiexec' run.sh | head -2 | tr '\n' ' ')" | tee -a "$LOG"
echo "==> cpuset: $(grep Cpus_allowed_list /proc/self/status 2>/dev/null)" | tee -a "$LOG"
echo "==> taskset: $(taskset -pc $$ 2>/dev/null)" | tee -a "$LOG"
echo | tee -a "$LOG"

# 计算节点无 build( .gitignore 不同步), 先编译一次
bash compile.sh 2>&1 | tail -2 | tee -a "$LOG"

for R in "$@"; do
  case "$R" in
    *:*) RANK="${R%%:*}"; OMP="${R##*:}" ;;
    *)   RANK="$R";       OMP=1 ;;
  esac
  sed -i "s/^MPI_processes    = .*/MPI_processes    = $RANK/" "$INPUT"
  sed -i "s/^OMP_threads      = .*/OMP_threads      = $OMP/" "$INPUT"
  echo "########## MPI=$RANK OMP=$OMP ##########" | tee -a "$LOG"
  ./run.sh --twop-cache 2>&1 | tee -a "$LOG" | grep -aE "Running ./ABE|Timestep #|Total Evolve Time|This Program Cost"
done

echo "==> 恢复后 Input.py: $(grep '^MPI_processes' "$INPUT")" | tee -a "$LOG"
echo "DONE (完整日志: $LOG)" | tee -a "$LOG"
