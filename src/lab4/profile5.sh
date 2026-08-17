#!/bin/bash
# profile5.sh — 30-rank perf record 探测 + 多 rank 函数级热点
#
# 目的: 之前 30 rank record 都失败(默认 -m 超预算, -m 1 丢样本)。
# 现在页大小=4KB, mlock 预算=516KB (全局总额, 不是每核一份)。
# 30 核时: 缓冲总量 = 30 × (-m × 4KB)
#   -m 1  → 30×4KB  = 120KB  < 516 (不超, 但太小 → 之前丢 100%)
#   -m 2  → 30×8KB  = 240KB  < 516 (不超, 稍大)
#   -m 4  → 30×16KB = 480KB  < 516 (不超, 4倍于 -m1, 可能够用)
#   -m 8  → 30×32KB = 960KB  > 516 (超 → Permission error)
# 所以唯一可能成功的中间值是 -m 4 (480KB)。本脚本试 -m 4, 失败则退 -m 2。
#
# 用 --twop-cache 跳过 TwoPuncture (初值不变), 让 profile 聚焦 ABE 演化 + MPI。
# 步数 4 (临时)。MPI=30。
set -uo pipefail
export LC_ALL=C
cd /home/h3250104945/HPC101/src/lab4 || exit 1
echo "==> pwd=$(pwd) paranoid=$(cat /proc/sys/kernel/perf_event_paranoid) mlock_kb=$(cat /proc/sys/kernel/perf_event_mlock_kb) pagesize=$(getconf PAGESIZE)"

# 确保 TwoPuncture cache 存在 (home 跨 job 保留)
export AMSS_CACHE_DIR=/home/h3250104945/HPC101/src/lab4/twopuncture_cache
echo "==> cache: $(ls $AMSS_CACHE_DIR/*/ 2>&1 | head -3)"

SCRATCH=/workspace/lab4; mkdir -p "$SCRATCH" 2>/dev/null || SCRATCH="$(pwd)"
PROF_DIR=perf_profiles; mkdir -p "$PROF_DIR"

# 编译带 -g
./compile.sh -DAMSS_OPT='-O3 -g -fno-omit-frame-pointer' 2>&1 | tail -3
echo

# 尝试 -m 4
PERFDATA="$SCRATCH/perf5_m4.data"; rm -f "$PERFDATA"
echo "############################################################"
echo "# 尝试1: 30 rank + -m 4 (480KB, 不超预算)"
echo "############################################################"
perf record -F 99 -m 4 --call-graph dwarf -o "$PERFDATA" -- ./run.sh --twop-cache 2>&1 | tee "$PROF_DIR/perf5_m4.log" | tail -10
echo "==> perf5_m4.data: $(ls -lh "$PERFDATA" 2>&1)"

if [[ -s "$PERFDATA" ]]; then
  DATA="$PERFDATA"; OUT=perf5_m4
else
  echo "!!! -m 4 失败, 退回 -m 2"
  PERFDATA="$SCRATCH/perf5_m2.data"; rm -f "$PERFDATA"
  perf record -F 99 -m 2 --call-graph dwarf -o "$PERFDATA" -- ./run.sh --twop-cache 2>&1 | tee "$PROF_DIR/perf5_m2.log" | tail -10
  echo "==> perf5_m2.data: $(ls -lh "$PERFDATA" 2>&1)"
  DATA="$PERFDATA"; OUT=perf5_m2
fi
echo

if [[ -s "$DATA" ]]; then
  echo "############################################################"
  echo "# 成功! 生成 report (原始, 计算节点可能不解析符号)"
  echo "############################################################"
  perf report -i "$DATA" --stdio --no-children --sort=symbol --percent-limit 0.3 2>&1 | tee "$PROF_DIR/${OUT}_symbols_raw.txt" | tail -40
  cp "$DATA" "$PROF_DIR/${OUT}.data" 2>/dev/null && echo "==> 已存 $PROF_DIR/${OUT}.data (回 devpod 用 grep -vE 'unwind|get_proc_name' 重解析)"
else
  echo "!!! 30 rank record 彻底失败 (-m 2/4 都不行), 确认 30 rank 死路"
fi
echo "DONE"
