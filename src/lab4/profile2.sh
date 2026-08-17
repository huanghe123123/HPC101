#!/bin/bash
# profile2.sh — perf record (fixed mmap) on AMSS-NCKU CPU baseline.
#
# Fix vs profile.sh: perf record died on "Permission error mapping pages".
# Causes: default per-CPU mmap too big for perf_event_mlock_kb. Fixes:
#   1. raise kernel.perf_event_mlock_kb (we run as root in the job container)
#   2. cap -m (mmap pages per cpu) explicitly
#   3. build TwoPuncture cache in PERSISTENT home, then perf-record a cache-hit
#      run so the profile is ABE evolution, not the TwoPuncture bicgstab solver.
# Also add -fno-omit-frame-pointer so FP callgraph works as a dwarf fallback.
set -uo pipefail
export LC_ALL=C

LAB=/home/h3250104945/HPC101/src/lab4
cd "$LAB" || { echo "FATAL: cd $LAB"; exit 1; }
echo "==> pwd=$(pwd)"

PROF_DIR=$LAB/perf_profiles
mkdir -p "$PROF_DIR"
# Persistent TwoPuncture cache in HOME (survives across jobs); /workspace is ephemeral.
CACHE_DIR=$LAB/twopuncture_cache
mkdir -p "$CACHE_DIR"
export AMSS_CACHE_DIR="$CACHE_DIR"

SCRATCH=/workspace/lab4
mkdir -p "$SCRATCH" 2>/dev/null || SCRATCH="$LAB"

echo "==> perf=$(perf --version 2>&1)"
echo "==> paranoid=$(cat /proc/sys/kernel/perf_event_paranoid 2>&1)"
echo "==> mlock_kb(before)=$(cat /proc/sys/kernel/perf_event_mlock_kb 2>&1)"

# (1) Raise mlock budget as root (best-effort; ignore if RO).
sysctl -w kernel.perf_event_mlock_kb=1048576 2>/dev/null || \
  echo "(could not raise mlock_kb; relying on -m cap)"
echo "==> mlock_kb(after)=$(cat /proc/sys/kernel/perf_event_mlock_kb 2>&1)"
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [1/4] compile -O3 -g -fno-omit-frame-pointer"
echo "############################################################"
./compile.sh -DAMSS_OPT='-O3 -g -fno-omit-frame-pointer' 2>&1 | tail -6
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [2/4] build TwoPuncture cache (plain run, no perf)"
echo "############################################################"
echo "    cache dir: $CACHE_DIR"
./run.sh --twop-cache 2>&1 | tail -8
echo "==> cache contents:"; ls -lh "$CACHE_DIR" 2>&1 | tail -6
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [3/4] perf record  (cache-hit run → ABE evolution profile)"
echo "############################################################"
PERFDATA="$SCRATCH/perf.data"
rm -f "$PERFDATA"
# -m 1024 = 1024 pages (4 MiB) per CPU; capped so mlock stays in budget.
# -F 99 Hz; --call-graph dwarf for stack traces.
perf record -F 99 -m 1024 --call-graph dwarf -o "$PERFDATA" -- \
    ./run.sh --twop-cache 2>&1 | tee "$PROF_DIR/perf_record2.log"
echo "==> perf.data: $(ls -lh "$PERFDATA" 2>&1)"
echo

# If record still produced nothing, try a minimal fallback.
if [[ ! -s "$PERFDATA" ]]; then
  echo "!!! perf.data empty — retry with smaller mmap + FP callgraph"
  rm -f "$PERFDATA"
  perf record -F 99 -m 256 --call-graph fp -o "$PERFDATA" -- \
      ./run.sh --twop-cache 2>&1 | tail -8
  echo "==> perf.data(retry): $(ls -lh "$PERFDATA" 2>&1)"
fi
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [4/4] perf report  (function-level hotspots)"
echo "############################################################"
if [[ -s "$PERFDATA" ]]; then
  echo "===================== SELF TIME (own work) ====================="
  perf report -i "$PERFDATA" --stdio --no-children --percent-limit 0.2 2>&1 | sed -n '1,160p' \
    | tee "$PROF_DIR/perf2_report_self.txt"

  echo "================= CUMULATIVE (call paths) ================="
  perf report -i "$PERFDATA" --stdio --percent-limit 0.5 2>&1 | sed -n '1,120p' \
    | tee "$PROF_DIR/perf2_report_total.txt"

  echo "================= TOP SYMBOLS (self) ================="
  perf report -i "$PERFDATA" --stdio --no-children --sort=symbol --percent-limit 0.2 2>&1 \
    | sed -n '1,80p' | tee "$PROF_DIR/perf2_report_symbols.txt"

  echo "================= BY COMM (process/rank split) ================="
  perf report -i "$PERFDATA" --stdio --sort=comm 2>&1 | sed -n '1,50p'

  cp "$PERFDATA" "$PROF_DIR/perf2.data" 2>/dev/null && \
    echo "==> perf2.data saved ($(ls -lh "$PROF_DIR/perf2.data"|awk '{print $5}'))"
else
  echo "!!! perf.data still empty — perf record unusable on this node"
fi

echo
echo "############################################################"
echo "# DONE."
echo "############################################################"
