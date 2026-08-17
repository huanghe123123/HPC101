#!/bin/bash
# profile4.sh — compute-node perf record, tuned to STOP sample loss.
#
# Diagnosis of profile3 failure: 542375 samples taken but "lost 100.00%".
# Root cause: -m 1 (1 page/8KiB per CPU ring) is far too small → the ring
# overflows between samples and the kernel drops them. Need a LARGER ring, but
# larger mmap hit the mlock permission error earlier. The key realization:
# the mlock limit (516 KiB) applies to the TOTAL across all perf users. With
# 60 ranks the per-CPU budget is tiny. Two fixes combined:
#   (a) shrink concurrency: profile a SINGLE-rank run (MPI_processes=1) so
#       only a few CPUs are active → mmap ring per active CPU fits in budget.
#   (b) use a moderate -m (e.g. 64 pages = 256 KiB) per CPU — fits under 516.
# Single-rank is still representative for finding compute hotspots (BSSN RHS,
# diff kernels); it just removes MPI-wait noise, which perf stat already
# characterized. NO twop cache so TwoPuncture stays visible.
set -uo pipefail
export LC_ALL=C
cd /home/h3250104945/HPC101/src/lab4 || exit 1
echo "==> pwd=$(pwd)  paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)  mlock_kb=$(cat /proc/sys/kernel/perf_event_mlock_kb)"

SCRATCH=/workspace/lab4; mkdir -p "$SCRATCH" 2>/dev/null || SCRATCH="$(pwd)"
PROF_DIR=perf_profiles; mkdir -p "$PROF_DIR"

# Temporarily set MPI to 1 for this profile run only.
python3 - <<'PY'
import re
p="AMSS_NCKU_Input.py"
s=open(p).read()
s2=re.sub(r'^MPI_processes\s*=\s*\d+', 'MPI_processes    = 1', s, count=1, flags=re.M)
if s2!=s: open(p,"w").write(s2); print("set MPI_processes=1 for profile4")
else: print("MPI_processes unchanged")
PY

./compile.sh -DAMSS_OPT='-O3 -g -fno-omit-frame-pointer' 2>&1 | tail -3
echo

PERFDATA="$SCRATCH/perf4.data"; rm -f "$PERFDATA"
echo "############################################################"
echo "# perf record -m 64 (256KiB/CPU), MPI=1, full chain no cache"
echo "############################################################"
perf record -F 99 -m 64 --call-graph dwarf -o "$PERFDATA" -- ./run.sh 2>&1 | tee "$PROF_DIR/perf4_record.log" | tail -8
echo "==> perf4.data: $(ls -lh "$PERFDATA" 2>&1)"

if [[ -s "$PERFDATA" ]]; then
  echo "############### REPORT: SELF TIME ###############"
  perf report -i "$PERFDATA" --stdio --no-children --percent-limit 0.1 2>&1 | sed -n '1,200p' | tee "$PROF_DIR/perf4_self.txt"
  echo "############### REPORT: TOP SYMBOLS ###############"
  perf report -i "$PERFDATA" --stdio --no-children --sort=symbol --percent-limit 0.1 2>&1 | sed -n '1,90p' | tee "$PROF_DIR/perf4_symbols.txt"
  echo "############### REPORT: CUMULATIVE ###############"
  perf report -i "$PERFDATA" --stdio --percent-limit 0.5 2>&1 | sed -n '1,120p' | tee "$PROF_DIR/perf4_total.txt"
  cp "$PERFDATA" "$PROF_DIR/perf4.data" 2>/dev/null && echo "==> saved perf4.data"
else
  echo "!!! perf4 empty"
fi
echo "DONE"
