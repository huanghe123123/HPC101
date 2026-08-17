#!/bin/bash
# profile3.sh — compute-node perf record retry with MINIMAL mmap (-m 1).
#
# Root cause of prior failures: compute node forbids raising
# perf_event_mlock_kb (RO /proc/sys in container), and -m 1024/256 exceeded
# the default mlock budget → "Permission error mapping pages".
# Fix: -m 1 = 1 page (8 KiB) per CPU mmap ring, far under the ~516 KiB default.
# We keep the full chain (NO twop cache) so TwoPuncture stays visible.
set -uo pipefail
export LC_ALL=C
cd /home/h3250104945/HPC101/src/lab4 || exit 1
echo "==> pwd=$(pwd)  paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)  mlock_kb=$(cat /proc/sys/kernel/perf_event_mlock_kb)"

SCRATCH=/workspace/lab4; mkdir -p "$SCRATCH" 2>/dev/null || SCRATCH="$(pwd)"
PROF_DIR=perf_profiles; mkdir -p "$PROF_DIR"

# Already compiled with -g in profile2 run, but rebuild to be safe.
./compile.sh -DAMSS_OPT='-O3 -g -fno-omit-frame-pointer' 2>&1 | tail -4
echo

PERFDATA="$SCRATCH/perf3.data"; rm -f "$PERFDATA"
echo "############################################################"
echo "# perf record -m 1 (minimal mmap) — full chain, no cache"
echo "############################################################"
perf record -F 99 -m 1 --call-graph dwarf -o "$PERFDATA" -- ./run.sh 2>&1 | tee "$PROF_DIR/perf3_record.log" | tail -6
echo "==> perf3.data: $(ls -lh "$PERFDATA" 2>&1)"

if [[ -s "$PERFDATA" ]]; then
  echo "############### REPORT: SELF TIME ###############"
  perf report -i "$PERFDATA" --stdio --no-children --percent-limit 0.1 2>&1 | sed -n '1,180p' | tee "$PROF_DIR/perf3_self.txt"
  echo "############### REPORT: TOP SYMBOLS ###############"
  perf report -i "$PERFDATA" --stdio --no-children --sort=symbol --percent-limit 0.1 2>&1 | sed -n '1,80p' | tee "$PROF_DIR/perf3_symbols.txt"
  echo "############### REPORT: BY COMM ###############"
  perf report -i "$PERFDATA" --stdio --sort=comm 2>&1 | sed -n '1,50p' | tee "$PROF_DIR/perf3_comm.txt"
  cp "$PERFDATA" "$PROF_DIR/perf3.data" 2>/dev/null && echo "==> saved perf3.data"
else
  echo "!!! -m 1 still empty — perf record truly blocked on compute node"
fi
echo "DONE"
