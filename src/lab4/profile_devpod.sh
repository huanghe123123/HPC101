#!/bin/bash
# profile_devpod.sh — full-chain perf on devpod (NO twop cache).
# Covers Python → TwoPuncture → mpirun → ABE so the profile reflects the
# real end-to-end breakdown, per the requirement that TwoPuncture must be
# visible (it may itself be a hotspot).
set -uo pipefail
export LC_ALL=C
cd /home/h3250104945/HPC101/src/lab4 || exit 1

PROF_DIR=perf_profiles
mkdir -p "$PROF_DIR"
echo "==> perf=$(perf --version 2>&1)  paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)  mlock_kb=$(cat /proc/sys/kernel/perf_event_mlock_kb)"
echo "==> Final_Evolution_Time steps (4) — NO twop cache, full chain"
echo

echo "############################################################"
echo "# [1/2] perf stat -dddd (full chain, no cache)"
echo "############################################################"
perf stat -d -d -d -o "$PROF_DIR/devpod_perf_stat.txt" -- ./run.sh 2>&1 | tail -40
echo "==> stat saved"
echo

echo "############################################################"
echo "# [2/2] perf record (full chain, no cache, dwarf call-graph)"
echo "############################################################"
rm -f "$PROF_DIR/devpod_perf.data"
# -m 512 per-cpu mmap (within devpod mlock budget); -F 99 Hz; dwarf stacks.
perf record -F 99 -m 512 --call-graph dwarf -o "$PROF_DIR/devpod_perf.data" -- ./run.sh 2>&1 | tee "$PROF_DIR/devpod_record.log" | tail -8
echo "==> perf.data: $(ls -lh "$PROF_DIR/devpod_perf.data" 2>&1)"
echo

if [[ -s "$PROF_DIR/devpod_perf.data" ]]; then
  echo "############################################################"
  echo "# REPORT: SELF TIME (function's own work)"
  echo "############################################################"
  perf report -i "$PROF_DIR/devpod_perf.data" --stdio --no-children --percent-limit 0.1 2>&1 | sed -n '1,180p' | tee "$PROF_DIR/devpod_report_self.txt"

  echo "############################################################"
  echo "# REPORT: TOP SYMBOLS (self, sorted)"
  echo "############################################################"
  perf report -i "$PROF_DIR/devpod_perf.data" --stdio --no-children --sort=symbol --percent-limit 0.1 2>&1 | sed -n '1,80p' | tee "$PROF_DIR/devpod_report_symbols.txt"

  echo "############################################################"
  echo "# REPORT: CUMULATIVE (call paths)"
  echo "############################################################"
  perf report -i "$PROF_DIR/devpod_perf.data" --stdio --percent-limit 0.5 2>&1 | sed -n '1,140p' | tee "$PROF_DIR/devpod_report_total.txt"

  echo "############################################################"
  echo "# REPORT: BY COMM (TwoPuncture vs ABE ranks vs python)"
  echo "############################################################"
  perf report -i "$PROF_DIR/devpod_perf.data" --stdio --sort=comm 2>&1 | sed -n '1,50p' | tee "$PROF_DIR/devpod_report_comm.txt"
else
  echo "!!! perf.data empty on devpod too"
fi
echo
echo "############################################################"
echo "# DONE — all reports in $PROF_DIR/"
echo "############################################################"
ls -lh "$PROF_DIR/"
