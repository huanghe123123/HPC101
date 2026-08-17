#!/bin/bash
# profile.sh — perf stat + perf record on the AMSS-NCKU CPU baseline.
#
# Strategy: profile the full `./run.sh --twop-cache` so perf follows the whole
# Python→TwoPuncture→mpirun→ABE tree (perf_event inheritance is on by default).
# - Run 1 (perf stat): also generates the TwoPuncture cache.
# - Run 2 (perf record): cache hit → skips TwoPuncture → clean ABE profile.
# All outputs go to stdout (captured to j<id>.out in home by the scheduler) AND
# to $PROF_DIR. perf.data is written to /workspace scratch (writable on compute
# nodes; /home may be read-only there) and best-effort copied to home.
set -uo pipefail
export LC_ALL=C

LAB=/home/h3250104945/HPC101/src/lab4
[[ -d "$LAB" ]] || LAB="$(pwd)"
cd "$LAB" || { echo "FATAL: cannot cd to $LAB"; exit 1; }
echo "==> LAB root : $LAB"
echo "==> pwd      : $(pwd)"

# Writable scratch for perf.data and cache (platform may mount /home RO on
# compute nodes; /workspace is the writable scratch it provides).
SCRATCH=/workspace/lab4
mkdir -p "$SCRATCH/twopuncture_cache" 2>/dev/null || SCRATCH="$LAB"
export AMSS_CACHE_DIR="$SCRATCH/twopuncture_cache"
# AMSS_BUILD_DIR / AMSS_OUTPUT_ROOT are injected by the platform → /workspace.

# Persistent output dir (home). Best-effort; reports also go to stdout.
PROF_DIR=/home/h3250104945/HPC101/src/lab4/perf_profiles
if ! mkdir -p "$PROF_DIR" 2>/dev/null || ! : >"$PROF_DIR/.w" 2>/dev/null; then
  echo "==> home perf dir not writable, using scratch"
  PROF_DIR="$SCRATCH/perf_profiles"
  mkdir -p "$PROF_DIR"
fi
rm -f "$PROF_DIR/.w" 2>/dev/null
echo "==> PROF_DIR : $PROF_DIR"
echo "==> SCRATCH  : $SCRATCH"
echo "==> perf     : $(perf --version 2>&1)"
echo "==> paranoid : $(cat /proc/sys/kernel/perf_event_paranoid 2>&1)"
echo "==> whoami   : $(whoami)"
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [1/4] compile with -g (symbols for profile attribution)"
echo "############################################################"
./compile.sh -DAMSS_OPT='-O3 -g' 2>&1 | tail -12
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [2/4] perf stat  (4-step run; also builds twop cache)"
echo "############################################################"
perf stat -d -d -d -- ./run.sh --twop-cache 2>&1 | tee "$PROF_DIR/perf_stat.txt"
echo "==> perf stat saved: $PROF_DIR/perf_stat.txt"
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [3/4] perf record (4-step, twop cache hit, call-graph=dwarf)"
echo "############################################################"
PERFDATA="$SCRATCH/perf.data"
rm -f "$PERFDATA"
perf record -F 49 --call-graph dwarf -o "$PERFDATA" -- ./run.sh --twop-cache 2>&1 \
  | tee "$PROF_DIR/perf_record.log"
echo "==> perf.data: $(ls -lh "$PERFDATA" 2>&1)"
echo

#---------------------------------------------------------------
echo "############################################################"
echo "# [4/4] perf report  (rendered to stdout so it lands in j<id>.out)"
echo "############################################################"
echo "===================== REPORT: SELF TIME ====================="
perf report -i "$PERFDATA" --stdio --no-children --percent-limit 0.3 2>&1 | sed -n '1,140p' \
  | tee "$PROF_DIR/perf_report_self.txt"

echo "================= REPORT: CUMULATIVE (call paths) ================="
perf report -i "$PERFDATA" --stdio --percent-limit 0.5 2>&1 | sed -n '1,100p' \
  | tee "$PROF_DIR/perf_report_total.txt"

echo "================= REPORT: BY COMM (per-process/rank) ================="
perf report -i "$PERFDATA" --stdio --sort=comm 2>&1 | sed -n '1,40p'

echo "================= REPORT: TOP SYMBOLS (self, sorted) ================="
perf report -i "$PERFDATA" --stdio --no-children --sort=symbol --percent-limit 0.3 2>&1 \
  | sed -n '1,60p' | tee "$PROF_DIR/perf_report_symbols.txt"

# Best-effort: persist perf.data to home for later flame-graph / GUI work.
if cp "$PERFDATA" "$PROF_DIR/perf.data" 2>/dev/null; then
  echo "==> perf.data copied to $PROF_DIR/perf.data ($(ls -lh "$PROF_DIR/perf.data" | awk '{print $5}'))"
else
  echo "==> perf.data left in ephemeral $PERFDATA (home not writable); reports are in stdout + $PROF_DIR/*.txt"
fi

echo
echo "############################################################"
echo "# DONE."
echo "############################################################"
