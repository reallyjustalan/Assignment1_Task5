#!/usr/bin/env bash
# Run from a terminal (or tmux) and screenshot the command plus completion output.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh
mode="${1:-}"
case "$mode" in
  annotated) resume=() ;;
  unannotated)
    [[ -s results/annotated/run_name.txt ]] || {
      echo 'Complete the annotated scenario first.' >&2; exit 1;
    }
    resume=(-resume "$(< results/annotated/run_name.txt)")
    ;;
  *) echo "Usage: bash run_scenario.sh annotated|unannotated" >&2; exit 1 ;;
esac
mkdir -p "results/$mode"
# Do not silently overwrite evidence from an earlier attempt.
if [[ -e "results/$mode/report.html" || -e "results/$mode/terminal.log" ]]; then
    echo "Existing evidence in results/$mode; inspect it before rerunning." >&2
    exit 1
fi
command=(nextflow run main.nf --mode "$mode" --outdir "$PWD/results"
         -work-dir "$PWD/work" -ansi-log false
         -with-report "results/$mode/report.html"
         -with-trace "results/$mode/trace.tsv"
         -with-timeline "results/$mode/timeline.html" "${resume[@]}")
{ printf 'COMMAND: '; printf '%q ' "${command[@]}"; printf '\n'; } | tee "results/$mode/terminal.log"
"${command[@]}" 2>&1 | tee -a "results/$mode/terminal.log"
# Store the exact run name so Scenario 2 cannot accidentally resume a smoke test.
nextflow log -q | tail -n 1 > "results/$mode/run_name.txt"
