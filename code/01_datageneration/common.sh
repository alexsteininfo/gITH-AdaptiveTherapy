# Shared functions of the run_*.sh scripts (sourced, not executed).
# Env: WORKERS (default 12), JULIA (default julia), TEST=1 for the smoke test.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
JULIA=${JULIA:-julia}
WORKERS=${WORKERS:-12}
TEST=${TEST:-0}
export TEST
mkdir -p data/raw/logs
LOG="data/raw/logs/$(basename "$0" .sh)-$(date +%Y-%m-%d_%H%M).log"
exec > >(tee -a "$LOG") 2>&1
echo "# $(basename "$0") on $(hostname), $(date), WORKERS=$WORKERS, TEST=$TEST, julia=$JULIA"

# processed block name (the smoke test writes to 00_test_<block>)
pblock() { if [ "$TEST" = "1" ]; then echo "00_test_$1"; else echo "$1"; fi; }

# simulate one block (extra arguments: cohort names), then collect and compute statistics
simulate_and_process() {
    local block=$1; shift
    "$JULIA" -p "$WORKERS" --project=. "code/01_datageneration/$block.jl" "$@"
    "$JULIA" -p "$WORKERS" --project=. code/02_processing/01_collect.jl "$(pblock "$block")"
    "$JULIA" -p "$WORKERS" --project=. code/02_processing/02_statistics.jl "$(pblock "$block")"
}

# summary script of code/02_processing for some blocks (skipped if the script does not exist yet)
summarise() {
    local script=$1; shift
    local args=()
    for b in "$@"; do args+=("$(pblock "$b")"); done
    if [ -f "code/02_processing/$script" ]; then
        "$JULIA" -p "$WORKERS" --project=. "code/02_processing/$script" "${args[@]}"
    else
        echo "skipping $script (not present)"
    fi
}
