# 01_datageneration — simulation specifications

Simulations for `.docs/theory/PLAN_v1.md`, run with **ClonalEvolutionTx.jl**
(`/Users/alexanderstein/Documents/GitHub/ClonalEvolutionTx.jl`).

**Design.** The simulations store **distributions**, not statistics: at selected time points the
site frequency spectrum (SFS), the single-cell mutational burden (scMB) distribution and the clone
sizes, for the entire tree and for subsamples of 100, 1000 and 10⁴ cells; plus simple population
counts over time. All statistics (Hill numbers, entropies, power-law exponents with any $k_\text{min}$,
detection thresholds $f_\text{min}$, …) are computed afterwards in `code/02_processing/` and can be
recomputed freely. Trees are not stored, except for 10 example tumours per cohort.

## Run

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'     # once (also on the cluster)
bash code/01_datageneration/run_part1.sh                 # Part 1: dynamics
bash code/01_datageneration/run_part2.sh                 # Part 2: trade-off + sensitivity
bash code/01_datageneration/run_part3.sh                 # Part 3: virtual trial
bash code/01_datageneration/run_supp.sh                  # supplement
WORKERS=32 JULIA=/path/to/julia bash code/01_datageneration/run_part2.sh   # cluster
TEST=1 bash code/01_datageneration/run_part1.sh          # smoke test (N_det = 10^4, 3 tumours)
```

Each bash script simulates its blocks (`WORKERS` worker processes, default 12) and then runs the
processing of those blocks (`code/02_processing/`), logging to `data/raw/logs/`. Runs are
resumable (finished tumours are skipped); cohort names can be passed to run a subset. Figures are
made interactively from `data/processed/` (`code/03_figures/`).

## Model (all cohorts unless stated)

| | Value |
|---|---|
| birth, death | $b=1$, $d_0=0.5$ (1 time unit = 1 cell cycle ≈ 4 days) |
| competition | death $+\,lN$ with $l=(b-d_0)/K$, $K=N_\text{det}/f_K$ |
| drug | death $+\,\alpha\, d_S\, m$, $d_S=4$, dose $m\in[0,1]$ |
| neutral mutations | Poisson($\mu=1$) per daughter cell |
| resistance mutations | probability $\nu$ per daughter cell; new clone with resistance factor $\alpha$ |
| detection | $N_\text{det}=10^6$, grown from one cell (restart if extinct) |
| progression | $N\ge N_\text{prog}=(N_\text{det}+K)/2$ (halfway to $K$) |
| follow-up | $t_\text{max}=250$ after detection (≈ 1000 days); later progression is censored |
| defaults | $f_K=N_\text{det}/K=0.8$, $\nu=10^{-4}$ |

Resistance classes: **full** $\alpha=0$; **partial** $\alpha=0.006$ (still progresses under MTD near
$K$, slower); **weak** $\alpha=0.8$ (curable by MTD).

## Schedules (from detection)

| arm | rule | levels $L$ (fraction of $N_\text{det}$) |
|---|---|---|
| `single` | full dose until $N\le0.1N_\text{det}$, then off until $N\ge N_\text{det}$ | – |
| `mtd` | full dose until progression or extinction | – |
| `hol_L` | **holidays**: full dose until $N\le L\,N_\text{det}$, off until $N\ge N_\text{det}$, repeat | 0.5, 0.6, 0.7, 0.8, 0.9 |
| `mod_L` | **modulation**: full dose until $N\le L\,N_\text{det}$, then the dose that keeps $E[\dot N]=0$ | 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 |
| `met` | metronomic: full dose for 1 of every 3 time units | – |

`hol` = `hol_0.5` and `mod` = `mod_0.8` are the representative settings used outside the scans.
The scans are the actual object: which containment rule works best, not a copy of one trial.

## All planned simulations

| # | block | cohort | $\nu$ | $\alpha$ | $f_K$ | tumours | arms | $t_\text{max}$ | purpose | CPU h |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `01_dynamics` | `neutral` | 0 | – | 0.8 | 500 | single, hol, mod | 40 | Figs 2–3: neutral dynamics | 8 |
| 2 | | `full` | $10^{-4}$ | 0 | 0.8 | 500 | single, mtd, hol, mod | 250 | Figs 1, 2, 4; within-cohort prediction (Fig 6) | 10 |
| 3 | | `partial` | $10^{-4}$ | 0.006 | 0.8 | 500 | single, mtd, hol, mod | 250 | Fig 4 | 12 |
| 4 | | `weak` | $10^{-4}$ | 0.8 | 0.8 | 500 | single, mtd, hol, mod | 60 | Fig 2 | 11 |
| 5 | `02_tradeoff` | `full` (tumours 1–300 of #2) | $10^{-4}$ | 0 | 0.8 | 300 | mtd, hol_L ×5, mod_L ×6 | 250 | Fig 5: trade-off curve | 28 |
| 6 | | `partial` (tumours 1–300 of #3) | $10^{-4}$ | 0.006 | 0.8 | 300 | same | 250 | Fig 5 | 29 |
| 7 | `03_sensitivity` | `fK0.5` | $10^{-4}$ | 0 | 0.5 | 200 | same | 250 | Fig 5: competition | 11 |
| 8 | | `fK0.7` | $10^{-4}$ | 0 | 0.7 | 200 | same | 250 | Fig 5 | 14 |
| 9 | | `fK0.9` | $10^{-4}$ | 0 | 0.9 | 200 | same | 250 | Fig 5 | 28 |
| 10 | | `nu1e-5` | $10^{-5}$ | 0 | 0.8 | 200 | same | 250 | Fig 5: rare resistance | 33 |
| 11 | | `nu1e-3` | $10^{-3}$ | 0 | 0.8 | 200 | same | 250 | Fig 5: frequent resistance | 8 |
| 12 | `04_trial` | `trial` | log-uniform $10^{-5}$–$10^{-3}$ | 0–0.006 | 0.5–0.9 | 500 (Latin hypercube) | mtd, hol, mod | 250 | Fig 6: does diversity predict course and AT benefit? | 9 |
| 13 | `05_metronomic` | `neutral`, `full` (same tumours as #1, #2) | | | 0.8 | 500 each | met | 40 / 250 | Fig 1, Supp 2 | 7 |
| 14 | `06_cost` | `cost0`, `cost0.05` | $10^{-4}$ | Beta(1,4), min rule; cost $\Delta s=-c(1-\alpha)$ | 0.8 | 300 each | single, mtd, hol, mod, mod_1.0 | 250 | Supp 3: cost of resistance | 30 |
| 15 | `07_validation` | `sfs_pb`, `sfs_bd`, `ld` | 0 / 0 / $10^{-3}$ | – / – / 0 | no competition ($d_0=0$ for `sfs_pb`) | 200 each, $N_\text{det}=10^5$ | growth only | – | Supp 1: checks vs theory | < 1 |

Rows with the same parameters share tumours (same seeds): #5 ⊂ #2, #6 ⊂ #3, #13 = #1/#2. All arms
of a tumour start from the same detection-time population with the same random seed afterwards
(paired comparison). **Total ≈ 240 CPU hours** (≈ 20 h on 12 local workers, ≈ 10–15 h on 32
cluster cores); disk ≈ 20–40 GB. Rows 5–11 (`run_part2.sh`, ≈ 150 CPU h) are the largest part.

## What is stored per tumour (`data/raw/<block>/<cohort>/tNNNN.jld2`)

| record | when | content |
|---|---|---|
| outcome | per arm | progressed / extinct / censored, time to progression, cumulative dose, cycles, block log |
| trajectory | every 0.25 | $N$, dose $m$, sensitive and resistant cells |
| counts | every 1 | number of resistant clones; resistant lineages (resistance arisen in a sensitive cell) before / after detection; largest resistant clone, mean burden |
| distributions | key time points | entire tree: SFS (sparse), scMB distribution, clone sizes with $\alpha$ and origin time, total burden $B$ split into pre-existing / new; subsamples $n=100,1000,10^4$: SFS, scMB, clone sizes |

Key time points: growth at $0.1N_\text{det}$ (block 01); detection; end of each of the first 12
treatment blocks (6 holiday cycles); every 10 time units; end. Scan arms and blocks 05–06:
detection and end only. Trial: detection, end of the first full-dose phase, end.
