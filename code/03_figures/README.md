# 03_figures — figure panels

One script per figure; each reads only tables in `data/processed/` and writes one PNG per panel
(300 dpi, 45 × 45 mm unless stated, no panel letters or titles) to `figures/<figure>/`. Panels
are assembled and lettered in Illustrator. Run interactively, e.g.

```bash
julia --project=. code/03_figures/fig_2/fig_2.jl
```

or from a Julia session (`include("code/03_figures/fig_2/fig_2.jl")`) to iterate on one panel.
Shared style — theme, colours, labels, panel size, `readproc`, `savepanel`, summaries over
tumours — is in `style.jl`.

**Test data.** `FIG_DATA=test` makes the scripts read the smoke-test tables
(`data/processed/00_test_<block>`, from `TEST=1 bash code/01_datageneration/run_*.sh`) and write
to `figures/00_test/`. `julia --project=. code/03_figures/smoke.jl [regex]` runs every script on
the test data (making missing test tables first) and checks that all panels are written.

| figure | script | tables (`data/processed/…`) | panels |
|---|---|---|---|
| 1 Model and schedules | `fig_1/fig_1.jl` | `01_dynamics/example_traj.csv`, `example_clones.csv` | `schematic.ai/.png` (a–c, copied from BirthDeathMutation_AT.jl), `d` MTD, `e` metronomic, `f` holidays, `g` modulation, `legend` |
| 2 One treatment round | `fig_2/fig_2.jl` | `01_dynamics/sfs_mean.csv`, `snapshots.csv` | single round, responders only, cohorts none/weak/partial/full: top row after treatment (0.1N), bottom row after regrowth (N): `a`/`d` mean SFS vs untreated, `b`/`e` passenger ᵠD treated/untreated, `c`/`f` resistant-clone ᵠD treated/untreated, `legend` |
| 3 Holidays and modulation (neutral) | `fig_3/fig_3.jl` | `01_dynamics/snapshots.csv`, `sfs_mean.csv` | `a` S per cycle, `b` γ̂ per cycle, `c` burden B pre-existing/new under modulation with 2μN ln(1+t), `d` SFS under modulation |
| 4 Composition under the schedules | `fig_4/fig_4.jl` | `01_dynamics/counts.csv.gz`, `snapshots.csv`, `km.csv` | per cohort (`_full`, `_partial`): `a` resistant clones, `b` pre-existing/new, `c` resistant ²D, `d` resistant fraction, `e` Kaplan–Meier; `legend` |
| 5 Trade-off | `fig_5/fig_5.jl` | `02_tradeoff/tradeoff.csv`, `03_sensitivity/tradeoff.csv` | `a`/`b` TTP vs resistant clones (full/partial), `c` TTP vs B, `d_hol`/`d_mod` gain vs level by N/K, `e_hol`/`e_mod` by ν, `legend` |
| 6 Virtual trial | `fig_6/fig_6.jl` | `04_trial/trial.csv`, `trial_cor.csv`, `trial_policy.csv` | `a` correlation heatmap (110 × 95 mm), `b` gain vs best predictor, `c` \|ρ\| vs sample size, `d` policies |
| S1 Validation | `supp_1/supp_1.jl` | `07_validation/validation_sfs.csv`, `validation_ld.csv`, `snapshots.csv` | `a`/`b` SFS vs theory (pure birth / birth–death), `c` resistant clones, `d` resistant cells |
| S2 Metronomic | `supp_2/supp_2.jl` | `05_metronomic/*`, `01_dynamics/*` (counts, km, outcomes) | `a` mean burden, `b` resistant clones, `c` Kaplan–Meier, `d` TTP metronomic vs MTD |
| S3 Cost of resistance | `supp_3/supp_3.jl` | `06_cost/counts.csv.gz`, `km.csv` | `a_cost0`, `a_cost0.02` resistant clones; `b_*` Kaplan–Meier; `legend` |
| S4 Hill definitions | `supp_4/supp_4.jl` | `01_dynamics/snapshots.csv` | `a` H₁, `b` ²D: normalised vs cell fraction |
| S5 Sampling | `supp_5/supp_5.jl` | `04_trial/samples.csv`, `snapshots.csv` | `a`–`c` ρ(sample, tumour) vs n for f = 1%, 0, 5%; `d` relative error of γ̂ |
| S6 Power law | `supp_6/supp_6.jl` | `01_dynamics/powerlaw.csv` | `a` γ̂ vs k_min with KS choice, `b` SE, `c` KS-selected k_min |

Conventions: time in cell cycles (≈ 4 days) since detection; schedule colours from ColorBrewer
Paired (MTD red, holidays blue, modulation green, metronomic purple, single round orange);
trajectories as in BirthDeathMutation_AT Figure 1 (total black, sensitive blue, resistant red,
treatment grey ∝ dose); time series show median and interquartile range over tumours while at
least half of them are still running; time to progression is restricted to the follow-up
(tumours that do not progress, including extinct ones, count as progression-free until `t_max`).
