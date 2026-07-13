# SymmetricRings benchmark report

## Overall summary

| Run health | Families | Cases | Improvements | Regressions | Stable | New |
|---|---:|---:|---:|---:|---:|---:|
| **clean** | 4 | 24 | 0 | 0 | 0 | 24 |

## Baseline comparison

### BasisConversion-extra

| Cases | Improvements | Regressions | Stable | New |
|---:|---:|---:|---:|---:|
| 5 | 0 | 0 | 0 | 5 |

| case_id | operation | tier | coefficient_ring | runs | current_median | best_valid_baseline | change_vs_best_percent | most_recent_baseline | change_vs_recent_percent | classification |
|---|---|---|---|---|---|---|---|---|---|---|
| extra-conv-p-combination-S | PowerSumCombinationToSchur | Medium | FracQQt | 3 | 0.001268 | NA | NA | NA | NA | new |
| extra-conv-S-to-h | DirectBasisConversion | Small | FracQQt | 3 | 0.000353 | NA | NA | NA | NA | new |
| extra-conv-p-parameter-S | ParameterPowerSumCombinationToSchur | Medium | FracQQt | 3 | 0.001402 | NA | NA | NA | NA | new |
| extra-conv-h-to-S | DirectBasisConversion | Small | FracQQt | 3 | 0.00054 | NA | NA | NA | NA | new |
| extra-conv-S-combination-p | SchurCombinationToPowerSums | Medium | FracQQt | 3 | 0.006775 | NA | NA | NA | NA | new |

### HallLittlewoodProduct-extra

| Cases | Improvements | Regressions | Stable | New |
|---:|---:|---:|---:|---:|
| 6 | 0 | 0 | 0 | 6 |

| case_id | operation | tier | coefficient_ring | runs | current_median | best_valid_baseline | change_vs_best_percent | most_recent_baseline | change_vs_recent_percent | classification |
|---|---|---|---|---|---|---|---|---|---|---|
| extra-hl-Pomega-retained | HallLittlewoodBasisProductRetained | Medium | FracQQt | 3 | 0.093682 | NA | NA | NA | NA | new |
| extra-hl-P-retained | HallLittlewoodBasisProductRetained | Medium | FracQQt | 3 | 0.093075 | NA | NA | NA | NA | new |
| extra-hl-B-retained | HallLittlewoodBasisProductRetained | Medium | FracQQt | 3 | 0.038547 | NA | NA | NA | NA | new |
| extra-hl-Pomega-ordinary | HallLittlewoodBasisProduct | Medium | FracQQt | 3 | 0.408745 | NA | NA | NA | NA | new |
| extra-hl-P-ordinary | HallLittlewoodBasisProduct | Medium | FracQQt | 3 | 0.404738 | NA | NA | NA | NA | new |
| extra-hl-B-ordinary | HallLittlewoodBasisProduct | Medium | FracQQt | 3 | 0.249952 | NA | NA | NA | NA | new |

### InnerProduct-extra

| Cases | Improvements | Regressions | Stable | New |
|---:|---:|---:|---:|---:|
| 10 | 0 | 0 | 0 | 10 |

| case_id | operation | tier | coefficient_ring | runs | current_median | best_valid_baseline | change_vs_best_percent | most_recent_baseline | change_vs_recent_percent | classification |
|---|---|---|---|---|---|---|---|---|---|---|
| extra-inner-S-e | SchurClassicalInnerProduct | Medium | QQ | 3 | 0.000227 | NA | NA | NA | NA | new |
| extra-inner-Q-P-diagonal | HallLittlewoodDiagonalInnerProduct | Medium | FracQQt | 3 | 0.000283 | NA | NA | NA | NA | new |
| extra-inner-b-ff-diagonal | GeneratorDualInnerProduct | Small | FracQQt | 3 | 0.000319 | NA | NA | NA | NA | new |
| extra-inner-Somega-h | SchurClassicalInnerProduct | Medium | QQ | 3 | 0.000315 | NA | NA | NA | NA | new |
| extra-inner-ordinary-expansions | OrdinaryExpansionInnerProduct | Medium | QQ | 3 | 0.001976 | NA | NA | NA | NA | new |
| extra-inner-S-h | SchurClassicalInnerProduct | Medium | QQ | 3 | 0.000233 | NA | NA | NA | NA | new |
| extra-inner-q-m-diagonal | GeneratorDualInnerProduct | Small | FracQQt | 3 | 0.000311 | NA | NA | NA | NA | new |
| extra-inner-B-Pomega-diagonal | HallLittlewoodDiagonalInnerProduct | Medium | FracQQt | 3 | 0.000276 | NA | NA | NA | NA | new |
| extra-inner-Somega-e | SchurClassicalInnerProduct | Medium | QQ | 3 | 0.000321 | NA | NA | NA | NA | new |
| extra-inner-targeted-Q-P | TargetedHallLittlewoodInnerProduct | Medium | FracQQt | 3 | 0.039034 | NA | NA | NA | NA | new |

### SchurProduct-extra

| Cases | Improvements | Regressions | Stable | New |
|---:|---:|---:|---:|---:|
| 3 | 0 | 0 | 0 | 3 |

| case_id | operation | tier | coefficient_ring | runs | current_median | best_valid_baseline | change_vs_best_percent | most_recent_baseline | change_vs_recent_percent | classification |
|---|---|---|---|---|---|---|---|---|---|---|
| extra-prod-multiply-to-S | SchurProductMultiplyToBasis | Medium | FracQQt | 3 | 0.00042 | NA | NA | NA | NA | new |
| extra-prod-post-plethysm | PlethysmSchurProductToSchur | Large | FracQQt | 3 | 0.034576 | NA | NA | NA | NA | new |
| extra-prod-direct-lr | SchurProductExpanded | Medium | FracQQt | 3 | 0.000494 | NA | NA | NA | NA | new |

## System

| Field | Value |
|---|---|
| captured_at | 2026-07-12 20:25:09 EDT |
| computer_model | MacBook Air |
| model_identifier | Mac16,13 |
| cpu_model | Apple M4 |
| physical_cores | 10 (4 performance and 6 efficiency) |
| logical_cores | 10 (4 performance and 6 efficiency) |
| memory | 16 GB |
| memory_bytes | unknown |
| os_name | macOS |
| os_version | 15.7.7 |
| os_build | 24G720 |
| kernel | Darwin 24.6.0 |
| architecture | arm64 |
| m2_version | 1.26.06-26-ga4df3f098a-dirty (stable) |
| m2_executable | /Users/johngraf/M2Dev/Project-SymFcns/M2/M2/BUILD/build/M2 |
| git_commit | afe03aab1b7eaf1aa3b00303d738a9e119abd947 |
| git_branch | stable |
| git_dirty | true |
| cmake_build_type | Release |

## Run conditions

Overall health: **clean**.

- **Clean:** conditions were stable and no acceptance concern was detected.
- **Warning:** moderate drift or resource pressure was detected; review before accepting baselines.
- **Compromised:** severe drift, thermal pressure, or a performance-limiting mode makes the run unsuitable as a baseline.

No run-condition warnings were detected.

Pageout activity: 81 pages (1.27 MiB).

### Benchmark selection

| Mode | Level | Random seed | New only | Reference run |
|---|---|---:|---|---|
| none | none | none | true | 20260712-190928 |

### Condition snapshots

| Metric | Before | After |
|---|---|---|
| captured_at | 2026-07-12 20:25:09 EDT | 2026-07-12 20:25:40 EDT |
| snapshot_epoch | unknown | unknown |
| power_source | AC Power | AC Power |
| battery_state | 40%; discharging; 3:57 remaining present: true | 40%; discharging; 3:47 remaining present: true |
| low_power_mode | disabled | disabled |
| thermal_state | nominal | nominal |
| memory_free_percent | 75 | 73 |
| pageouts | 1131541 | 1131622 |
| page_size_bytes | unknown | unknown |
| load_averages | 2.10 2.09 2.01 | 2.14 2.10 2.02 |

### Family calibration probes

The fixed `conv-S-three-p` probe runs in a fresh process immediately before and after each family; no monitor runs concurrently with timed cases.

| Family | Before CPU (s) | After CPU (s) | Drift | Before wall (s) | After wall (s) |
|---|---:|---:|---:|---:|---:|
| BasisConversion-extra | .18651 | .193873 | +3.95% | .186498 | .193872 |
| SchurProduct-extra | .192972 | .190993 | -1.03% | .192992 | .191004 |
| HallLittlewoodProduct-extra | .190463 | .195666 | +2.73% | .190492 | .195695 |
| InnerProduct-extra | .195008 | .190827 | -2.14% | .194999 | .190828 |

