# Biological Ranking Pilot

This is a real current-pipeline pilot, not a trained-model result. Four
TM-align + Rosetta structures generated for `1FGNH/1TFHA` were compared with
native `1ahw.pdb`, mapping model chains `A,H` to native chains `B,C`.

| Candidate | Baseline score | iRMSD | DockQ |
| --- | ---: | ---: | --- |
| `1h5bAB_o1` | 0.611420 | 16.937 | 0.031 |
| `3lqmAB_o1` | 0.553352 | 27.240 | 0.004 |
| `3lqmAB_o2` | 0.561563 | 31.949 | 0.006 |
| `2f0xEH_o2` | 0.458931 | 33.521 | 0.004 |

The deterministic baseline selected `1h5bAB_o1`, which also had the lowest
iRMSD and highest DockQ in this case. However, every candidate is below the
native-like threshold (`DockQ >= 0.23`), so the ranking improved relative
quality without producing a successful native-like prediction. This is an
observation for one native complex, not proof of generalization. The baseline
scores and iRMSD values are preserved in:

- `tmp/agent/20260711-biological-ranking/current-1ahw-ranked.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-irmsd.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-dockq-original.csv`
- `tmp/agent/20260711-biological-ranking/current-1ahw-labeled.csv`

The next valid experiment is to repeat this workflow over independent native
complexes, split by complex/sequence cluster, and compare the deterministic
baseline with the tabular and contact-graph models.
