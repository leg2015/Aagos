# 2020-04-16 Experiments

NOTES

- These runs happened prior to systematics tracking and prior to capacity to change tournament size
  during phase two of evolution.

Goal: reproduce qualitative results from our 2018/2019 experiments with refactored code.

- Re-do full sweep of change rates and mutation rates for NK fitness model
- Perform full sweep of change rates and mutation rates for gradient fitness model

Because each individual run is fairly short (minutes) and we have so many conditions to run (@ 50 reps each!),
the submission process is as follows: one slurm script that submits job array, 1 job per condition;
each job runs a python script that maps the array ID to a particular run configuration; the python
script sequentially runs 50 replicates (each with unique random number seed) of each condition.

- [submit.sb](submit.sb) - HPCC slurm submission script (point of entry)
- [run.py](run.py) - python script to manage running multiple replicates of each condition
- [Aagos.cfg](Aagos.cfg) - base configuration file for these runs

Random seeds used = [1001:10800]

## Outcome

Gradient model produced qualitatively similar results to original NK model. See associated analyses.

## Plan

### NK Fitness Model Conditions

- Change magnitudes (14): 0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096
- Change frequency: 1
- Bit flip mutation rates (7): 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1
- ~~Tournament size (i.e., pure drift?): 1 (pure drift), 8~~ (already too many combinations, need to scale back)
- Replicates: 30 or 50?
- Phase 1 Generations = 50,000
  - -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001
- Phase 2 Generations = 10,000, locked down architecture
  - -GENE_MOVE_PROB 0 -BIT_INS_PROB 0 -BIT_DEL_PROB 0

### Gradient Fitness Model Conditions

- Change:
  - (nk chg mag = 0) magnitude = 0; frequency = 0
  - (nk chg mag = 1) magnitude = 1; frequency = 256
  - (nk chg mag = 2) magnitude = 1; frequency = 128
  - (nk chg mag = 4) magnitude = 1; frequency = 64
  - (nk chg mag = 8) magnitude = 1; frequency = 32
  - (nk chg mag = 16) magnitude = 1; frequency = 16
  - (nk chg mag = 32) magnitude = 1; frequency = 8
  - (nk chg mag = 64) magnitude = 1; frequency = 4
  - (nk chg mag = 128) magnitude = 1; frequency = 2
  - (nk chg mag = 256) magnitude = 1; frequency = 1
  - (nk chg mag = 512) magnitude = 2; frequency = 1
  - (nk chg mag = 1024) magnitude = 4; frequency = 1
  - (nk chg mag = 2048) magnitude = 8; frequency = 1
  - (nk chg mag = 4096) magnitude = 16; frequency = 1
- Bit flip mutation rates (7): 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1
- Replicates: 30 or 50?
- Phase 1 Generations = 50,000
  - -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001
- Phase 2 Generations = 10,000, locked down architecture, randomized static environment
  - -GENE_MOVE_PROB 0 -BIT_INS_PROB 0 -BIT_DEL_PROB 0