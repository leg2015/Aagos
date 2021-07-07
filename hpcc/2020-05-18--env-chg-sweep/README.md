# 2020-05-18 - Environment Change Rate Sweep Experiments

Goals outlined in [notes/2020-05-18--redesign.md](../../notes/2020-05-18--redesign.md) -
this is experiment (1) on the list of experiments.

## Outcome

## Plan

- (1) Test hypothesis: a changing environment promotes gene spread
  - Measure gene overlap/number of coding sites at a range of environment change rates
  - Experiment (15 x 2 = 30 conditions):
    - Change rate (approx) (x17): 0, 1/256, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2/1, 4/1, 8/1, 16/1, drift
    - Fitness model (x2): NK, Gradient
    - Selection: tournament size 8, 1 (for drift control => drift requires no change because environment is meaningless)
    - Bit flip rate: 0.003
    - Comparisons
      - (1) Number of coding sites across change rates (+ drift)
        - Expectation: higher change rates (to a point) have less overlap
      - (2) Fitnesses after environment lock down across change rates (+ drift)
        - Expectation: more modular genomes have higher fitness
      - (3) Fitness (after lock-down) vs. modularity comparison

## Directory guide

Backbone configuration is given in Aagos.cfg. Note that configuration in Aagos.cfg is overridden by
commandline arguments. Treatments will differ by commandline arguments.

gen-chg-rate-sub-script.py generates the HPCC submission script for submitting everything.
