# 2020-05-20 - Experiments - Mutation rate sweep

- (2) Test hypothesis: high mutation rates (relative to lower mutation rates) promote gene
  overlap
  - Measure gene overlap/number of coding sites at a range of mutation rates
  - Experiment (6 x 5 x 2 = 60 conditions):
    - Bit flip rate (x6): 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1
    - Change rate (approx) (x5): 0, 1/128, 1/4, 4, drift
    - Fitness model (x2): NK, gradient
    - Selection: tournament size 8, 1 (for drift)