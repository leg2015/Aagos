# 2020-04-23 Experiments

- Pivot mutation rate: 0.003
- Pivot change rates: 64, 2, 0
- Fitness models: gradient model & nk model
- At anchor points, gene counts: 8, 16, 32 (shift starting genome size appropriately)
- Random drift treatment for each mutation rate, gene count, & starting genome size
  - + turn selection on after 50k generations of drift! TODO
- Large replicate count: 100

## Outcome

Development of the web visualization revealed a bug in how we computed the fitness contribution of genes
that overlap the edge of a genome. Specifically, when genomes are any size other than the starting genome
size, any part of a gene that wraps around to the beginning of the gene was mis-read, resulting in
an inaccurate fitness contribution calculation.

As a result, we need to re-run our experiments, and we are using this as an opportunity to reduce the
number of total runs we need to do. The plans described below was taking upward of a month to finish,
which is significantly longer than previous experiments because of the addition of phylogeny tracking.

## Plans

- Runs:
  - **Vary mutation rate, static environment (14 conditions)**
    - BIT_FLIP_PROB = [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1]
    - Environment = "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0"
    - Genome configuration = "-NUM_BITS 128 -NUM_GENES 16 -MAX_SIZE 1024"
    - Fitness models = [gradient, nk]
  - **Varying mutation rates, change rate of 2 (14 conditions)**
    - BIT_FLIP_PROB = [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1]
    - Environment
      - gradient = "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 128",
      - nk = "-CHANGE_MAGNITUDE 2 -CHANGE_FREQUENCY 1",
    - Genome configuration = "-NUM_BITS 128 -NUM_GENES 16 -MAX_SIZE 1024"
    - Fitness models = [gradient, nk]
  - **Varying bit flip mutation, change rate of 64 (14 conditions)**
    - BIT_FLIP_PROB (x7) = [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1]
    - Environment
      - gradient = "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 4",
      - nk = "-CHANGE_MAGNITUDE 64 -CHANGE_FREQUENCY 1",
    - Genome configuration = "-NUM_BITS 128 -NUM_GENES 16 -MAX_SIZE 1024"
    - Fitness models (x2) = [gradient, nk]
  - **Mutation rate = 0.003, varying change rate (2x14x3=84 conditions)**
    - BIT_FLIP_PROB = [0.003]
    - Fitness models (x2) = [gradient, nk]
    - Environment
      - gradient (x14) = ["-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 256", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 128", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 64", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 32", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 16", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 8", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 4", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 2", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 2 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 4 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 8 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 16 -CHANGE_FREQUENCY 1"]
      - nk (x14) = ["-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0", "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 2 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 4 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 8 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 16 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 32 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 64 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 128 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 256 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 512 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 1024 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 2048 -CHANGE_FREQUENCY 1", "-CHANGE_MAGNITUDE 4096 -CHANGE_FREQUENCY 1"]
    - Genome configuration (x3)
      - Num genes = [8, 16, 32]
  - **Genetic drift controls (42 conditions)**
    - BIT_FLIP_PROB = [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1]
    - Environment: "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0"
      - note that this doesn't matter for random drift runs
    - SELECTION: "-TOURNAMENT_SIZE 1 -PHASE_2_TOURNAMENT_SIZE 8"
    - Fitness model: ["-GRADIENT_MODEL 0", "-GRADIENT_MODEL 1"]
    - Genome configuration (x3)
      - Num genes = [8, 16, 32]
