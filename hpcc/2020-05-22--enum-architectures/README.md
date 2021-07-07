# Experiment - Enumerated genetic architectures

Configuration

- NUM_GENES=4
- GENE_SIZE=4
- GENOME_SIZE=16

- Bit flip rates: 0.003, 0.03, 0.3
- Locked genetic architectures (x65)

- 65 * 2 * 3

- Hand-designed architectures
  - Enumerate all possible architectures (for small gene count, gene size, and genome size)
    - Characterize distribution of graph properties
    - Evolve each for 10k generations in random static environment (x10)
    - Look at landscape (graph properties vs. fitness)

- Initial runs @ 30 replicates