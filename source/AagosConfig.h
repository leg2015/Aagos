#ifndef AAGOS_CONFIG_H
#define AAGOS_CONFIG_H

#include "emp/config/config.hpp"

EMP_BUILD_CONFIG(AagosConfig,
  GROUP(WORLD_STRUCTURE, "How should each organism's genome be setup?"),
    VALUE(CHANGE_MAGNITUDE, size_t, 0, "How many changes to fitness tables each generation?"),
    VALUE(CHANGE_FREQUENCY, size_t, 1, "How many generations elapse between environment changes? Frequency = 0 means no changes."),
    VALUE(POP_SIZE, size_t, 1000, "How many organisms should be in the population?"),
    VALUE(MAX_GENS, size_t, 50000, "How many generations should the runs go for?"),
    VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
    VALUE(TOURNAMENT_SIZE, size_t, 2, "How many organisms should be chosen for each tournament?"),
    VALUE(GRADIENT_MODEL, bool, false, "Whether the current experiment uses a gradient model for fitness or trad. fitness"),
    VALUE(LOAD_ANCESTOR, bool, false, "Should we initialize population with ancestor genotype from file?"),
    VALUE(LOAD_ANCESTOR_FILE, std::string, "ancestor.csv", "File to load ancestor genotype from"),
    VALUE(RANDOMIZE_LOAD_ANCESTOR_BITS, bool, false, "Should we randomize the bit values for loaded ancestor?"),
    VALUE(LOAD_ENV_FROM_FILE, bool, false, "Should we load the environment from a file?"),
    VALUE(LOAD_ENV_FILE, std::string, "environment.env", "File to load environment from (if configured to load)"),

  GROUP(RUN_SECOND_PHASE, "Will run have a second phase with new configuration parameters? (limited set of things can change)"),
    VALUE(PHASE_2_ACTIVE, bool, false, "Should run continue to a second phase with new parameters?"),
    VALUE(PHASE_2_LOAD_ENV_FROM_FILE, bool, false, "Should we load initial phase 2 environment from a file?"),
    VALUE(PHASE_2_ENV_FILE, std::string, "environment.env", "File to load environment from (if configured to load phase 2 environment from file)"),
    VALUE(PHASE_2_CHANGE_MAGNITUDE, size_t, 0, "Change magnitude for second phase of evolution"),
    VALUE(PHASE_2_CHANGE_FREQUENCY, size_t, 1, "Change frequency for the second phase of evolution"),
    VALUE(PHASE_2_MAX_GENS, size_t, 100, "Number of generations in second phase of evolution"),
    VALUE(PHASE_2_TOURNAMENT_SIZE, size_t, 2, "How many organisms should be chosen for each tournament during phase 2 of evolution?"),
    VALUE(PHASE_2_GENE_MOVE_PROB, double, 0.01, "GENE_MOVE_PROB for second phase of evolution"),
    VALUE(PHASE_2_BIT_FLIP_PROB, double, 0.01, "BIT_FLIP_PROB for second phase of evolution"),
    VALUE(PHASE_2_BIT_INS_PROB, double, 0.01, "BIT_INS_PROB for second phase of evolution"),
    VALUE(PHASE_2_BIT_DEL_PROB, double, 0.01, "BIT_DEL_PROB for second phase of evolution"),

  GROUP(GENOME_STRUCTURE, "How should each organism's genome be setup?"),
    VALUE(NUM_BITS, size_t, 128, "Starting number of bits in each organism"),
    VALUE(NUM_GENES, size_t, 16, "Number of genes in each organism"),
    VALUE(GENE_SIZE, size_t, 8, "Size of each gene in each organism"),
    VALUE(MAX_SIZE, size_t, 1024, "maxiumum size of a genome"),
    VALUE(MIN_SIZE, size_t, 8, "minimum size of a genome"),

  GROUP(MUTATIONS, "Various mutation rates for Aagos"),
    VALUE(APPLY_BIT_MUTS_PER_GENE, bool, false, "Control. apply bit mutations to coding sites at a per-gene-per-site rate instead of just per-site"),
    VALUE(GENE_MOVE_PROB, double, 0.01, "Probability of each gene moving each generation"),
    VALUE(BIT_FLIP_PROB, double, 0.01, "Probability of each bit toggling"),
    VALUE(BIT_INS_PROB, double, 0.01, "Probability of a single bit being inserted."),
    VALUE(BIT_DEL_PROB, double, 0.01, "Probability of a single bit being removed."),

  GROUP(OUTPUT, "Output rates for Aagos"),
    VALUE(PRINT_INTERVAL, size_t, 1000, "How many updates between prints?"),
    VALUE(SUMMARY_INTERVAL, size_t, 1000, "How many updates between statistic gathering?"),
    VALUE(SNAPSHOT_INTERVAL, size_t, 10000, "How many updates between snapshots?"),
    VALUE(PHYLOGENY_TRACKING, bool, true, "Should we collect phylogeny data?"),
    VALUE(DATA_FILEPATH, std::string, "./output/", "what directory should all data files be written to?")
)

#endif