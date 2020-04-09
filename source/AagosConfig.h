#ifndef AAGOS_CONFIG_H
#define AAGOS_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG(AagosConfig,
  GROUP(WORLD_STRUCTURE, "How should each organism's genome be setup?"),
    VALUE(CHANGE_RATE, size_t, 0, "How many changes to fitness tables each generation?"),
    VALUE(POP_SIZE, size_t, 1000, "How many organisms should be in the population?"),
    VALUE(MAX_GENS, size_t, 50000, "How many generations should the runs go for?"),
    VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
    VALUE(ELITE_COUNT, size_t, 0, "How many organisms should be selected via elite selection?"),
    VALUE(TOURNAMENT_SIZE, size_t, 2, "How many organisms should be chosen for each tournament?"),
    VALUE(GRADIENT_MODEL, bool, false, "Whether the current experiment uses a gradient model for fitness or trad. fitness"),

  GROUP(GENOME_STRUCTURE, "How should each organism's genome be setup?"),
    VALUE(NUM_BITS, size_t, 128, "Starting number of bits in each organism"),
    VALUE(NUM_GENES, size_t, 16, "Number of genes in each organism"),
    VALUE(GENE_SIZE, size_t, 8, "Size of each gene in each organism"),
    VALUE(MAX_SIZE, size_t, 1024, "maxiumum size of a genome"),
    VALUE(MIN_SIZE, size_t, 8, "minimum size of a genome"),

  GROUP(MUTATIONS, "Various mutation rates for Aagos"),
    VALUE(GENE_MOVE_PROB, double, 0.01, "Probability of each gene moving each generation"),
    VALUE(BIT_FLIP_PROB, double, 0.01, "Probability of each bit toggling"),
    VALUE(BIT_INS_PROB, double, 0.01, "Probability of a single bit being inserted."),
    VALUE(BIT_DEL_PROB, double, 0.01, "Probability of a single bit being removed."),

  GROUP(OUTPUT, "Output rates for Aagos"),
    VALUE(PRINT_INTERVAL, size_t, 1000, "How many updates between prints?"),
    VALUE(SUMMARY_INTERVAL, size_t, 1000, "How many updates between statistic gathering?"),
    VALUE(SNAPSHOT_INTERVAL, size_t, 10000, "How many updates between snapshots?"),
    VALUE(DATA_FILEPATH, std::string, "", "what directory should all data files be written to?")
)

#endif