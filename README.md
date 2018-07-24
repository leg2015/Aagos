# Aagos

## Motivation

* In static systems, high mutation rates create a pressure for gene overlap and smaller genomes
* However, more modular, separated genetic architectures have been shown to be more evolvable in previous work
![](https://github.com/leg2015/Aagos/blob/master/Background.png "Background")
* Furthermore, in evolutionary biology, we see a general trend of larger essential genomes in more complex organisms
* Therefore, _what environmental pressures have selected for modular, evolvable genetic architectures in biological systems?_

## System Model
![](https://github.com/leg2015/Aagos/blob/master/NKLandscape.png "Nk-Inspried Landscape")
* NK-Inspired landscape
* Variable length genomes


## Results

* High mutation rates in a static landscape do indeed lead to smaller minimal genomes
* Moderate environmental change leads to larger minimal genomes, however high environmental change leads to meltdown
* Even with a high mutation rate, in the presence of moderate environmental change, minimal genome size is significantly larger than in static environments

## Result Replication

To run these experiments for yourself, simply download the [Empirical library](https://github.com/devosoft/Empirical) and [this repository](https://github.com/leg2015/Aagos) and enter the commands:

`make`

`./Aagos -[parameters]`

**Environmental Parameters**
  * CHANGE_RATE, default 0, How many changes to fitness tables each generation
  * POP_SIZE, default 1000, How many organisms should be in the population
  * MAX_GENS, default 50000, How many generations should the runs go for
  * SEED, default 0, Random number seed (0 for based on time)
  * ELITE_COUNT, default 0, How many organisms should be selected via elite selection
  * TOURNAMENT_SIZE, default 2, How many organisms should be chosen for each tournament

**Genomic Structure**

  * NUM_BITS, default 128, Starting number of bits in each organism
  * NUM_GENES, default 16, Number of genes in each organism
  * GENE_SIZE, default 8, Size of each gene in each organism
  * MAX_SIZE, default 1024, maxiumum size of a genome
  * MIN_SIZE, default 8, minimum size of a genome

**Mutations**

  * GENE_MOVE_PROB, default 0.01, Probability of each gene moving each generation
  * BIT_FLIP_PROB, default 0.01, Probability of each bit toggling
  * BIT_INS_PROB, default 0.01, Probability of a single bit being inserted
  * BIT_DEL_PROB, default 0.01, Probability of a single bit being removed
  
**Output**

  * PRINT_INTERVAL, default 1000, How many updates between prints?
                 
