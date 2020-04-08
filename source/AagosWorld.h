/// @todo Need to add cyclic environments!

#ifndef AAGOS_WORLD_H
#define AAGOS_WORLD_H

#include "Evolve/NK.h"
#include "Evolve/World.h"
#include "tools/Distribution.h"
#include "tools/math.h"
#include "tools/Range.h"
#include "tools/stats.h"
#include "tools/string_utils.h"

#include <sstream>
#include <string>

#include "AagosOrg.h"
#include "AagosConfig.h"

class AagosMutator {
protected:
  size_t num_genes;
  emp::Range<size_t> genome_size_constraints;
  double prob_gene_moves;
  double prob_bit_flip;
  double prob_bit_ins;
  double prob_bit_del;

  // Used for mutation?
  emp::Binomial gene_moves_binomial;
  emp::vector<emp::Binomial> bit_flips_binomials;
  emp::vector<emp::Binomial> inserts_binomials;
  emp::vector<emp::Binomial> deletes_binomials;

public:
  AagosMutator(size_t n_genes, const emp::Range<size_t> & genome_size,
               double p_gene_moves, double p_bit_flip, double p_bit_ins, double p_bit_del)
    : num_genes(n_genes),
      genome_size_constraints(genome_size),
      prob_gene_moves(p_gene_moves),
      prob_bit_flip(p_bit_flip),
      prob_bit_ins(p_bit_ins),
      prob_bit_del(p_bit_del),
      gene_moves_binomial(prob_gene_moves, num_genes)
  {
    for (size_t i = genome_size_constraints.GetLower(); i <= genome_size_constraints.GetUpper(); ++i) {
      bit_flips_binomials.emplace_back(prob_bit_flip, i);
      inserts_binomials.emplace_back(prob_bit_ins, i);
      deletes_binomials.emplace_back(prob_bit_del, i);
    }
  }

  /// Apply gene moves, single-bit substitutions, insertions, and deletions to org's genome.
  /// TODO - test that this works as expected
  size_t ApplyMutations(AagosOrg & org, emp::Random & random) {
    const size_t min_genome_size = genome_size_constraints.GetLower();
    const size_t max_genome_size = genome_size_constraints.GetUpper();
    emp_assert(org.GetNumBits() >= min_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");
    emp_assert(org.GetNumBits() <= max_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");

    auto & genome = org.GetGenome();
    size_t bin_array_offset = org.GetNumBits() - min_genome_size; // offset is num bits - min size of genome
    // Do gene moves
    const size_t num_moves = gene_moves_binomial.PickRandom(random);
    for (size_t m = 0; m < num_moves; ++m) {
      const size_t gene_id = random.GetUInt(0, num_genes);              // Pick a random gene
      genome.gene_starts[gene_id] = random.GetUInt(genome.bits.GetSize()); // Pick a random new location
    }
    // Do bit flips
    const size_t num_flips = bit_flips_binomials[bin_array_offset].PickRandom(random);
    for (size_t m = 0; m < num_flips; ++m) {
      const size_t pos = random.GetUInt(genome.bits.GetSize());
      genome.bits[pos] ^= 1;
    }
    // Do insertions and deletions.
    int num_insert = (int)inserts_binomials[bin_array_offset].PickRandom(random);
    int num_delete = (int)deletes_binomials[bin_array_offset].PickRandom(random);
    const int proj_size = (int)genome.bits.GetSize() + num_insert - num_delete;
    // Check gene size is within range.
    if (proj_size > (int)max_genome_size) {
      num_insert -= (proj_size - (int)max_genome_size);
    } else if (proj_size < (int)min_genome_size) {
      num_delete -= ((int)min_genome_size - proj_size);
    }
    // Assert that we'll be in size limitations.
    emp_assert((int)genome.bits.GetSize() + num_insert - num_delete >= (int)min_genome_size);
    emp_assert((int)genome.bits.GetSize() + num_insert - num_delete <= (int)max_genome_size);
    // Do insertions
    for (int i = 0; i < num_insert; ++i) {
      const size_t pos = random.GetUInt(org.GetNumBits()); // Figure out the position for insertion.
      genome.bits.Resize(genome.bits.GetSize() + 1);       // Increase genome size to make room for insertion.
      emp::BitVector mask(pos, 1);                         // Setup a mask to perserve early bits.
      mask.Resize(genome.bits.GetSize());                     // Align mask size.
      // Now build the new string!
      genome.bits = (mask & genome.bits) | ((genome.bits << 1) & ~mask);
      genome.bits[pos] = random.P(0.5); // Randomize the new bit.
      // Shift any genes that started at pos or later.
      for (auto & x : genome.gene_starts) {
        x += ((size_t)x >= pos);
      }
    }
    // Do deletions
    for (int i = 0; i < num_delete; ++i) {
      size_t pos = random.GetUInt(org.GetNumBits());
      emp::BitVector mask(pos, 1);
      mask.Resize(genome.bits.GetSize());
      genome.bits = (mask & genome.bits) | ((genome.bits >> 1) & ~mask);  // Build the new string!
      genome.bits.Resize(genome.bits.GetSize() - 1);                      // Decrease the size to account for deletion
      // Shift any genes that started at pos or later.
      if (pos == 0) {
        pos = 1; // Adjust position if beginning was deleted (don't want to subtract 1 if gene start = 0)
      }
      for (auto & x : genome.gene_starts) {
        x -= ((size_t)x >= pos);
      }
    }
    // Compute number of mutations, update organism's mutation-related tracking.
    const int num_muts = (int)num_moves + (int)num_flips + num_insert + num_delete;
    emp_assert(num_muts >= 0);
    if (num_muts > 0) {
      org.ResetHistogram();
    }
    return (size_t)num_muts;
  }

  // TODO - accessors
  // TODO - record mutations functionality

};


class AagosWorld : public emp::World<AagosOrg> {
public:
  // Convenience aliases
  using base_t = emp::World<AagosOrg>;
  using org_t = AagosOrg;
  using config_t = AagosConfig;
  using genome_t = AagosOrg::Genome;

  /// Fitness model for gradient fitness evaluation
  struct GradientFitnessModel {
    size_t num_genes;
    size_t gene_size;
    emp::vector<emp::BitVector> targets;

    GradientFitnessModel(emp::Random & rand, size_t n_genes, size_t g_size)
      : num_genes(n_genes), gene_size(g_size)
    {
      for (size_t i = 0; i < num_genes; ++i) {
        targets.emplace_back(emp::RandomBitVector(rand, gene_size));
        emp_assert(targets.back().GetSize() == gene_size);
      }
      emp_assert(targets.size() == num_genes);
    }

    const emp::BitVector & GetTarget(size_t id) const { return targets[id]; }
  };

  /// Fitness model for NK fitness evaluation
  struct NKFitnessModel {
    size_t num_genes;
    size_t gene_size;
    emp::NKLandscape landscape;

    NKFitnessModel(emp::Random & rand, size_t n_genes, size_t g_size)
      : num_genes(n_genes), gene_size(g_size)
    {
      landscape.Config(num_genes, gene_size-1, rand);
    }
  };

protected:
  const config_t & config;    ///< World configuration.

  emp::Ptr<NKFitnessModel> fitness_model_nk;
  emp::Ptr<GradientFitnessModel> fitness_model_gradient;
  std::function<void(org_t &)> evaluate_org;

  emp::Ptr<AagosMutator> mutator;

  // todo - data collection

  size_t gene_mask;
  size_t most_fit_id;

  void InitFitnessEval();
  void InitPop();

public:
  AagosWorld(emp::Random & random, const config_t & cfg)
    : base_t(random), config(cfg)
  {
    std::cout << "-- Constructing AagosWorld -- " << std::endl;

    // Basic setup
    gene_mask = emp::MaskLow<size_t>(config.GENE_SIZE());
    most_fit_id = 0;

    // Initialize fitness evaluation.
    std::cout << "Setting up fitness evaluation." << std::endl;
    InitFitnessEval();

    // Configure mutator
    std::cout << "Constructing mutator..." << std::endl;
    mutator = emp::NewPtr<AagosMutator>(config.NUM_GENES(), emp::Range<size_t>(config.MIN_SIZE(), config.MAX_SIZE()),
                                        config.GENE_MOVE_PROB(), config.BIT_FLIP_PROB(),
                                        config.BIT_INS_PROB(), config.BIT_DEL_PROB());
    std::cout << "  ...done constructing mutator." << std::endl;
    SetMutFun([this](org_t & org, emp::Random & rnd) {
      // NOTE - here's where we would intercept mutation-type distributions (with some extra infrastructure
      //        built into the mutator)!
      return mutator->ApplyMutations(org, rnd);
    });

    // TODO - initialize population
    std::cout << "Initialize the population" << std::endl;
    InitPop();

    SetPopStruct_Mixed(true);
    SetAutoMutate(config.ELITE_COUNT());          // Configure world to auto-mutate organisms (if id > elite count)

    // TODO - setup data tracking!
  }

  ~AagosWorld() {
    if (config.GRADIENT_MODEL()) fitness_model_gradient.Delete();
    mutator.Delete();
  }

  /// Advance world by a single time step (generation).
  void RunStep();

  /// Run world for configured number of generations.
  void Run();

};

void AagosWorld::RunStep() {
  // (1) evaluate population, (2) select parents, (3) update the world
  // == Do evaluation ==
  most_fit_id = 0;
  for (size_t org_id = 0; org_id < this->GetSize(); ++org_id) {
    emp_assert(IsOccupied(org_id));
    std::cout << "-- Evaluating org_id " << org_id << " --" << std::endl;
    evaluate_org(GetOrg(org_id));
    if (CalcFitnessID(org_id) > CalcFitnessID(most_fit_id)) {
      most_fit_id = org_id;
    }
    // NOTE - if we wanted to add phenotype tracking to systematics, here's where we could intercept
    //        the necessary information.
  }

  // == Do selection ==
  if (config.ELITE_COUNT()) emp::EliteSelect(*this, config.ELITE_COUNT(), 1);
  // Run a tournament for the rest...
  emp::TournamentSelect(*this, config.TOURNAMENT_SIZE(), config.POP_SIZE() - config.ELITE_COUNT());

  // == Do update ==
  // If it's a generation to print to console, do so
  const size_t update = GetUpdate();
  if (update % config.PRINT_INTERVAL() == 0) {
    std::cout << update
              << ": max fitness=" << CalcFitnessID(most_fit_id)
              << "; size=" << GetOrg(most_fit_id).GetNumBits()
              << std::endl;
    GetOrg(most_fit_id).Print();
    std::cout << std::endl;
  }

  Update();
  ClearCache();
}

void AagosWorld::Run() {
  for (size_t gen = 0; gen <= config.MAX_GENS(); ++gen) {
    RunStep();
  }
}

void AagosWorld::InitPop() {
  // Initialize population randomly (for now).
  for (size_t i = 0; i < config.POP_SIZE(); ++i) {
    genome_t genome(config.NUM_BITS(), config.NUM_GENES(), config.GENE_SIZE());
    genome.Randomize(*random_ptr);
    Inject(genome);
  }
  emp_assert(this->GetSize() == config.POP_SIZE());
}

void AagosWorld::InitFitnessEval() {
  // Fitness evaluation depends on configured fitness model.
  // Current model options: gradient, no gradient
  if (config.GRADIENT_MODEL()) {
    std::cout << "Initializing gradient model of fitness." << std::endl;
    fitness_model_gradient = emp::NewPtr<GradientFitnessModel>(*random_ptr, config.NUM_GENES(), config.GENE_SIZE());
    // Print out the gene targets
    std::cout << "Initial gene targets:" << std::endl;
    const auto & targets = fitness_model_gradient->targets;
    for (size_t i = 0; i < targets.size(); ++i) {
      std::cout << "  Target " << i << ": ";
      targets[i].Print();
      std::cout << std::endl;
    }

    evaluate_org = [this](org_t & org) {
      std::cout << "== Evaluating org ==" << std::endl;

      const size_t num_genes = config.NUM_GENES();
      const size_t num_bits = config.NUM_BITS();
      const size_t gene_size = config.GENE_SIZE();

      // Grab reference to and reset organism's phenotype.
      auto & phen = org.GetPhenotype();
      phen.Reset();
      // Calculate fitness contribution of each gene independently.
      double fitness = 0.0;
      const auto & gene_starts = org.GetGeneStarts();
      for (size_t gene_id = 0; gene_id < num_genes; ++gene_id) {
        emp_assert(gene_id < gene_starts.size());
        const size_t gene_start = gene_starts[gene_id];
        uint32_t gene_val = org.GetBits().GetUIntAtBit(gene_start) & gene_mask;
        const size_t tail_bits = num_bits - gene_start;
        // If a gene runs off the end of the bitstring, loop around to the beginning.
        if (tail_bits < gene_size) {
          gene_val |= (org.GetBits().GetUInt(0) << tail_bits) & gene_mask;
        }
        // Compute fitness contribution of this gene
        // - Remember, we assume the first index of gene_starts maps to the first index of the target bitstring.
        emp_assert(gene_starts.size() == fitness_model_gradient->targets.size());
        const emp::BitVector & target = fitness_model_gradient->GetTarget(gene_id);
        emp::BitVector gene = emp::BitVector(gene_size);
        gene.SetUInt(0, gene_val);
        const double fitness_contribution = (double)target.EQU(gene).count() / (double)gene_size;
        phen.gene_fitness_contributions[gene_id] = fitness_contribution;
        fitness += fitness_contribution;
      }
      phen.fitness = fitness;
    };
  } else {
    // todo
    std::cout << "TODO" << std::endl;
  }
  // Note that this assumes that this organism has been evaluated.
  SetFitFun([this](org_t & org) {
    return org.GetPhenotype().fitness;
  });
}


#endif
