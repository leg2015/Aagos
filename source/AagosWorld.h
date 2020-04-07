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
    std::cout << "Generating binomials for mutation...";
    for (size_t i = genome_size_constraints.GetLower(); i <= genome_size_constraints.GetUpper(); ++i) {
      bit_flips_binomials.emplace_back(prob_bit_flip, i);
      inserts_binomials.emplace_back(prob_bit_ins, i);
      deletes_binomials.emplace_back(prob_bit_del, i);
    }
    std::cout << " ...done" << std::endl;
  }

  // TODO - apply mutations function
  // TODO - accessors
  // TODO - record mutations functionality

};


class AagosWorld : public emp::World<AagosOrg> {
public:
  // Convenience aliases
  using base_t = emp::World<AagosOrg>;
  using org_t = AagosOrg;
  using config_t = AagosConfig;

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
  std::function<double(const org_t &)> calc_fitness;

  emp::Ptr<AagosMutator> mutator;

  // todo - data collection

  size_t gene_mask;
  size_t most_fit_id;

  void InitFitnessEval();

public:
  AagosWorld(emp::Random & random, const config_t & cfg)
    : base_t(random), config(cfg)
  {
    std::cout << "Constructing AagosWorld" << std::endl;
    // Basic setup.
    // Why gene size - 1? Because K is the number of other sites this site depends on.

    // todo - configure mutator
    mutator = emp::NewPtr<AagosMutator>(config.NUM_GENES(), emp::Range<size_t>(config.MIN_SIZE(), config.MAX_SIZE()),
                                        config.GENE_MOVE_PROB(), config.BIT_FLIP_PROB(),
                                        config.BIT_INS_PROB(), config.BIT_DEL_PROB());
    // Initialize fitness evaluation.
    InitFitnessEval();
  }

  ~AagosWorld() {
    if (config.GRADIENT_MODEL()) fitness_model_gradient.Delete();
    mutator.Delete();
  }

};

void AagosWorld::InitFitnessEval() {
  // Fitness evaluation depends on configured fitness model.
  // Current model options: gradient, no gradient
  if (config.GRADIENT_MODEL()) {
    std::cout << "Initializing gradient model of fitness." << std::endl;
    fitness_model_gradient = emp::NewPtr<GradientFitnessModel>(GetRandom(), config.NUM_GENES(), config.GENE_SIZE());
    calc_fitness = [this](const org_t & org) {
      const size_t num_genes = config.NUM_GENES();
      double fitness = 0.0;
      // Calculate fitness contribution of each gene independently.
      for (size_t gene_id = 0; gene_id < num_genes; ++gene_id) {
        // const size_t gene_pos = org.gene_starts
      }
      return fitness;
    };
  } else {
    // todo
    std::cout << "TODO" << std::endl;
  }
}


#endif
