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
#include "control/Signal.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>

#include "AagosOrg.h"
#include "AagosConfig.h"

class AagosMutator {
public:
  enum class MUTATION_TYPES {
    BIT_FLIPS=0,
    BIT_INSERTIONS,
    BIT_DELETIONS,
    GENE_MOVES,
  };

protected:
  const size_t num_genes;
  const emp::Range<size_t> genome_size_constraints;
  const double prob_gene_moves;
  const double prob_bit_flip;
  const double prob_bit_ins;
  const double prob_bit_del;

  // Used for mutation?
  emp::Binomial gene_moves_binomial;
  emp::vector<emp::Binomial> bit_flips_binomials;
  emp::vector<emp::Binomial> inserts_binomials;
  emp::vector<emp::Binomial> deletes_binomials;

  // Mutation tracking
  std::unordered_map<MUTATION_TYPES, int> last_mutation_tracker;

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
      const size_t gene_id = random.GetUInt(0, num_genes);                 // Pick a random gene
      genome.gene_starts[gene_id] = random.GetUInt(genome.bits.GetSize()); // Pick a random new location
    }
    last_mutation_tracker[MUTATION_TYPES::GENE_MOVES] = (int)num_moves;

    // Do bit flips
    const size_t num_flips = bit_flips_binomials[bin_array_offset].PickRandom(random);
    for (size_t m = 0; m < num_flips; ++m) {
      const size_t pos = random.GetUInt(genome.bits.GetSize());
      genome.bits[pos] ^= 1;
    }
    last_mutation_tracker[MUTATION_TYPES::BIT_FLIPS] = (int)num_flips;

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
    last_mutation_tracker[MUTATION_TYPES::BIT_INSERTIONS] = num_insert;
    last_mutation_tracker[MUTATION_TYPES::BIT_DELETIONS] = num_delete;

    // Compute number of mutations, update organism's mutation-related tracking.
    const int num_muts = (int)num_moves + (int)num_flips + num_insert + num_delete;
    emp_assert(num_muts >= 0);
    if (num_muts > 0) {
      org.ResetHistogram();
    }
    return (size_t)num_muts;
  }

  // TODO - accessors
  std::unordered_map<MUTATION_TYPES, int> & GetLastMutations() { return last_mutation_tracker; }

  void ResetLastMutationTracker() {
    last_mutation_tracker[MUTATION_TYPES::BIT_FLIPS] = 0;
    last_mutation_tracker[MUTATION_TYPES::BIT_INSERTIONS] = 0;
    last_mutation_tracker[MUTATION_TYPES::BIT_DELETIONS] = 0;
    last_mutation_tracker[MUTATION_TYPES::GENE_MOVES] = 0;
  }

};


class AagosWorld : public emp::World<AagosOrg> {
public:
  // Convenience aliases
  using base_t = emp::World<AagosOrg>;
  using org_t = AagosOrg;
  using config_t = AagosConfig;
  using genome_t = AagosOrg::Genome;
  using phenotype_t = AagosOrg::Phenotype;
  using mutator_t = AagosMutator;

  using mut_landscape_t = emp::datastruct::mut_landscape_info<phenotype_t>;
  using systematics_t = emp::Systematics<org_t, genome_t, mut_landscape_t>;
  using taxon_t = typename systematics_t::taxon_t;

  // todo - make different fitness models derive from a base class that defines functions like randomize,
  //        change, load_from_file, etc
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

    /// Mutate a number of target bits equal to bit cnt.
    void RandomizeTargetBits(emp::Random & rand, size_t bit_cnt) {
      for (size_t i = 0; i < bit_cnt; ++i) {
        // Select a random target sequence.
        const size_t target_id = rand.GetUInt(targets.size());
        emp::BitVector & target = targets[target_id];
        const size_t target_pos = rand.GetUInt(target.GetSize());
        target.Set(target_pos, !target.Get(target_pos));
      }
    }

    /// Randomize a number of targets equal to cnt.
    void RandomizeTargets(emp::Random & rand, size_t cnt) {
      // Change a number of targets (= cnt).
      emp::vector<size_t> target_ids(targets.size());
      std::iota(target_ids.begin(), target_ids.end(), 0);
      emp::Shuffle(rand, target_ids);
      cnt = emp::Min(cnt, target_ids.size()); // No need to randomize more targets than exist.
      for (size_t i = 0; i < cnt; ++i) {
        // Select a random target sequence.
        emp_assert(i < target_ids.size());
        const size_t target_id = target_ids[i];
        emp::BitVector & target = targets[target_id];
        emp::RandomizeBitVector(target, rand);
      }
    }

    void PrintTargets(std::ostream & out=std::cout) {
      // lazily outsource to emp::to_string
      out << emp::to_string(targets);
    }

    /// Load targets from file (specified by given path)
    /// First non-commented line of file should give what emp::to_string would output.
    ///  - [ bits bits bits bits ]
    ///  - Loaded environment must be consistent with num_genes and gene_size
    bool LoadTargets(const std::string & path) {
      std::ifstream environment_fstream(path);
      if (!environment_fstream.is_open()) {
        std::cout << "Failed to open environment file (" << path << "). Exiting..." << std::endl;
        exit(-1);
      }
      std::string cur_line;
      emp::vector<std::string> line_components;
      bool success = false;
      while (!environment_fstream.eof()) {
        std::getline(environment_fstream, cur_line);
        emp::left_justify(cur_line); // Remove leading whitespace
        if (cur_line == emp::empty_string()) continue; // Skip empty line.
        else if (cur_line[0] == '#') continue; // Treat '#' as a commented line
        else if (cur_line[0] == '[') {
          // Attempt to read environments state.
          emp::remove_chars(cur_line, "[]"); // Remove brackets
          emp::left_justify(cur_line);       // Remove leading whitespace
          emp::right_justify(cur_line);      // Remove trailing whitespace
          line_components.clear();
          emp::slice(cur_line, line_components, ' ');
          // Each slice should be a gene target, check to make sure number of targets is correct.
          if (line_components.size() != num_genes) return false;

          for (size_t i = 0; i < line_components.size(); ++i) {
            emp_assert(i < targets.size());
            auto & target_str = line_components[i];
            // Target string should be correct size.
            if (target_str.size() != gene_size) return false;
            for (size_t bit = 0; bit < target_str.size(); ++bit) {
              emp_assert(bit < targets[i].GetSize());
              // if (component == '1')
              targets[i].Set(targets[i].GetSize() - bit - 1, target_str[bit] == '1');
            }
          }
          success = true;
          break;
        }
      }
      return success;
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

    emp::NKLandscape & GetLandscape() { return landscape; }

    void RandomizeLandscapeBits(emp::Random & rand, size_t cnt) {
      landscape.RandomizeStates(rand, cnt);
    }

    void PrintLandscape(std::ostream & out=std::cout) {
      out << emp::to_string(landscape.GetLandscape());
    }

    /// Load NK landscape from file (specified by given path)
    /// First non-commented line of file should give what emp::to_string would output.
    ///  - [ [ fitness fitness fitness ... ] [ fitness fitness fitness ... ] ... ]
    ///  - Loaded landscape must be consistent with num_genes and gene_size
    bool LoadLandscape(const std::string & path) {
      std::ifstream environment_fstream(path);
      if (!environment_fstream.is_open()) {
        std::cout << "Failed to open environment file (" << path << "). Exiting..." << std::endl;
        exit(-1);
      }
      std::string cur_line;
      emp::vector<std::string> line_components;
      bool success = false;
      while (!environment_fstream.eof()) {
        std::getline(environment_fstream, cur_line);
        emp::left_justify(cur_line); // Remove leading whitespace
        if (cur_line == emp::empty_string()) continue; // Skip empty lines
        else if (cur_line[0] == '#') continue; // Treat '#' as a commented line
        else if (cur_line[0] == '[') {
          // Attempt to read landscape state.
          emp::remove_chars(cur_line, "["); // Remove opening brackets. Using closing brackets to delimit gene-associated states.
          emp::left_justify(cur_line);
          emp::right_justify(cur_line);
          line_components.clear();
          emp::slice(cur_line, line_components, ']');
          // Last line component is blank space.
          line_components.pop_back();
          if (landscape.GetN() != line_components.size()) return false;
          for (size_t n = 0; n < line_components.size(); ++n) {
            std::string & values_str = line_components[n];
            emp::vector<std::string> values;
            emp::left_justify(values_str);
            emp::right_justify(values_str);
            emp::slice(values_str, values, ' ');
            if (landscape.GetStateCount() != values.size()) return false;
            for (size_t state = 0; state < values.size(); ++state) {
              double value = emp::from_string<double>(values[state]);
              landscape.SetState(n, state, value);
            }
          }
          success = true;
          break;
        }
      }
      return success;
    }
  };

protected:
  size_t TOTAL_GENS;
  size_t CUR_CHANGE_MAGNITUDE;
  size_t CUR_CHANGE_FREQUENCY;
  size_t CUR_TOURNAMENT_SIZE;
  double CUR_GENE_MOVE_PROB;
  double CUR_BIT_FLIP_PROB;
  double CUR_BIT_INS_PROB;
  double CUR_BIT_DEL_PROB;
  size_t cur_phase=0;

  config_t & config;    ///< World configuration.
  std::string output_path;
  bool setup=false;

  emp::Ptr<NKFitnessModel> fitness_model_nk;
  emp::Ptr<GradientFitnessModel> fitness_model_gradient;
  std::function<void(org_t &)> evaluate_org;
  std::function<void()> change_environment;
  std::function<void()> randomize_environment;
  std::function<bool(const std::string &)> load_environment_from_file;

  emp::Signal<void(size_t)> after_eval_sig; ///< Triggered after organism (ID given by size_t argument) evaluation.

  emp::Ptr<AagosMutator> mutator;

  emp::Ptr<systematics_t> sys_ptr; ///< Shortcut pointer to the correctly-typed systematics manager.
                                   ///< NOTE: The base world class will be responsible for memory management.

  // todo - data collection
  emp::DataManager<double, emp::data::Log, emp::data::Stats, emp::data::Pull> manager;
  emp::Ptr<emp::DataFile> gene_stats_file;
  emp::Ptr<emp::DataFile> representative_org_file;
  emp::Ptr<emp::DataFile> env_file;

  size_t gene_mask;
  size_t most_fit_id;

  void InitFitnessEval();
  void InitEnvironment(); // todo - test environment change
  void InitPop();
  void InitDataTracking();

  void InitLocalConfigs();     ///< Localize paramters that may change for phase two.
  void ActivateEvoPhaseTwo();  ///< Do all the work to transition world into phase two of evolution

  void SetupStatsFile();
  void SetupRepresentativeFile();
  void SetupEnvironmentFile();
  void SetupSystematics();
  void DoPopulationSnapshot();
  void DoConfigSnapshot();
  // TODO - setup environment tracking file?

  /// Shortcut for computing organism's coding sites.
  size_t ComputeCodingSites(org_t & org) {
    size_t count = 0;
    const emp::vector<size_t> & bins = org.GetGeneOccupancyHistogram().GetHistCounts();
    for (size_t i = 1; i < bins.size(); ++i) {
      count += bins[i];
    }
    return count;
  }

  /// Short cut for computing organism's neutral sites.
  size_t ComputeNeutralSites(org_t & org) {
    return org.GetGeneOccupancyHistogram().GetHistCount(0);
  }

public:
  // AagosWorld(emp::Random & random, config_t & cfg) : base_t(random), config(cfg) { Setup(); }
  AagosWorld(config_t & cfg) : config(cfg) { }

  ~AagosWorld() {
    if (config.GRADIENT_MODEL()) fitness_model_gradient.Delete();
    else fitness_model_nk.Delete();
    mutator.Delete();
    representative_org_file.Delete();
    gene_stats_file.Delete();
    env_file.Delete();
  }

  /// Advance world by a single time step (generation).
  void RunStep(bool auto_advance=true);

  ///
  void AdvanceWorld();

  /// Run world for configured number of generations.
  void Run();

  void Setup();

  size_t GetMostFitID() const { return most_fit_id; }
  bool IsSetup() const { return setup; }
  const config_t & GetConfig() const { return config; }

  const NKFitnessModel & GetNKFitnessModel() const { emp_assert(!config.GRADIENT_MODEL()); return *fitness_model_nk; }
  const GradientFitnessModel & GetGradientFitnessModel() const { emp_assert(config.GRADIENT_MODEL()); return *fitness_model_gradient; }

};

// auto update is a concession to the web interface...
void AagosWorld::RunStep(bool auto_advance/*=true*/) {
  emp_assert(setup);
  // (1) evaluate population, (2) select parents, (3) update the world
  // == Do evaluation ==
  most_fit_id = 0;
  for (size_t org_id = 0; org_id < this->GetSize(); ++org_id) {
    emp_assert(IsOccupied(org_id));
    // std::cout << "-- Evaluating org_id " << org_id << " --" << std::endl;
    evaluate_org(GetOrg(org_id));
    if (CalcFitnessID(org_id) > CalcFitnessID(most_fit_id)) {
      most_fit_id = org_id;
    }
    after_eval_sig.Trigger(org_id); // Record phenotype information for this organism's taxon.
    // NOTE - if we wanted to add phenotype tracking to systematics, here's where we could intercept
    //        the necessary information.
  }

  // == Do selection ==
  // if (config.ELITE_COUNT()) emp::EliteSelect(*this, config.ELITE_COUNT(), 1);
  // Run a tournament for the rest...
  // emp::TournamentSelect(*this, config.TOURNAMENT_SIZE(), config.POP_SIZE() - config.ELITE_COUNT());
  emp::TournamentSelect(*this, CUR_TOURNAMENT_SIZE, config.POP_SIZE());

  // == Do update ==
  // If it's a generation to print to console, do so
  const size_t u = GetUpdate();
  if (u % config.PRINT_INTERVAL() == 0) {
    std::cout << u
              << ": max fitness=" << CalcFitnessID(most_fit_id)
              << "; size=" << GetOrg(most_fit_id).GetNumBits();
              // << "; genome=";
    // GetOrg(most_fit_id).Print();
    std::cout << std::endl;
  }

  // Handle managed output files.
  if (config.SUMMARY_INTERVAL()) {
    if ( !(u % config.SUMMARY_INTERVAL()) || (u == config.MAX_GENS()) || (u == TOTAL_GENS) ) {
      gene_stats_file->Update();
      representative_org_file->Update();
    }
  }
  if (config.SNAPSHOT_INTERVAL()) {
    if ( !(u % config.SNAPSHOT_INTERVAL()) || (u == config.MAX_GENS()) ||  (u == TOTAL_GENS) ) {
      DoPopulationSnapshot();
      if (u) sys_ptr->Snapshot(output_path + "phylo_" + emp::to_string(u) + ".csv"); // Don't snapshot phylo at update 0
      env_file->Update();
    }
  }

  if (auto_advance) {
    AdvanceWorld();   // Web interface needs to manage when world update gets called...
  }

}

void AagosWorld::AdvanceWorld() {
  // Should the environment change?
  if (CUR_CHANGE_FREQUENCY) {
    if (!(GetUpdate() % CUR_CHANGE_FREQUENCY)) {
      change_environment();
    }
  }
  Update();
  ClearCache();
}

void AagosWorld::Run() {
  emp_assert(setup);
  for (size_t gen = 0; gen <= config.MAX_GENS(); ++gen) {
    RunStep();
  }
  // Transition?
  if (!config.PHASE_2_ACTIVE()) return;
  // Transition run into phase 2 of evolution
  ActivateEvoPhaseTwo();
  // Run phase of evolution
  for (size_t gen = 0; gen <= config.PHASE_2_MAX_GENS(); ++gen) {
    RunStep();
  }
}

// todo - make callable multiple times?
void AagosWorld::Setup() {
  std::cout << "-- Setting up AagosWorld -- " << std::endl;

  Reset(); // Reset the world

  // Reset world's random number seed.
  random_ptr->ResetSeed(config.SEED());

  // Asserts
  emp_assert(config.NUM_GENES() > 0);

  // Localize phase-one-specific configs
  InitLocalConfigs();

  // Basic setup
  gene_mask = emp::MaskLow<size_t>(config.GENE_SIZE());
  most_fit_id = 0;
  output_path = config.DATA_FILEPATH();
  SetPopStruct_Mixed(true);

  // Initialize fitness evaluation.
  std::cout << "Setting up fitness evaluation." << std::endl;
  InitFitnessEval();
  InitEnvironment();

  // Configure mutator
  std::cout << "Constructing mutator..." << std::endl;
  if (mutator != nullptr) mutator.Delete();
  mutator = emp::NewPtr<AagosMutator>(config.NUM_GENES(), emp::Range<size_t>(config.MIN_SIZE(), config.MAX_SIZE()),
                                      CUR_GENE_MOVE_PROB, CUR_BIT_FLIP_PROB,
                                      CUR_BIT_INS_PROB, CUR_BIT_DEL_PROB);
  std::cout << "  ...done constructing mutator." << std::endl;
  SetMutFun([this](org_t & org, emp::Random & rnd) {
    // NOTE - here's where we would intercept mutation-type distributions (with some extra infrastructure
    //        built into the mutator)!
    org.ResetMutations();
    mutator->ResetLastMutationTracker();
    const size_t mut_cnt = mutator->ApplyMutations(org, rnd);
    auto & mut_dist = mutator->GetLastMutations();
    auto & org_mut_tracker = org.GetMutations();
    org_mut_tracker["bit_flips"] = mut_dist[mutator_t::MUTATION_TYPES::BIT_FLIPS];
    org_mut_tracker["bit_insertions"] = mut_dist[mutator_t::MUTATION_TYPES::BIT_INSERTIONS];
    org_mut_tracker["bit_deletions"] = mut_dist[mutator_t::MUTATION_TYPES::BIT_DELETIONS];
    org_mut_tracker["gene_moves"] = mut_dist[mutator_t::MUTATION_TYPES::GENE_MOVES];
    return mut_cnt;
  });

  // Configure data tracking
  if (!setup) InitDataTracking();

  // Initialize population
  std::cout << "Initialize the population" << std::endl;
  InitPop();

  // Configure world to auto-mutate organisms (if id > elite count)
  // - mutations occur on_before_placement (right before organism added to systematics)
  // SetAutoMutate(config.ELITE_COUNT());
  SetAutoMutate();

  DoConfigSnapshot(); // Snapshot run settings

  setup = true;
}

// todo - add total_gens to config snapshot
void AagosWorld::InitLocalConfigs() {
  CUR_CHANGE_MAGNITUDE = config.CHANGE_MAGNITUDE();
  CUR_CHANGE_FREQUENCY = config.CHANGE_FREQUENCY();
  CUR_GENE_MOVE_PROB = config.GENE_MOVE_PROB();
  CUR_BIT_FLIP_PROB = config.BIT_FLIP_PROB();
  CUR_BIT_INS_PROB = config.BIT_INS_PROB();
  CUR_BIT_DEL_PROB = config.BIT_DEL_PROB();
  CUR_TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();

  TOTAL_GENS = (config.PHASE_2_ACTIVE()) ?  config.MAX_GENS() + config.PHASE_2_MAX_GENS() : config.MAX_GENS();
  cur_phase = 0;
}

void AagosWorld::ActivateEvoPhaseTwo() {
  std::cout << "==> Transitioning to evolution phase two <==" << std::endl;
  // todo

  // Update localized configs as appropriate.
  CUR_CHANGE_MAGNITUDE = config.PHASE_2_CHANGE_MAGNITUDE();
  CUR_CHANGE_FREQUENCY = config.PHASE_2_CHANGE_FREQUENCY();
  CUR_GENE_MOVE_PROB = config.PHASE_2_GENE_MOVE_PROB();
  CUR_BIT_FLIP_PROB = config.PHASE_2_BIT_FLIP_PROB();
  CUR_BIT_INS_PROB = config.PHASE_2_BIT_INS_PROB();
  CUR_BIT_DEL_PROB = config.PHASE_2_BIT_DEL_PROB();
  CUR_TOURNAMENT_SIZE = config.PHASE_2_TOURNAMENT_SIZE();

  // Destruct and re-make mutator for phase two. No need to change the world's mutation function because
  // we're still using the same mutator pointer.
  emp_assert(mutator != nullptr);
  mutator.Delete();
  std::cout << "  Constructing mutator..." << std::endl;
  mutator = emp::NewPtr<AagosMutator>(config.NUM_GENES(), emp::Range<size_t>(config.MIN_SIZE(), config.MAX_SIZE()),
                                      CUR_GENE_MOVE_PROB, CUR_BIT_FLIP_PROB,
                                      CUR_BIT_INS_PROB, CUR_BIT_DEL_PROB);
  std::cout << "    ...done constructing mutator." << std::endl;

  if (config.PHASE_2_LOAD_ENV_FROM_FILE()) {
    // Load the environment from a file.
    const bool success = load_environment_from_file(config.PHASE_2_ENV_FILE());
    if (!success) {
      std::cout << "Failed to load environment from file (" << config.PHASE_2_ENV_FILE() << "). Exiting..." << std::endl;
      exit(-1);
    }
  } else {
    // Randomize the environment.
    randomize_environment();
  }


  ++cur_phase;
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
    if (fitness_model_gradient != nullptr) fitness_model_gradient.Delete();
    fitness_model_gradient = emp::NewPtr<GradientFitnessModel>(*random_ptr, config.NUM_GENES(), config.GENE_SIZE());
    // Print out the gene targets
    std::cout << "Initial gene targets:" << std::endl;
    const auto & targets = fitness_model_gradient->targets;
    for (size_t i = 0; i < targets.size(); ++i) {
      std::cout << "  Target " << i << ": ";
      targets[i].Print();
      std::cout << std::endl;
    }
    // Configure the organism evaluation function.
    evaluate_org = [this](org_t & org) {
      const size_t num_genes = config.NUM_GENES();
      // const size_t num_bits = config.NUM_BITS();
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
        // const size_t tail_bits = num_bits - gene_start; // original?
        const size_t tail_bits = org.GetBits().GetSize() - gene_start; // fix?

        // If a gene runs off the end of the bitstring, loop around to the beginning.
        if (tail_bits < gene_size) {
          gene_val |= (org.GetBits().GetUIntAtBit(0) << tail_bits) & gene_mask;
        }
        // Compute fitness contribution of this gene
        // - Remember, we assume the first index of gene_starts maps to the first index of the target bitstring.
        emp_assert(gene_starts.size() == fitness_model_gradient->targets.size());
        const emp::BitVector & target = fitness_model_gradient->GetTarget(gene_id);
        emp::BitVector gene = emp::BitVector(gene_size);
        gene.SetUIntAtBit(0, gene_val);
        const double fitness_contribution = (double)target.EQU(gene).count() / (double)gene_size;
        phen.gene_fitness_contributions[gene_id] = fitness_contribution;
        fitness += fitness_contribution;
      }
      phen.fitness = fitness;
      // phen.coding_sites = ComputeCodingSites(org);
      // phen.neutral_sites = ComputeNeutralSites(org);
    };
  } else {
    std::cout << "Initializing NK model of fitness." << std::endl;
    if (fitness_model_nk != nullptr) fitness_model_nk.Delete();
    fitness_model_nk = emp::NewPtr<NKFitnessModel>(*random_ptr, config.NUM_GENES(), config.GENE_SIZE());
    // Configure the organism evaluation function.
    evaluate_org = [this](org_t & org) {
      const size_t num_genes = config.NUM_GENES();
      // const size_t num_bits = config.NUM_BITS();
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
        // const size_t tail_bits = num_bits - gene_start; // original?
        const size_t tail_bits = org.GetBits().GetSize() - gene_start; // fix?

        // If a gene runs off the end of the bitstring, loop around to the beginning.
        if (tail_bits < gene_size) {
          gene_val |= (org.GetBits().GetUIntAtBit(0) << tail_bits) & gene_mask;
        }
        // Compute fitness contribution of this gene using nk landscape
        const double fitness_contribution = fitness_model_nk->GetLandscape().GetFitness(gene_id, gene_val);
        phen.gene_fitness_contributions[gene_id] = fitness_contribution;
        fitness += fitness_contribution;
      }
      phen.fitness = fitness;
      // phen.coding_sites = ComputeCodingSites(org);
      // phen.neutral_sites = ComputeNeutralSites(org);
    };
  }
  // Note that this assumes that this organism has been evaluated.
  SetFitFun([](org_t & org) {
    return org.GetPhenotype().fitness;
  });
}

void AagosWorld::InitEnvironment() {
  if (config.GRADIENT_MODEL()) {
    // Configure environment change for gradient fitness model.
    change_environment = [this]() {
      fitness_model_gradient->RandomizeTargetBits(*random_ptr, CUR_CHANGE_MAGNITUDE);
    };
    randomize_environment = [this]() {
       fitness_model_gradient->RandomizeTargets(*random_ptr, config.NUM_GENES());
    };
    load_environment_from_file = [this](const std::string & path) {
      return fitness_model_gradient->LoadTargets(path);
    };
  } else {
    // Configure environment change for nk landscape fitness model.
    change_environment = [this]() {
      fitness_model_nk->RandomizeLandscapeBits(*random_ptr, CUR_CHANGE_MAGNITUDE);
    };
    randomize_environment = [this]() {
      fitness_model_nk->GetLandscape().Reset(*random_ptr);
    };
    load_environment_from_file = [this](const std::string & path) {
      return fitness_model_nk->LoadLandscape(path);
    };
  }
}

void AagosWorld::InitDataTracking() {
  // Create output directory
  #ifndef EMSCRIPTEN
  mkdir(output_path.c_str(), ACCESSPERMS);
  #endif

  if(output_path.back() != '/') {
      output_path += '/';
  }

  SetupFitnessFile(output_path + "fitness.csv").SetTimingRepeat(config.SUMMARY_INTERVAL());
  SetupStatsFile();
  SetupRepresentativeFile();
  SetupEnvironmentFile();
  SetupSystematics();
}

/// Setup data tracking nodes for general statistics about the population.
void AagosWorld::SetupStatsFile() {
  // emp::DataFile & gene_stats_file = SetupFile(output_path + "gene_stats.csv");
  gene_stats_file = emp::NewPtr<emp::DataFile>(output_path + "gene_stats.csv");
  gene_stats_file->AddVar(update, "update", "current generation");
  gene_stats_file->AddVar(cur_phase, "evo_phase", "Current phase of evolution");

  // data node to track number of neutral sites
  // num neutral sites is the size of 0 bin for each org
  auto & neutral_sites_node = manager.New("neutral_sites");
  neutral_sites_node.AddPullSet([this]() {
    emp::vector<double> pop_neut;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_neut.emplace_back(ComputeNeutralSites(*org_ptr));
    }
    return pop_neut;
  });
  // data node to track number of single gene sites
  // size of 1 bin for each org
  auto & single_gene_sites_node = manager.New("single_gene_sites");
  single_gene_sites_node.AddPullSet([this]() {
    emp::vector<double> pop_one;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_one.emplace_back(org_ptr->GetGeneOccupancyHistogram().GetHistCount(1));
    }
    return pop_one;
  });

  // todo - add tracking for each value in histogram?

  // node to track number of multiple overlap sites
  // all bins of size > 1
  auto & multi_gene_sites_node = manager.New("multi_gene_sites");
  multi_gene_sites_node.AddPullSet([this]() {
    emp::vector<double> pop_multi;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      size_t count = 0;
      const emp::vector<size_t> & bins = org_ptr->GetGeneOccupancyHistogram().GetHistCounts(); // get all bins
      for (size_t i = 2; i < bins.size(); ++i) {
        count += bins[i]; // assuming bins are in order, sum all bins
      }
      pop_multi.emplace_back(count);
    }
    return pop_multi;
  });

  // Node to track the number of sites with at least one gene corresponding to it
  auto & coding_sites_node = manager.New("coding_sites");
  coding_sites_node.AddPullSet([this]() {
    emp::vector<double> pop_coding;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_coding.emplace_back(ComputeCodingSites(*org_ptr));
    }
    return pop_coding;
  });

  // Node to track the gene length of each organism
  auto & genome_len_node = manager.New("genome_length");
  genome_len_node.AddPullSet([this]() {
    emp::vector<double> pop_len;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_len.emplace_back(org_ptr->GetNumBits());
    }
    return pop_len;
  });

  // Avg occupancy => average number of genes per site
  auto & occupancy_node = manager.New("avg_occupancy");
  occupancy_node.AddPullSet([this]() {
    emp::vector<double> pop_occupancy;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_occupancy.emplace_back(org_ptr->GetGeneOccupancyHistogram().GetMean());
    }
    return pop_occupancy;
  });

  // Avg gene neighbors
  auto & neighbor_node = manager.New("avg_num_neighbors");
  neighbor_node.AddPullSet([this]() {
    emp::vector<double> pop_neighbor;
    for (emp::Ptr<org_t> org_ptr : pop) {
      if (!org_ptr) continue;
      pop_neighbor.emplace_back(emp::Mean(org_ptr->GetGeneNeighbors()));
    }
    return pop_neighbor;
  });

  // Add all data nodes to the stats file
  gene_stats_file->AddStats(neutral_sites_node, "neutral_sites", "sites with no genes associated with them", true, true);
  gene_stats_file->AddStats(single_gene_sites_node, "single_gene_sites", "sites with exactly one gene associated with them", true, true);
  gene_stats_file->AddStats(multi_gene_sites_node, "multi_gene_sites", "sites with more thone one genes associated with them", true, true);
  gene_stats_file->AddStats(occupancy_node, "site_occupancy", "Average number of genes occupying each site", true, true);
  gene_stats_file->AddStats(neighbor_node, "neighbor_genes", "Average number of other genes each gene overlaps with", true, true);
  gene_stats_file->AddStats(coding_sites_node, "coding_sites", "Number of genome sites with at least one corresponding gene", true, true);
  gene_stats_file->AddStats(genome_len_node, "genome_length", "Length of genome", true, true);

  gene_stats_file->PrintHeaderKeys();
}

/// Setup data tracking for representative organism
void AagosWorld::SetupRepresentativeFile() {
  // emp::DataFile & representative_file = SetupFile(output_path + "representative_org.csv");
  representative_org_file = emp::NewPtr<emp::DataFile>(output_path + "representative_org.csv");
  representative_org_file->AddVar(update, "update", "Current generation");
  representative_org_file->AddVar(cur_phase, "evo_phase", "Current phase of evolution");
  const size_t num_genes = config.NUM_GENES();

  // Fitness
  std::function<double()> fitness_fun = [this]() {
    return CalcFitnessID(most_fit_id);
  };
  representative_org_file->AddFun(fitness_fun, "fitness", "Organism fitness (at this update)");

  // Genome length
  std::function<size_t()> genome_length_fun = [this]() {
    const org_t & org = GetOrg(most_fit_id);
    return org.GetNumBits();
  };
  representative_org_file->AddFun(genome_length_fun, "genome_length", "How many bits in genome?");

  // Number of coding sites for representative organism.
  std::function<size_t()> coding_sites_fun = [this]() {
    org_t & org = GetOrg(most_fit_id);
    return ComputeCodingSites(org);
  };
  representative_org_file->AddFun(coding_sites_fun, "coding_sites", "How many sites in this organism's genome are coding?");

  // Number of neutral sites
  std::function<size_t()> neutral_sites_fun = [this]() {
    org_t & org = GetOrg(most_fit_id);
    return ComputeNeutralSites(org);
  };
  representative_org_file->AddFun(neutral_sites_fun, "neutral_sites", "How many sites in this organim's genome are neutral?");


  // Gene start locations in representative org
  std::function<std::string()> gene_starts_fun = [this]() {
    std::ostringstream stream;
    stream << "\"[";
    const org_t & org = GetOrg(most_fit_id);
    for (size_t i = 0; i < org.GetGeneStarts().size(); ++i) {
      if (i) stream << ",";
      stream << org.GetGeneStarts()[i];
    }
    stream << "]\"";
    return stream.str();
  };
  representative_org_file->AddFun(gene_starts_fun, "gene_starts", "Starting positions for each gene");

  // Organism genome bits
  std::function<std::string()> genome_bits_fun = [this]() {
    std::ostringstream stream;
    const org_t & org = GetOrg(most_fit_id);
    org.GetGenome().bits.Print(stream);
    return stream.str();
  };
  representative_org_file->AddFun(genome_bits_fun, "genome_bitstring", "Bitstring component of genome");

  // Organism gene size.
  std::function<size_t()> genome_gene_size = [this]() {
    const org_t & org = GetOrg(most_fit_id);
    return org.GetGenome().GetGeneSize();
  };
  representative_org_file->AddFun(genome_gene_size, "gene_size", "How many bits is each gene?");

  // Per-gene neighbors
  std::function<std::string()> per_gene_neighbors_fun = [this]() {
    org_t & org = GetOrg(most_fit_id);
    const auto & gene_neighbors = org.GetGeneNeighbors();
    std::ostringstream stream;
    stream << "\"[";
    for (size_t i = 0; i < gene_neighbors.size(); ++i) {
      if (i) stream << ",";
      stream << gene_neighbors[i];
    }
    stream << "]\"";
    return stream.str();
  };
  representative_org_file->AddFun(per_gene_neighbors_fun, "gene_neighbors", "Per-gene neighbors");

  // Mean per-gene neighbors
  std::function<double()> mean_gene_neighbors = [this]() {
    org_t & org = GetOrg(most_fit_id);
    const auto & gene_neighbors = org.GetGeneNeighbors();
    return emp::Mean(gene_neighbors);
  };
  representative_org_file->AddFun(mean_gene_neighbors, "avg_gene_neighbors", "Average per-gene neighbors");

  // Per-site gene occupancy counts
  // For each level of site occupancy, add function that returns the number of sites with that occupancy level.
  for (size_t i = 0; i < num_genes + 1; ++i) {
    std::function<double()> gene_occupancy_fun = [this, i]() {
      emp_assert(i < GetOrg(most_fit_id).GetGeneOccupancyHistogram().GetHistCounts().size());
      emp_assert(most_fit_id < this->GetSize());
      emp_assert(IsOccupied(most_fit_id));
      return (double)GetOrg(most_fit_id).GetGeneOccupancyHistogram().GetHistCount(i);
    };
    representative_org_file->AddFun(gene_occupancy_fun, "site_cnt_" + emp::to_string(i) + "_gene_occupancy", "The number of sites with a particular occupancy level.");
  }

  // representative_file.SetTimingRepeat(config.SUMMARY_INTERVAL());
  representative_org_file->PrintHeaderKeys();

}

void AagosWorld::SetupEnvironmentFile() {
  // environment file should get updated at every snapshot/summary interval
  env_file = emp::NewPtr<emp::DataFile>(output_path + "environment.csv");
  env_file->AddVar(update, "update", "Current generation");
  env_file->AddVar(cur_phase, "evo_phase", "Current phase of evolution");

  std::function<std::string()> get_env_state;
  if (config.GRADIENT_MODEL()) {
    get_env_state = [this]() {
      std::ostringstream stream;
      stream << "\"";
      fitness_model_gradient->PrintTargets(stream);
      stream << "\"";
      return stream.str();
    };
  } else {
    get_env_state = [this]() {
      std::ostringstream stream;
      stream << "\"";
      fitness_model_nk->PrintLandscape(stream);
      stream << "\"";
      return stream.str();
    };
  }
  env_file->AddFun(get_env_state, "env_state", "Current state of the environment");
  env_file->PrintHeaderKeys();
}

void AagosWorld::SetupSystematics() {
  sys_ptr = emp::NewPtr<systematics_t>([](const org_t & o) { return o.GetGenome(); });
  // We want to record phenotype information immediately after an organism is evaluated.
  after_eval_sig.AddAction([this](size_t pop_id) {
    emp::Ptr<taxon_t> taxon = sys_ptr->GetTaxonAt(pop_id);
    taxon->GetData().RecordFitness(this->CalcFitnessID(pop_id));
    taxon->GetData().RecordPhenotype(this->GetOrg(pop_id).GetPhenotype());
  });
  // We want to record mutations when an organism is added to the population
  // - because mutations are applied automatically by this->DoBirth => this->AddOrgAt => sys->OnNew
  std::function<void(emp::Ptr<taxon_t>, org_t&)> record_taxon_mut_data =
    [](emp::Ptr<taxon_t> taxon, org_t & org) {
      taxon->GetData().RecordMutation(org.GetMutations()); // TODO - add mutation tracking to organism!
    };
  sys_ptr->OnNew(record_taxon_mut_data); // Mutations safely happen right before this is triggered
  // Add snapshot functions
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
    return emp::to_string(taxon.GetData().GetFitness());
  }, "mean_fitness", "Taxon fitness");
  // - coding sites
  // sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
  //   return emp::to_string(taxon.GetData().GetPhenotype().coding_sites);
  // }, "coding_sites", "Number of coding sites in taxon genotype.");
  // // - neutral sites
  // sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
  //   return emp::to_string(taxon.GetData().GetPhenotype().neutral_sites);
  // }, "neutral_sites", "Number of neutral sites in taxon genotype.");
  // - mutations from parent
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) -> std::string {
    if (taxon.GetData().HasMutationType("gene_moves")) {
      return emp::to_string(taxon.GetData().GetMutationCount("gene_moves"));
    } else {
      return "0";
    }
  }, "gene_move_muts", "Mutation count");
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) -> std::string {
    if (taxon.GetData().HasMutationType("bit_flips")) {
      return emp::to_string(taxon.GetData().GetMutationCount("bit_flips"));
    } else {
      return "0";
    }
  }, "bit_flip_muts", "Mutation count");
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) -> std::string {
    if (taxon.GetData().HasMutationType("bit_insertions")) {
      return emp::to_string(taxon.GetData().GetMutationCount("bit_insertions"));
    } else {
      return "0";
    }
  }, "bit_ins_muts", "Mutation count");
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) -> std::string {
    if (taxon.GetData().HasMutationType("bit_deletions")) {
      return emp::to_string(taxon.GetData().GetMutationCount("bit_deletions"));
    } else {
      return "0";
    }
  }, "bit_del_muts", "Mutation count");
  // - genome length
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
    return emp::to_string(taxon.GetInfo().bits.GetSize());
  }, "genome_length", "Number of bits in taxon genotype.");
  // - gene starts
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
    const genome_t & taxon_genome = taxon.GetInfo();
    std::ostringstream stream;
    stream << "\"[";
    for (size_t i = 0; i < taxon_genome.gene_starts.size(); ++i) {
      if (i) stream << ",";
      stream << taxon_genome.gene_starts[i];
    }
    stream << "]\"";
    return stream.str();
  }, "gene_starts", "Starting position of each gene.");
  // - genome
  sys_ptr->AddSnapshotFun([](const taxon_t & taxon) {
    std::ostringstream stream;
    taxon.GetInfo().bits.Print(stream);
    return stream.str();
  }, "genome_bitstring", "Bitstring component of taxon genotype.");

  AddSystematics(sys_ptr);
  SetupSystematicsFile(0, output_path + "systematics.csv").SetTimingRepeat(config.SUMMARY_INTERVAL());
}

/// Setup population snapshotting
void AagosWorld::DoPopulationSnapshot() {
  emp::DataFile snapshot_file(output_path + "pop_" + emp::to_string((int)GetUpdate()) + ".csv");
  const size_t num_genes = config.NUM_GENES();
  size_t cur_org_id = 0;
  // Add functions
  snapshot_file.AddVar(update, "update", "Current generation");
  snapshot_file.AddVar(cur_phase, "evo_phase", "Current phase of evolution");

  // Organism ID
  std::function<size_t()> org_id_fun = [&cur_org_id]() {
    return cur_org_id;
  };
  snapshot_file.AddFun(org_id_fun, "org_id", "Organism id");

  // Fitness
  std::function<double()> fitness_fun = [this, &cur_org_id]() {
    return CalcFitnessID(cur_org_id);
  };
  snapshot_file.AddFun(fitness_fun, "fitness", "Organism fitness (at this update)");

  // Genome length
  std::function<size_t()> genome_length_fun = [this, &cur_org_id]() {
    const org_t & org = GetOrg(cur_org_id);
    return org.GetNumBits();
  };
  snapshot_file.AddFun(genome_length_fun, "genome_length", "How many bits in genome?");

  // Number of coding sites for representative organism.
  std::function<size_t()> coding_sites_fun = [this, &cur_org_id]() {
    org_t & org = GetOrg(cur_org_id);
    const emp::vector<size_t> & bins = org.GetGeneOccupancyHistogram().GetHistCounts();
    size_t count = 0;
    for (size_t i = 1; i < bins.size(); ++i) {
      count += bins[i];
    }
    return count;
  };
  snapshot_file.AddFun(coding_sites_fun, "coding_sites", "How many sites in this organism's genome are coding?");

  // Number of neutral sites
  std::function<size_t()> neutral_sites_fun = [this, &cur_org_id]() {
    org_t & org = GetOrg(cur_org_id);
    return org.GetGeneOccupancyHistogram().GetHistCount(0);
  };
  snapshot_file.AddFun(neutral_sites_fun, "neutral_sites", "How many sites in this organim's genome are neutral?");

  // Gene start locations in representative org
  std::function<std::string()> gene_starts_fun = [this, &cur_org_id]() {
    const org_t & org = GetOrg(cur_org_id);
    std::ostringstream stream;
    stream << "\"[";
    for (size_t i = 0; i < org.GetGeneStarts().size(); ++i) {
      if (i) stream << ",";
      stream << org.GetGeneStarts()[i];
    }
    stream << "]\"";
    return stream.str();
  };
  snapshot_file.AddFun(gene_starts_fun, "gene_starts", "Starting positions for each gene");

  // Organism genome bits
  std::function<std::string()> genome_bits_fun = [this, &cur_org_id]() {
    const org_t & org = GetOrg(cur_org_id);
    std::ostringstream stream;
    org.GetGenome().bits.Print(stream);
    return stream.str();
  };
  snapshot_file.AddFun(genome_bits_fun, "genome_bitstring", "Bitstring component of genome");

  // Organism gene size.
  std::function<size_t()> genome_gene_size = [this, &cur_org_id]() {
    const org_t & org = GetOrg(cur_org_id);
    return org.GetGenome().GetGeneSize();
  };
  snapshot_file.AddFun(genome_gene_size, "gene_size", "How many bits is each gene?");

  // Per-gene neighbors
  std::function<std::string()> per_gene_neighbors_fun = [this, &cur_org_id]() {
    org_t & org = GetOrg(cur_org_id);
    const auto & gene_neighbors = org.GetGeneNeighbors();
    std::ostringstream stream;
    stream << "\"[";
    for (size_t i = 0; i < gene_neighbors.size(); ++i) {
      if (i) stream << ",";
      stream << gene_neighbors[i];
    }
    stream << "]\"";
    return stream.str();
  };
  snapshot_file.AddFun(per_gene_neighbors_fun, "gene_neighbors", "Per-gene neighbors");

  // Mean per-gene neighbors
  std::function<double()> mean_gene_neighbors = [this, &cur_org_id]() {
    org_t & org = GetOrg(cur_org_id);
    const auto & gene_neighbors = org.GetGeneNeighbors();
    return emp::Mean(gene_neighbors);
  };
  snapshot_file.AddFun(mean_gene_neighbors, "avg_gene_neighbors", "Average per-gene neighbors");

  // Per-site gene occupancy counts
  // For each level of site occupancy, add function that returns the number of sites with that occupancy level.
  for (size_t i = 0; i < num_genes + 1; ++i) {
    std::function<double()> gene_occupancy_fun = [this, i, &cur_org_id]() {
      emp_assert(i < GetOrg(cur_org_id).GetGeneOccupancyHistogram().GetHistCounts().size());
      emp_assert(cur_org_id < this->GetSize());
      emp_assert(IsOccupied(cur_org_id));
      return (double)GetOrg(cur_org_id).GetGeneOccupancyHistogram().GetHistCount(i);
    };
    snapshot_file.AddFun(gene_occupancy_fun, "site_cnt_" + emp::to_string(i) + "_gene_occupancy", "The number of sites with a particular occupancy level.");
  }

  snapshot_file.PrintHeaderKeys();
  for (cur_org_id = 0; cur_org_id < GetSize(); ++cur_org_id) {
    emp_assert(IsOccupied(cur_org_id));
    snapshot_file.Update();
  }
}

/// Take a snapshot of the configuration settings
void AagosWorld::DoConfigSnapshot() {
  emp::DataFile snapshot_file(output_path + "run_config.csv");
  std::function<std::string()> get_cur_param = []() { return ""; };
  std::function<std::string()> get_cur_value = []() { return ""; };
  snapshot_file.template AddFun<std::string>([&get_cur_param]() -> std::string { return get_cur_param(); }, "parameter");
  snapshot_file.template AddFun<std::string>([&get_cur_value]() -> std::string { return get_cur_value(); }, "value");
  snapshot_file.PrintHeaderKeys();
  get_cur_param = []() { return "TOTAL_GENS"; };
  get_cur_value = [this]() { return emp::to_string(TOTAL_GENS); };
  snapshot_file.Update();
  for (const auto & entry : config) {
    get_cur_param = [&entry]() { return entry.first; };
    get_cur_value = [&entry]() { return emp::to_string(entry.second->GetValue()); };
    snapshot_file.Update();
  }
}


#endif
