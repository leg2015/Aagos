/// @todo Need to add cyclic environments!

#ifndef AAGOS_WORLD_H
#define AAGOS_WORLD_H

#include "Evolve/NK.h"
#include "Evolve/World.h"
#include "tools/math.h"
#include "tools/stats.h"
#include "tools/string_utils.h"
#include <sstream>

#include "AagosOrg.h"

EMP_BUILD_CONFIG(AagosConfig,
                 GROUP(WORLD_STRUCTURE, "How should each organism's genome be setup?"),
                 VALUE(CHANGE_RATE, size_t, 0, "How many changes to fitness tables each generation?"),
                 VALUE(POP_SIZE, size_t, 1000, "How many organisms should be in the population?"),
                 VALUE(MAX_GENS, size_t, 50000, "How many generations should the runs go for?"),
                 VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
                 VALUE(ELITE_COUNT, size_t, 0, "How many organisms should be selected via elite selection?"),
                 VALUE(TOURNAMENT_SIZE, size_t, 2, "How many organisms should be chosen for each tournament?"),

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
                 VALUE(STATISTICS_INTERVAL, size_t, 1000, "How many updates between statistic gathering?"),
                 VALUE(SNAPSHOT_INTERVAL, size_t, 10000, "How many updates between snapshots?"))

class AagosWorld : public emp::World<AagosOrg>
{
private:
  using base_t = emp::World<AagosOrg>;

  AagosConfig &config;
  emp::NKLandscape landscape;
  // need a node manager for data tracking since so many different data points to draw
  emp::DataManager<double, emp::data::Log, emp::data::Stats, emp::data::Pull> manager; 
  emp::Ptr<emp::ContainerDataFile<emp::vector<emp::Ptr<AagosOrg>>>> snapshot_file;

  // Configured values
  size_t num_bits;
  size_t num_genes;
  size_t gene_size;

  // Calculated values
  size_t gene_mask;
  int fittest_id;

public:
  AagosWorld(AagosConfig &_config, const std::string &world_name = "AagosWorld")
      : emp::World<AagosOrg>(world_name)
      , config(_config)
      , landscape(config.NUM_GENES() , config.GENE_SIZE() - 1, GetRandom())
     // , manager()
      , num_bits(config.NUM_BITS())
      , num_genes(config.NUM_GENES())
      , gene_size(config.GENE_SIZE())
      , gene_mask(emp::MaskLow<size_t>(config.GENE_SIZE()))
      , fittest_id(-1) // set to -1 to indicate fittest individual hasn't been calc yet

  {
    emp_assert(config.MIN_SIZE() >= config.GENE_SIZE(), "BitSet can't handle a genome smaller than gene_size");
    // fitness function for aagos orgs
    auto fit_fun = [this](AagosOrg &org) {
      double fitness = 0.0;
      for (size_t gene_id = 0; gene_id < num_genes; gene_id++)
      {
        const size_t gene_pos = org.gene_starts[gene_id];
        size_t gene_val = org.bits.GetUIntAtBit(gene_pos) & gene_mask;
        const size_t tail_bits = num_bits - gene_pos;

        // If a gene runs off the end of the bitstring, loop around to the beginning.
        if (tail_bits < gene_size)
        {
          gene_val |= (org.bits.GetUInt(0) << tail_bits) & gene_mask;
        }
        fitness += landscape.GetFitness(gene_id, gene_val);
      }
      return fitness;
    };
    SetFitFun(fit_fun);

    // Setup the mutation function. Per site.
    std::function<size_t(AagosOrg &, emp::Random &)> mut_fun =
        [this](AagosOrg &org, emp::Random &random) {
          // Do gene moves.
          size_t num_moves = random.GetRandBinomial(org.GetNumGenes(), config.GENE_MOVE_PROB());
          for (size_t m = 0; m < num_moves; m++)
          {
            size_t gene_id = random.GetUInt(org.GetNumGenes());
            org.gene_starts[gene_id] = random.GetUInt(org.GetNumBits());
          }

          // Do bit flips mutations
          size_t num_flips = random.GetRandBinomial(org.GetNumBits(), config.BIT_FLIP_PROB());
          for (size_t m = 0; m < num_flips; m++)
          {
            const size_t pos = random.GetUInt(org.GetNumBits());
            org.bits[pos] ^= 1;
          }

          // Get num of insertions and deletions.
          int num_insert = random.GetRandBinomial(org.GetNumBits(), config.BIT_INS_PROB());
          int num_delete = random.GetRandBinomial(org.GetNumBits(), config.BIT_DEL_PROB());
          const int proj_size = (int)org.bits.GetSize() + num_insert - num_delete;

          // checks gene size is within range
          if (proj_size > config.MAX_SIZE())
          { // if size of genome larger than max, restrict to max size
            num_insert -= proj_size - config.MAX_SIZE();
          }
          else if (proj_size < config.MIN_SIZE())
          { // else if size of genome smaller than min, restrict as well
            num_delete -= config.MIN_SIZE() - proj_size;
          }

          // asserts size limitations
          emp_assert((int)org.bits.GetSize() + num_insert - num_delete >= config.MIN_SIZE(),
                       "the genome size can't be smaller than the genome length, else BitSet breaks");
          emp_assert((int)org.bits.GetSize() + num_insert - num_delete <= config.MAX_SIZE(), "some limit on bloat of program");

          // Do insertions
          for (int i = 0; i < num_insert; i++) // For each insertion that occurs,
          {
            const size_t pos = random.GetUInt(org.GetNumBits()); // Figure out position for insertion.
            org.bits.Resize(org.bits.GetSize() + 1);             // Increase size to make room for insertion.
            emp::BitVector mask(pos, 1);                         // Setup a mask to preserve early bits.
            mask.Resize(org.bits.GetSize());                     // Align mask size.

            // Now build the new string!
            org.bits = (mask & org.bits) | ((org.bits << 1) & ~mask);
            org.bits[pos] = random.P(0.5); // Randomize the new bit.

            // Shift any genes that started at pos or later.
            for (auto &x : org.gene_starts)
              if (x >= pos)
                x++;
          }

          // Do deletions
          for (int i = 0; i < num_delete; i++) // For each deletion that occurs,
          {
            size_t pos = random.GetUInt(org.GetNumBits()); // Figure out position to delete.
            emp::BitVector mask(pos, 1);                   // Setup a mask to preserve early bits.
            mask.Resize(org.bits.GetSize());               // Align mask size.

            org.bits = (mask & org.bits) | ((org.bits >> 1) & ~mask); // Build the new string!
            org.bits.Resize(org.bits.GetSize() - 1);                  // Decrease size to account for deletion

            // Shift any genes that started at pos or later.
            if (pos == 0)
              pos = 1; // Adjust position if beginning was deleted.
            for (auto &x : org.gene_starts)
              if (x >= pos)
                x--;
          }

          return num_moves + num_flips + num_insert + num_delete; // Returns total num mutations
        };
    SetMutFun(mut_fun);       // set mutation function of world to above
    SetPopStruct_Mixed(true); // uses well-mixed population structure
    SetDataTracking();        // sets up data tracking
    
  }

  ~AagosWorld() { ; }

  // Finds fittest individual in the curr population
  void FindFittest()
  {
    if (fittest_id == -1) // only run if fittest individual not yet calculated
    {
      for (size_t i = 0; i < GetSize(); i++)
      {
        if (!pop[i])
          continue;
        if (CalcFitnessID(i) > CalcFitnessID(fittest_id))
          fittest_id = i;
      }
    }
  }

  // gets all orgs that are valid in curr gen
  // given ids of all valid orgs
  emp::vector<emp::Ptr<AagosOrg>> GetValidOrgs(emp::vector<size_t> valid_org_ids)
  {
    emp::vector<emp::Ptr<AagosOrg>> valid_orgs;
    for (size_t curr : valid_org_ids)
    {
      valid_orgs.push_back(pop[curr]);
    }
    return valid_orgs;
  }

  // sets up all data tracking for world
  // includes both snapshots and statistics
  void SetDataTracking()
  {
    SetStatsFile();          // sets up all data files for generals stats
    SetRepresentativeFile(); // sets up all data files for representative pop member (stats runs only)
    SetSnapshotFile();       // sets up all data files for snapshots
  }

  // sets up data tracking nodes for general statistics about population
  void SetStatsFile()
  {
    SetupFitnessFile().SetTimingRepeat(config.STATISTICS_INTERVAL()); // set timing to interval
    auto gene_stats_file = SetupFile("gene_stats.csv");
    gene_stats_file.AddVar(update, "update", "update of current gen"); // tracks which update stats calc on

    // data node to track number of neutral sites
    // num neutral sites is the size of 0 bin for each org
    auto neutral_node = manager.New("neutral_sites");
    neutral_node.AddPullSet([this]() {
      emp::vector<double> pop_neut;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        pop_neut.emplace_back(org->GetHistogram().GetHistCount(0));
      }
      return pop_neut;
    });

    // data node to track number of single gene sites
    // size of 1 bin for each org
    auto one_gene_node = manager.New("one_gene_sites");
    one_gene_node.AddPullSet([this]() {
      emp::vector<double> pop_one;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        pop_one.emplace_back(org->GetHistogram().GetHistCount(1));
      }
      return pop_one;
    });

    // node to track number of multiple overlap sites
    // all bins of size > 1
    auto multi_gene_node = manager.New("multi_gene_sites");
    multi_gene_node.AddPullSet([this]() {
      emp::vector<double> pop_multi;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        int count = 0;
        const emp::vector<size_t> &bins = org->GetHistogram().GetHistCounts(); // get all bins
        for (size_t i = 2; i < bins.size(); i++)                               // check all bins that are > 1
        {
          count += bins[i]; // assuming bins are in order, sum all bins
        }
        pop_multi.emplace_back(count);
      }
      return pop_multi;
    });

    // avg overlap is average number of genes per site
    // measure amount of overlap in a genome
    // calculated as mean of histogram
    auto overlap_node = manager.New("avg_overlap");
    overlap_node.AddPullSet([this]() { //TODO: could addpullset be why node fn not getting triggered ever?
      emp::vector<double> pop_overlap;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        pop_overlap.emplace_back(org->GetHistogram().GetMean());
      }
      return pop_overlap;
    });

    // neighbor node gets mean of gene neighbors for each org
    auto neighbor_node = manager.New("avg_num_neighbors");
    neighbor_node.AddPullSet([this]() {
      emp::vector<double> pop_neighbor;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        pop_neighbor.emplace_back(emp::Mean(org->GetGeneNeighbors()));
      }
      return pop_neighbor;
    });

    // add all data nodes to stats data file
    gene_stats_file.AddStats(neutral_node, "Neutral_Sites", "sites with no genes associated with them", true, true);
    gene_stats_file.AddStats(one_gene_node, "One_Gene_Sites", "sites with exactly one gene associated with them", true, true);
    gene_stats_file.AddStats(multi_gene_node, "Multi_Gene_Sites", "sites with more thone one genes associated with them", true, true);
    gene_stats_file.AddStats(overlap_node, "Overlap", "Average number of genes per site", true, true);
    gene_stats_file.AddStats(neighbor_node, "Neighbor_Genes", "Number of genes overlapping each other gene", true, true);
    // set calc update timing
    gene_stats_file.SetTimingRepeat(config.STATISTICS_INTERVAL());
    gene_stats_file.PrintHeaderKeys();
  }

  // sets up data file for representative org
  void SetRepresentativeFile()
  {

    auto representative_file = SetupFile("representative_org.csv");
    representative_file.AddVar(update, "update", "update of current gen");

    // function param to add_fun in data node is a std::fn, not a lambda fn
    // so all fn had to be moved to outside scope to work properly

    // gets bin value for each bin in histogram of representative org
    std::function<double()> gene_overlap_fun;
    for (size_t b = 0; b < config.NUM_GENES() + 2; b++) // loops through each bin in histogram
    {
      // fn for current bin of histogram
      gene_overlap_fun = [this, b]() {
        FindFittest(); // since order not guaranteed, must look for fittest ind. in each fn call
        return pop[fittest_id]->GetHistogram().GetHistCount(b);
      };
      // add current function to file
      representative_file.AddFun(gene_overlap_fun, emp::to_string(b) + "_gene_overlap_frequency",
           "statistics for representative population member");
    }

    // gets gene start locations for representative org
    std::function<std::string()> gene_starts_fun = [this]() {
      FindFittest();
      return emp::to_string(pop[fittest_id]->GetGeneStarts());
    };
    representative_file.AddFun(gene_starts_fun, "gene_starts", 
          "all gene starts for the representative organism in the population");

    // gets genome size for representative org
    std::function<double()> genome_size_fun = [this]() {
      FindFittest();
      return pop[fittest_id]->GetNumBits();
    };
    representative_file.AddFun(genome_size_fun, "genome_size",
           "genome size of representative organism");

    // gets fitness of representative org
    std::function<double()> fitness_fun = [this]() {
      FindFittest();
      return CalcFitnessID(fittest_id);
    };
    representative_file.AddFun(fitness_fun, "fitness", "fitness of representative org");

    // gets mean of neighbors of rep. org.
    std::function<double()> genome_neighbor_fun = [this]() {
      FindFittest();
      return emp::Mean(pop[fittest_id]->GetGeneNeighbors());
    };
    representative_file.AddFun(genome_neighbor_fun, "gene_neighbors",
           "gene neighbors of representative org");
    // make sure file updated with each statistics pull
    representative_file.SetTimingRepeat(config.STATISTICS_INTERVAL());
    representative_file.PrintHeaderKeys();
  }

  // sets up collection file for each snapshot
  // aggregates all snapshot orgsi into a csv separated by col
  void SetSnapshotFile()
  {

    // fn that gets all valid orgs for snapshot
    std::function<emp::vector<emp::Ptr<AagosOrg>>()> snapshot_fun = [this]() {
      return GetValidOrgs(GetValidOrgIDs()); // gets ids of all valid orgs
    };
    //create snapshot file
    auto temp_file = emp::MakeContainerDataFile(snapshot_fun, "snapshot.csv");
    snapshot_file.New(temp_file);
    // lists which update file created on
    snapshot_file->AddVar(update, "update", "update of current gen");

    // gets full histogram of each org
    std::function<int(emp::Ptr<AagosOrg>)> snap_gene_overlap_fun;
    for (size_t b = 0; b < config.NUM_GENES() + 2; b++) // loop through each bin of hist & get val of bin
    {
      snap_gene_overlap_fun = [this, b](emp::Ptr<AagosOrg> org) {
        return org->GetHistogram().GetHistCount(b);
      };
      // add fn for each bin to file
      snapshot_file->AddContainerFun(snap_gene_overlap_fun, emp::to_string(b) + "_gene_overlap_frequency",
             "histogram statistics for current population member");
    }

    // gets gene start pos of each org as a string
    std::function<std::string(emp::Ptr<AagosOrg>)> snap_gene_starts_fun = [this](emp::Ptr<AagosOrg> org) {
      return emp::to_string(org->GetGeneStarts());
    };
    snapshot_file->AddContainerFun(snap_gene_starts_fun, "gene_starts",
           "all gene starts for the representative organism in the population");

    // gets genome size of each org
    std::function<double(emp::Ptr<AagosOrg>)> snap_genome_size_fun = [this](emp::Ptr<AagosOrg> org) {
      return org->GetNumBits();
    };
    snapshot_file->AddContainerFun(snap_genome_size_fun, "genome_size", "genome size of representative organism");

    // gets fitness of each org
    std::function<double(emp::Ptr<AagosOrg>)> snap_fitness_fun = [this](emp::Ptr<AagosOrg> org) {
      return CalcFitnessOrg(*org);
    };
    snapshot_file->AddContainerFun(snap_fitness_fun, "fitness", "fitness of representative org");

    // gets mean of gene neighbors of each org
    std::function<double(emp::Ptr<AagosOrg>)> snap_genome_neighbor_fun = [this](emp::Ptr<AagosOrg> org) {
      return emp::Mean(org->GetGeneNeighbors());
    };
    snapshot_file->AddContainerFun(snap_genome_neighbor_fun, "gene_neighbors", "gene neighbors of representative org");

    // gets full bitstring genome of each org
    std::function<std::string(emp::Ptr<AagosOrg>)> snap_bitstring_fun = [this](emp::Ptr<AagosOrg> org) {
      std::stringstream bit_tostring;          // need a bitstream to print bitvector
      org->GetBits().PrintArray(bit_tostring); // converts bitvector to string in bistream
      return bit_tostring.str();               // must extract bitstring from bistream
    };
    snapshot_file->AddContainerFun(snap_bitstring_fun, "genome", "genome of current org");
    snapshot_file->SetTimingRepeat(config.SNAPSHOT_INTERVAL());
    snapshot_file->PrintHeaderKeys();
    AddDataFile(snapshot_file);
  }
  // updates world
  void Update()
  {
    landscape.RandomizeStates(GetRandom(), config.CHANGE_RATE());
    base_t::Update();
    fittest_id = -1; // reset fittest id flag
  }
};

#endif
