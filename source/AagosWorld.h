/// @todo Need to add cyclic environments!

#ifndef AAGOS_WORLD_H
#define AAGOS_WORLD_H

#include "Evolve/NK.h"
#include "Evolve/World.h"
#include "tools/Binomial.h"
#include "tools/math.h"
#include "tools/stats.h"
#include "tools/string_utils.h"

#include <sstream>
#include <string>

#include "AagosOrg.h"

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
                 VALUE(STATISTICS_INTERVAL, size_t, 1000, "How many updates between statistic gathering?"),
                 VALUE(SNAPSHOT_INTERVAL, size_t, 10000, "How many updates between snapshots?"),
                 VALUE(DATA_FILEPATH, std::string, "", "what directory should all data files be written to?"))

class AagosWorld : public emp::World<AagosOrg>
{
private:
  using base_t = emp::World<AagosOrg>;

  AagosConfig &config;
  emp::NKLandscape landscape;
  // need a node manager for data tracking since so many different data points to draw
  emp::DataManager<double, emp::data::Log, emp::data::Stats, emp::data::Pull> manager; 
  emp::Ptr<emp::ContainerDataFile<emp::vector<emp::Ptr<AagosOrg>>>> snapshot_file;
  emp::vector<emp::BitVector> target_bits; // vector of target bitstrings for gradient version of model

  // Configured values
  size_t num_bits;
  size_t num_genes;
  size_t gene_size;
  size_t num_bins;
  bool gradient;

  std::string data_filepath;
  emp::Binomial gene_moves_binomial;
  emp::vector<emp::Binomial> bit_flips_binomials;
  emp::vector<emp::Binomial> inserts_binomials;
  emp::vector<emp::Binomial> deletes_binomials;

  // Calculated values
  size_t gene_mask;
  int fittest_id;

public:
  AagosWorld(emp::Random &rand, AagosConfig &_config, const std::string &world_name = "AagosWorld")
      : emp::World<AagosOrg>(rand, world_name), config(_config), landscape(config.NUM_GENES(), config.GENE_SIZE() - 1, GetRandom())
        // , manager()
        ,
        num_bits(config.NUM_BITS()), num_genes(config.NUM_GENES()), gene_size(config.GENE_SIZE()), num_bins(config.NUM_GENES() + 1), data_filepath(config.DATA_FILEPATH()) // TODO: only works if subdir is made before runs start... TODO: wouldn't work if subdir not created, runs wouldn't be stored
        ,
        gene_moves_binomial(config.GENE_MOVE_PROB(), config.NUM_GENES()) // since num genes doesn't evolve, can calculate 1 dist
        ,
        gene_mask(emp::MaskLow<size_t>(config.GENE_SIZE())), fittest_id(-1) // set to -1 to indicate fittest individual hasn't been calc yet
        ,
        gradient(config.GRADIENT_MODEL())

  {
    emp_assert(config.MIN_SIZE() >= config.GENE_SIZE(), "BitSet can't handle a genome smaller than gene_size");
    // for each possible length of genome, calculate the bin dist for that length
    // start at smallest possible gene length
    for (size_t i = config.MIN_SIZE(); i <= config.MAX_SIZE(); i++) {
      bit_flips_binomials.emplace_back(config.BIT_FLIP_PROB(), i);
      inserts_binomials.emplace_back(config.BIT_INS_PROB(), i);
      deletes_binomials.emplace_back(config.BIT_DEL_PROB(), i);
    }

  // if using gradient model, initialize target bitstrings
  if(gradient) {
    for(size_t i = 0; i < num_genes; i++) {
      auto &rand = GetRandom();//TODO: is this bad?
      target_bits.emplace_back(emp::RandomBitVector(rand, gene_size)); 
    }
    // TODO: will break if the number of genes is allowed to evolve ever
    emp_assert(target_bits.size() == num_genes "there should be the same number of target bitstrings as genes in genomes"); 
    }
    // fitness function for aagos orgs
    auto fit_fun = [this](AagosOrg &org) { //TODO: change to prportion of matching bits
      double fitness = 0.0; // TODO: use hamming distance to compare bistrings - UES
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
        // calculate fitness
        if(gradient) { // remember that we're assuming here that 1st index of gene_starts maps to 1st index in target bitstring
          //TODO: need to convert 32 bit unit that gene_val is to a bitvector, going to complete rest of code assuming this
          fitness += target_bits[gene_id].EQU(gene_val).count() / gene_size;
          //calcs hamming dist between target and curr gene & adds up # of matches
          // divide by num bits in gene so fitness range is (0, 1)
          } else {
          fitness += landscape.GetFitness(gene_id, gene_val);
        }
      }
      return fitness;
    };
    SetFitFun(fit_fun);

    // Setup the mutation function. Per site.
    std::function<size_t(AagosOrg &, emp::Random &)> mut_fun =
        [this](AagosOrg &org, emp::Random &random) {
          size_t bin_array_offset = org.GetNumBits() - config.MIN_SIZE(); // offset is num bits - min size of genome
          emp_assert(bin_array_offset >= 0, "index of bin dist cannot be negative!!");
          // Do gene moves.
          size_t num_moves = gene_moves_binomial.PickRandom(random);
          for (size_t m = 0; m < num_moves; m++)
          {
            size_t gene_id = random.GetUInt(org.GetNumGenes()); // get random gene
            org.gene_starts[gene_id] = random.GetUInt(org.GetNumBits()); // change its start to a random location
          }

          size_t num_flips = bit_flips_binomials[bin_array_offset].PickRandom(random);
          for (size_t m = 0; m < num_flips; m++)
          {
            const size_t pos = random.GetUInt(org.GetNumBits());
            org.bits[pos] ^= 1;
          }

          // Get num of insertions and deletions.
          size_t num_insert = inserts_binomials[bin_array_offset].PickRandom(random);
          size_t num_delete = deletes_binomials[bin_array_offset].PickRandom(random);
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

          int num_muts = num_moves + num_flips + num_insert + num_delete;
          if (num_muts > 0) {
            org.ResetHistogram();
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
      fittest_id = 0;
      for (size_t i = 1; i < GetSize(); i++)
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
    SetupFitnessFile(data_filepath + "fitness.csv").SetTimingRepeat(config.STATISTICS_INTERVAL()); // set timing to interval
    SetStatsFile();          // sets up all data files for generals stats
    SetRepresentativeFile(); // sets up all data files for representative pop member (stats runs only)
    SetSnapshotFile();       // sets up all data files for snapshots
  }

  // sets up data tracking nodes for general statistics about population
  void SetStatsFile()
  {
    emp::DataFile & gene_stats_file = SetupFile(data_filepath + "gene_stats.csv");
    gene_stats_file.AddVar(update, "update", "update of current gen"); // tracks which update stats calc on

    // data node to track number of neutral sites
    // num neutral sites is the size of 0 bin for each org
    auto & neutral_node = manager.New("neutral_sites");
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
    auto & one_gene_node = manager.New("one_gene_sites");
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
    auto & multi_gene_node = manager.New("multi_gene_sites");
    multi_gene_node.AddPullSet([this]() {
      emp::vector<double> pop_multi;
      for (emp::Ptr<AagosOrg> org : pop)
      {
        if (!org)
          continue;
        int count = 0;
        const emp::vector<size_t> &bins = org->GetHistogram().GetHistCounts(); // get all bins
        for (size_t i = 2; i < bins.size(); i++)                              // check all bins that are > 1
        {
          count += bins[i]; // assuming bins are in order, sum all bins
        }
        pop_multi.emplace_back(count);
      }
      return pop_multi;
    });
  // node to track the number of sites with at least one gene corresponding to it
  auto & coding_sites_node = manager.New("coding_sites");
  coding_sites_node.AddPullSet([this]() {
    emp::vector<double> pop_coding;
    for ( emp::Ptr<AagosOrg> org : pop) 
    {
      if (!org)
        continue;
      int count = 0;
      const emp::vector<size_t> &bins = org->GetHistogram().GetHistCounts();
      for (size_t i = 1; i < bins.size(); i++) // start with bin corresponding to one gene
      {
        count += bins[i]; // sum all bin counts
      }
      pop_coding.emplace_back(count);
    }
    return pop_coding;
  });
  // node to track the gene length of each organism
  auto & gene_len_node = manager.New("gene_len");
  gene_len_node.AddPullSet([this]() {
        emp::vector<double> pop_len;
    for ( emp::Ptr<AagosOrg> org : pop) 
    {
      if (!org)
        continue;
      int count = 0;
      const emp::vector<size_t> &bins = org->GetHistogram().GetHistCounts();
      for (size_t i = 0; i < bins.size(); i++) // start with bin corresponding to no gene
      {
        count += bins[i]; // sum all bin counts
      }
      pop_len.emplace_back(count);
    }
    return pop_len;
  });
    // avg overlap is average number of genes per site
    // measure amount of overlap in a genome
    // calculated as mean of histogram
    auto & overlap_node = manager.New("avg_overlap");
    overlap_node.AddPullSet([this]() { 
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
    auto & neighbor_node = manager.New("avg_num_neighbors");
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
    gene_stats_file.AddStats(neutral_node, "neutral_sites", "sites with no genes associated with them", true, true);
    gene_stats_file.AddStats(one_gene_node, "one_gene_sites", "sites with exactly one gene associated with them", true, true);
    gene_stats_file.AddStats(multi_gene_node, "multi_gene_sites", "sites with more thone one genes associated with them", true, true);
    gene_stats_file.AddStats(overlap_node, "overlap", "Average number of genes per site", true, true);
    gene_stats_file.AddStats(neighbor_node, "neighbor_genes", "Number of genes overlapping each other gene", true, true);
    gene_stats_file.AddStats(coding_sites_node, "coding_sites", "Number of genome sites with at least one corresponding gene", true, true);
    gene_stats_file.AddStats(gene_len_node, "gene_length", "Length of genome", true, true);
    // set calc update timing
    gene_stats_file.SetTimingRepeat(config.STATISTICS_INTERVAL());
    gene_stats_file.PrintHeaderKeys();
  }

  // sets up data file for representative org
  void SetRepresentativeFile()
  {

    emp::DataFile & representative_file = SetupFile( data_filepath + "representative_org.csv");
    representative_file.AddVar(update, "update", "update of current gen");

    // function param to add_fun in data node is a std::fn, not a lambda fn
    // so all fn had to be moved to outside scope to work properly

    // gets bin value for each bin in histogram of representative org
    // assumes each gene is same length across genome and across organisms
    std::function<double()> gene_overlap_fun;
    for (size_t b = 0; b < num_bins; b++) // loops through each bin in histogram
    {
      // fn for current bin of histogram
      gene_overlap_fun = [this, b]() {
        FindFittest(); // since order not guaranteed, must look for fittest ind. in each fn call
        return pop[fittest_id]->GetHistogram().GetHistCount(b);
      };
      // add current function to file
      representative_file.AddFun(gene_overlap_fun, emp::to_string(b) + "_gene_overlap",
           "statistics for representative population member");
    }

    // gets gene start locations for representative org
    std::function<std::string()> gene_starts_fun = [this]() {
      FindFittest();
      return emp::to_string(pop[fittest_id]->GetGeneStarts());
    };
    representative_file.AddFun(gene_starts_fun, "gene_starts", 
          "all gene starts for the representative organism in the population");

    // gets number of coding sites for representative org
    std::function<double()> coding_sites_fun = [this]() {
        FindFittest();
        const emp::vector<size_t> &bins = pop[fittest_id]->GetHistogram().GetHistCounts();
        int count = 0;
        for (size_t i = 1; i < bins.size(); i++) // start with bin corresponding to one gene
          {
            count += bins[i]; // sum all bin counts
          }
          return count;
    };
    representative_file.AddFun(coding_sites_fun, "coding_sites", 
          "number of coding sites for representative organism");


    // gets genome length for representative org
    std::function<double()> genome_size_fun = [this]() {
      FindFittest();
      return pop[fittest_id]->GetNumBits();
    };
    representative_file.AddFun(genome_size_fun, "genome_length",
           "genome length of representative organism");

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
    auto temp_file = emp::MakeContainerDataFile(snapshot_fun, data_filepath + "snapshot.csv");
    snapshot_file.New(temp_file);
    // lists which update file created on
    snapshot_file->AddVar(update, "update", "update of current gen");

    // gets full histogram of each org
    // assumes each gene is same length across genome and across organisms
    std::function<int(emp::Ptr<AagosOrg>)> snap_gene_overlap_fun;
    for (size_t b = 0; b < num_bins; b++) // loop through each bin of hist & get val of bin
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
    // do environmental change
    if(gradient) {
      for(size_t i = 0; i < config.CHANGE_RATE(); i++) {
        auto &rand = GetRandom();
        // grad a randomly chosen target sequence and assign to a new randomly generated target sequence
        auto to_change = target_bits[rand.GetUInt(target_bits.size())];
        to_change.Set(rand.GetUInt(to_change.size()), !to_change.Get(i)); //TODO: call set function on bitvector  call getBit on bitvector to flip
        // random has a getUInt, and also use getUInt to get the bit within the bit vector
        // bitvector.Set(i, !bitvector.Get(i))
        // bitvector.Set(random->GetUInt(bitvector.size()), !bitvector.Get(i)) its okay b/c bit NEEDS TO FLIP for env change
      }
    } else { // default
      landscape.RandomizeStates(GetRandom(), config.CHANGE_RATE());
    }
    base_t::Update();
    fittest_id = -1; // reset fittest id flag
  }
};

#endif
