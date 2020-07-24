#ifndef AAGOS_ORG_H
#define AAGOS_ORG_H

#include <algorithm>
#include <unordered_map>
#include "tools/BitVector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/string_utils.h"
#include "data/DataNode.h"

class AagosOrg {
public:

  /// AagosOrg genomes comprise a bitsequence and the starting positions of each gene in the bitsequence.
  struct Genome {
    size_t ancestral_id=0; ///< Only used for loaded genomes. Used to identify which loaded genome this genome descends from (when phylo tracking is off).
    size_t gene_size; ///< size of each gene in the genome
    size_t num_genes; ///< number of genes in the genome
    emp::BitVector bits;              ///< Bit sequence
    emp::vector<size_t> gene_starts;  ///< Starting positions of all genes.

    Genome(size_t _num_bits, size_t _num_genes, size_t _gene_size)
      : gene_size(_gene_size), num_genes(_num_genes),
        bits(_num_bits), gene_starts(_num_genes, 0) { }
    Genome(const Genome &) = default;
    Genome(Genome &&) = default;

    bool operator==(const Genome & other) const {
      return std::tie(gene_size, num_genes, bits, gene_starts, ancestral_id)
              == std::tie(other.gene_size, other.num_genes, other.bits, other.gene_starts, other.ancestral_id);
    }

    bool operator!=(const Genome & other) const {
      return !(*this == other);
    }

    bool operator<(const Genome & other) const {
      return std::tie(gene_size, num_genes, bits, gene_starts, ancestral_id)
              < std::tie(other.gene_size, other.num_genes, other.bits, other.gene_starts, other.ancestral_id);
    }

    // Randomize genome and gene starts
    void Randomize(emp::Random & random) {
      emp::RandomizeBitVector(bits, random);
      emp::RandomizeVector<size_t>(gene_starts, random, 0, bits.size());
    }

    size_t GetNumBits() const { return bits.size(); }
    size_t GetGeneSize() const { return gene_size; }
    size_t GetNumGenes() const { return num_genes; }
    size_t GetAncestralID() const { return ancestral_id; }

  };

  struct Phenotype {
    double fitness=0.0;
    emp::vector<double> gene_fitness_contributions;
    bool evaluated=false;
    // -- things that we want to know about phenotype for systematics tracking --
    // - coding_sites
    // size_t coding_sites=0;
    // - neutral sites
    // size_t neutral_sites=0;

    Phenotype(size_t num_genes=0) : gene_fitness_contributions(num_genes, 0.0) { }
    Phenotype(const Phenotype &) = default;
    Phenotype(Phenotype &&) = default;

    Phenotype & operator=(const Phenotype &) = default;

    void Reset() {
      evaluated=false;
      fitness = 0.0;
      // coding_sites=0;
      // neutral_sites=0;
      std::fill(gene_fitness_contributions.begin(), gene_fitness_contributions.end(), 0);
    }

    bool operator==(const Phenotype & other) const {
      return std::tie(fitness, gene_fitness_contributions, evaluated/*, coding_sites, neutral_sites*/)
              == std::tie(other.fitness, other.gene_fitness_contributions, other.evaluated/*, other.coding_sites, other.neutral_sites*/);
    }

    bool operator!=(const Phenotype & other) const {
      return !(*this == other);
    }

    bool operator<(const Phenotype & other) const {
      return std::tie(fitness, gene_fitness_contributions, evaluated/*, coding_sites, neutral_sites*/)
              < std::tie(other.fitness, other.gene_fitness_contributions, other.evaluated/*, other.coding_sites, other.neutral_sites*/);
    }

    bool IsEvaluated() const { return evaluated; }
  };

  using histogram_t = emp::DataNode<int, emp::data::Histogram, emp::data::Stats>;

protected:
  Genome genome;    ///< Genotype
  Phenotype phenotype;

  /// # neighbors (per gene measurement) - the number of neighbors each gene has where a neighbor is
  /// another gene that overlaps the focal gene by at least one bit.
  emp::vector<size_t> gene_neighbors;

  /// Histogram object that stores the number of overlapped genes at each bit in the genome.
  histogram_t occupancy_histogram;
  bool occupancy_histogram_initialized=false; ///< Has occupancy histogram been initialized?

  /// Mutations from parent
  std::unordered_map<std::string, int> mutations;

  // ==== Internal genome analysis computations ====
  /// Calculates the histogram of number of genes overlapping at each bit
  /// for the current organism. The histogram gives the number of netural
  /// sites (bin 0), the number of bits with only one gene overlap (bin 1),
  /// the average amount of overlap at each bit (mean of dist.), and finally
  /// it gives the amount of multi-overlap (all bins > 1).
  /// Takes the size of the genome as an argument
  void HistogramCalc();
  void NeighborCalc();

public:
  AagosOrg(size_t _num_bits=64, size_t _num_genes=64, size_t _gene_size=8)
    : genome(_num_bits, _num_genes, _gene_size),
      phenotype(_num_genes),
      gene_neighbors(_num_genes)
  {
    emp_assert(genome.bits.size() > 0, genome.bits.size());
    emp_assert(genome.num_genes > 0, genome.num_genes);
    emp_assert(genome.gene_size > 0, genome.gene_size);
    emp_assert(!occupancy_histogram_initialized, occupancy_histogram_initialized);
  }

  AagosOrg(const Genome & g)
    : genome(g), phenotype(g.num_genes), gene_neighbors(g.num_genes)
  {
    emp_assert(genome.bits.size() > 0, genome.bits.size());
    emp_assert(genome.num_genes > 0, genome.num_genes);
    emp_assert(genome.gene_size > 0, genome.gene_size);
    emp_assert(!occupancy_histogram_initialized, occupancy_histogram_initialized);
  }

  AagosOrg(const AagosOrg &) = default;
  AagosOrg(AagosOrg &&) = default;
  ~AagosOrg() { ; }

  size_t GetNumBits() const { return genome.bits.size(); }
  size_t GetNumGenes() const { return genome.num_genes; }
  size_t GetGeneSize() const { return genome.gene_size; }

  emp::BitVector & GetBits() { return genome.bits; }
  const emp::BitVector & GetBits() const { return genome.bits; }
  const emp::vector<size_t> & GetGeneStarts() const { return genome.gene_starts; }

  Genome & GetGenome() { return genome; }
  const Genome & GetGenome() const { return genome; }

  void RandomizeGenome(emp::Random & random) {
    genome.Randomize(random);
  }

  // Phenotype information
  Phenotype & GetPhenotype() { return phenotype; }
  const Phenotype & GetPhenotype() const { return phenotype; }

  // Mutation information
  std::unordered_map<std::string, int> & GetMutations() { return mutations; }
  const std::unordered_map<std::string, int> & GetMutations() const { return mutations; }

  void ResetMutations() {
    for (auto & pair : mutations) {
      pair.second = 0;
    }
  }

  void ResetHistogram() {
    occupancy_histogram.Reset();
    occupancy_histogram_initialized = false;
  }

  /// Get number of neighbors for each gene. If it has not been calculated, calculate it.
  const emp::vector<size_t> & GetGeneNeighbors() {
    // If the histogram has yet to be setup, compute histogram.
    if (!occupancy_histogram_initialized) {
      StatsCalc();
    }
    return gene_neighbors;
  }

  /// Get gene occupancy histogram
  const histogram_t & GetGeneOccupancyHistogram() {
    // if the histogram has yet to be setup, compute histogram.
    if (!occupancy_histogram_initialized) {
      StatsCalc();
    }
    return occupancy_histogram;
  }

  /// Calculates histogram and gene neighbors for the current organism
  /// only called when a snapshot or statistics need to be taken for a pop
  /// b/c GetHistogram and GetNeighbors only called when snapshot and stats calc
  void StatsCalc() {
    // set sentinel
    occupancy_histogram_initialized = true;
    // set up histogram
    HistogramCalc();
    NeighborCalc();
  }

  /// Print function for aagos organism
  void Print(std::ostream & os = std::cout) const {
    os << "(" << emp::to_string(genome.gene_starts) << ")";
    os << "[" << genome.bits << "]";
  }

};

// todo - write test!
void AagosOrg::HistogramCalc() {
  // histogram bins ranges from 0 (no overlap) to num_genes, b/c worst case all
  // genes overlap the same bit. Num bins is then num_genes + 1 b/c need a
  // bin for no overlap.
  const size_t num_bits = GetNumBits();
  const size_t num_genes = GetNumGenes();
  const size_t gene_size = GetGeneSize();
  const size_t num_bins = num_genes + 1;
  const auto & gene_starts = genome.gene_starts;
  // Configure histogram bins.
  occupancy_histogram.SetupBins(0, (int)num_bins, num_bins);
  // adds the number of genes associated with each bit to histogram data node
  // first loops through each bit
  for (size_t i = 0; i < num_bits; ++i) {
    // counts number of genes that overlap at given bit
    size_t overlap = 0;
    // max range that a gene can start in order to overlap current bit
    const int diff = (int)i - ((int)gene_size - 1);  // off by one error
    // loops through each gene to see if it overlaps current bit
    for (size_t j = 0; j < num_genes; ++j) {
      // if current gene starts within diff of the current bit
      // and also if the current gene starts before the current bit
      // then the curr gene must overlap the curr bit
      const bool overlap_pos = (diff <= (int)gene_starts[j] && (int)gene_starts[j] <= (int)i)
                                || (diff <= ((int)gene_starts[j] - (int)num_bits));
      overlap += (size_t)overlap_pos;

    }
    // once total overlap for given bit is calc, add to hist.
    occupancy_histogram.Add((int)overlap);
  }
}

// TODO - write test
void AagosOrg::NeighborCalc() {
  std::fill(gene_neighbors.begin(), gene_neighbors.end(), 0); // Reset gene neighbor counts.
  const auto & gene_starts = genome.gene_starts;
  const size_t num_bits = GetNumBits();
  const size_t num_genes = GetNumGenes();
  const size_t gene_size = GetGeneSize();
  emp_assert(gene_neighbors.size() == num_genes);
  emp_assert(gene_starts.size() == num_genes);
  for (size_t i = 0; i < num_genes; ++i) {
    for (size_t j = i+1; j < num_genes; ++j) {
      // if the current gene starts w/in gene_size on either side of gene in question, must overlap
      // the second check catches genes that overlap only by comparing the ends of both genes to each other modded
      const bool neighbors = abs((int)gene_starts[i] - (int)gene_starts[j]) < gene_size ||
                              abs(((int)gene_starts[i] + (int)gene_size) % (int)num_bits - ((int)gene_starts[j] + (int)gene_size) % (int)num_bits) < (int)gene_size;
      gene_neighbors[i] += (size_t)neighbors;
      gene_neighbors[j] += (size_t)neighbors;
    }
  }
}

#endif
