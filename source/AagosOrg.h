#ifndef AAGOS_ORG_H
#define AAGOS_ORG_H

#include "tools/BitVector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/string_utils.h"
#include "data/DataNode.h"

class AagosOrg
{
  friend class AagosWorld;

private:
  // genome of organism - bitstring
  emp::BitVector bits;
  // starting locations of all genes
  emp::vector<size_t> gene_starts;
  // size of each gene in genome
  int gene_size;
  // number of genes in genome
  int num_genes;
  size_t num_bins;
  // number of neighbors each gene in genome has
  // neighbor is defined as a gene that overlaps the current gene
  // at at least one bit
  emp::vector<int> gene_neighbors;
  // Histogram object that stores num overlapped genes at each bit in genome
  emp::DataNode<int, emp::data::Histogram, emp::data::Stats> histogram;
  // bool flag to check if histogram has been initialized yet
  bool initialized;

public:
  AagosOrg(size_t num_bits = 64, size_t num_genes = 64, size_t in_gene_size = 8)
      : bits(num_bits)
      , gene_starts(num_genes, 0)
      , gene_size(in_gene_size)
      , num_genes(num_genes)
      , num_bins(num_genes + 1)
      , gene_neighbors(num_genes)
      , initialized(false)
  {
    emp_assert(num_bits > 0, num_bits);
    emp_assert(num_genes > 0, num_genes);
    emp_assert(gene_size > 0, gene_size);
    emp_assert(!initialized, initialized);
  }

  // default copy constructor
  AagosOrg(const AagosOrg &) = default;
  AagosOrg(AagosOrg &&) = default;
  ~AagosOrg() { ; }

  // default comparision
  AagosOrg &operator=(const AagosOrg &) = default;
  AagosOrg &operator=(AagosOrg &&) = default;

  // getter function for size of organism genome
  size_t GetNumBits() const { return bits.size(); }
  // getter function for number of genes
  size_t GetNumGenes() const { return gene_starts.size(); }
  // getter for organism genome
  const emp::BitVector &GetBits() const { return bits; }
  // getter for gene start locations
  const emp::vector<size_t> &GetGeneStarts() const { return gene_starts; }
  // getter for number of bins in histogram
  const size_t GetNumBins() { return num_bins; }

  // randomizes genome and gene starts
  void Randomize(emp::Random &random)
  {
    emp::RandomizeBitVector(bits, random);
    emp::RandomizeVector<size_t>(gene_starts, random, 0, bits.size());
  }

  // print override for aagos organism
  void Print(std::ostream &is = std::cout) const
  {
    is << "Bits: " << bits << '\n';
    is << "Gene Starts: " << emp::to_string(gene_starts) << std::endl;
  }

  // getter function for gene neighbors
  emp::vector<int> &GetGeneNeighbors()
  {
    // if the histogram hasn't been set up, calculate
    if (!initialized)
      StatsCalc();
    return gene_neighbors;
  }

  // getter function for gene overlap histogram
  emp::DataNode<int, emp::data::Histogram, emp::data::Stats> &GetHistogram()
  {
    // if the histogram hasn't been set up, calculate
    if (!initialized)
      StatsCalc();
    return histogram;
  }

  // calculates histogram and gene neighbors for the current organism
  // only called when a snapshot or statistics need to be taken for a pop
  // b/c GetHistogram and GetNeighbors only called when snapshot and stats calc
  void StatsCalc()
  {
    // set sentinel
    initialized = true;
    int num_bits = GetNumBits();
    // set up histogram
    HistogramCalc(num_bits);
    NeighborCalc(num_bits);
  }

  // Calculates the histogram of number of genes overlapping at each bit
  // for the current organism. The histogram gives the number of netural
  // sites (bin 0), the number of bits with only one gene overlap (bin 1),
  // the average amount of overlap at each bit (mean of dist.), and finally
  // it gives the amount of multi-overlap (all bins > 1).
  // Taks the size of the genome as an argument
  void HistogramCalc(int num_bits)
  {

    // histogram bins ranges from 0 (no overlap) to num_genes, b/c worst case all
    // genes overlap the same bit. Num bins is then num_genes + 1 b/c need a
    // bin for no overlap.
    histogram.SetupBins(0, num_bins, num_bins);

    // adds the number of genes associated with each bit to histogram data node
    // first loops through each bit
    for (int i = 0; i < (int) bits.size(); i++)
    {
      // counts number of genes that overlap at given bit
      int overlap = 0;
      // max range that a gene can start in order to overlap current bit
      int diff = i - (gene_size - 1);// off by one error

      // loops through each gene to see if it overlaps current bit
      for (size_t j = 0; j < num_genes; j++)
      {
        // if current gene starts within diff of the current bit
        // and also if the current gene starts before the current bit
        // then the curr gene must overlap the curr bit
        if ((diff <= (int) gene_starts[j] && (int) gene_starts[j] <= i)
            // else if curr gene loops around to front of genome,
            // then if it loops around enough to overlap current bit
            // then diff will be more negative than the
            //  gene start - the number of bits in genome
            || ((diff <= ((int) gene_starts[j] - num_bits)))) 
        {
          overlap++;
        }
      }
      // once total overlap for given bit is calc, add to hist.
      histogram.Add(overlap);
    }
  }

  // calculates all neighboring genes for each gene in genome.
  // Neighbor genes are defined as genes that share at least 1
  // bit with the curr gene.
  // takes the size of the genome as its argument
  void NeighborCalc(int num_bits)
  {
    // loops through each gene to get its neighbor count
    for (size_t i = 0; i < num_genes; i++)
    {
      size_t count = 0;
      // loops through all remaining genes to see if neighbors
      for (size_t j = 0; j < num_genes; j++)
      {
        // clearly, the curr gene should not count towards its own neighbor count
        if (i != j)
        {
          // if the current gene starts w/in gene_size on either side of gene in question, must overlap
          if (abs((int)gene_starts[i] - (int)gene_starts[j]) < gene_size || // todo clean up these casts
              // this second check catches genes that overlap only by comparing the ends of both genes to each other modded
              abs(((int)gene_starts[i] + gene_size) % num_bits - ((int)gene_starts[j] + gene_size) % num_bits) < gene_size)
          {
            count++;
          }
        }
      }
      // save given gene's neighbor count
      gene_neighbors[i] = count;
    }
  }
};
#endif
