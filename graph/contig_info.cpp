#include "graph/contig_info.h"

#include <fstream>
#include <istream>
#include <ostream>
#include <string>

#include "graph/bit_edges.h"

using namespace std;

istream &operator>>(istream &is, ContigInfo &contig_info) {
    is.read((char *)&contig_info.in_edges_, sizeof(BitEdges));
    is.read((char *)&contig_info.out_edges_, sizeof(BitEdges));
    is.read((char *)&contig_info.kmer_size_, sizeof(uint16_t));
    is.read((char *)&contig_info.kmer_count_, sizeof(uint32_t));

    int size = 0;
    if (!is.read((char *)&size, sizeof(int)))
        return is;

    contig_info.counts_.resize(size);
    for (int i = 0; i < size; ++i)
        is.read((char *)&contig_info.counts_[i], sizeof(SequenceCountUnitType));

    return is;
}

ostream &operator<<(ostream &os, const ContigInfo &contig_info) {
    os.write((char *)&contig_info.in_edges_, sizeof(BitEdges));
    os.write((char *)&contig_info.out_edges_, sizeof(BitEdges));
    os.write((char *)&contig_info.kmer_size_, sizeof(uint16_t));
    os.write((char *)&contig_info.kmer_count_, sizeof(uint32_t));

    int size = contig_info.counts_.size();
    os.write((char *)&size, sizeof(int));
    for (int i = 0; i < size; ++i)
        os.write((char *)&contig_info.counts_[i], sizeof(SequenceCountUnitType));

    return os;
}
