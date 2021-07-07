#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include "assembly/assembly_utility.h"
#include "assembly/local_assembler.h"
#include "basic/bit_operation.h"
#include "basic/histgram.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "graph/scaffold_graph.h"
#include "misc/hash_aligner.h"
#include "misc/log.h"
#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/read_library.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"
#include "sequence/short_sequence.h"

using namespace std;

string directory = "out";

string kmer_file() { return directory + "/kmer"; }
string contig_file(int kmer_size) { return directory + FormatString("/contig-%d.fa", kmer_size); }
string graph_file(int kmer_size) { return directory + FormatString("/graph-%d.fa", kmer_size); }
int min_count = 2;
int prefix_length = 3;

void Assemble(HashGraph &hash_graph) {
    cout << "kmers " << hash_graph.num_vertices() << " " << hash_graph.num_edges() << endl;

    int kmer_size = hash_graph.kmer_size();

    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);
    hash_graph.clear();
    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(hash_graph);
    }

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    deque<Sequence> multi_contigs;
    deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    PrintN50(multi_contigs);
    WriteContig(contig_file(kmer_size), multi_contigs, multi_contig_infos, FormatString("contig-%d", kmer_size));
}

void BuildHashGraph(AssemblyInfo& assembly_info, int kmer_size) {
  cout << "building kmers " << endl;
  BuildKmerFile(assembly_info, kmer_size, min_count, prefix_length, kmer_file());

  HashGraph hash_graph(kmer_size);
  cout << "loading kmers into hashgraph" << endl;
  ReadKmerFile(kmer_file(), hash_graph);

  cout << "refresh edges" << endl;
  hash_graph.RefreshEdges();
  cout << "insert internal kmers" << endl;
  InsertInternalKmers(assembly_info, hash_graph, min_count);

  Assemble(hash_graph);
}

int main(int argc, char *argv[]) {
  AssemblyInfo assembly_info;
  ReadInput(argv[1], "", assembly_info);
  // 4.23user 0.48system 0:04.74elapsed 99%CPU (0avgtext+0avgdata 242100maxresident)k
  // 0inputs+0outputs (2major+60652minor)pagefaults 0swaps
  // $ ls -lh test.fa S2.fq.gz
  //   -rw-r--r--  1 fredyr  staff   151M Jun  9 08:35 S2.fq.gz
  //   -rw-r--r--  1 fredyr  staff   515M Jun  9 08:37 test.fa

  BuildHashGraph(assembly_info, atoi(argv[2]));

  // kmer=31
  // building kmers
  // loading kmers into hashgraph
  // refresh edges
  // insert internal kmers
  // kmers 5517097 5518100
  // contigs: 8615 n50: 4480 max: 20153 mean: 670 total length: 5775547 n80: 1789
  // 207.77user 4.25system 1:05.84elapsed 322%CPU (0avgtext+0avgdata 588864maxresident)k
  // 0inputs+0outputs (79major+158421minor)pagefaults 0swaps
  return 0;

}

/*
  kmer=31
  using ecoli samples below

  building kmers
  loading kmers into hashgraph
  refresh edges
  insert internal kmers
  kmers 8536897 8580057
  contigs: 397659 n50: 55 max: 1986 mean: 51 total length: 20466667 n80: 33
  378.69user 5.81system 2:00.36elapsed 319%CPU (0avgtext+0avgdata 824396maxresident)k
  0inputs+0outputs (45major+238546minor)pagefaults 0swaps


N=1 microbenchmarking of paired-end E. coli sample SRR1770413, kmer size 31
| Software                | Time (real) |
|-------------------------+-------------|
| Kmer histogram only     |             |
|                         |             |
| jellyfish (count+histo) | 1m.21.70s   |
| lh3/kmer-cnt/kc-c4      | 21.67s      |
| lh3/kmer-cnt/yak-count  | 27.47s      |
|                         |             |
| De Brujin graph         |             |
|                         |             |
| mccortex                | 5m58.980s   |
| fq2fa+IDBA_UD           | 5m50.73s    |
| GATB/bcalm              | 2m34.07s    |
|                         |             |
| DBA (wo fa conversion)  | 2m00.36s    |

The E. coli sample files were 113M and 127M respectively.

Notable is that Jellyfish 81.70/21.67 = 3.77x slower than kc-c4, and even though
bcalm is twice as fast as mccortex and IDBA, it is still 154.07/21.67 = 7.1x
slower than kc-c4 or 154.07/27.47 = 5.6x slower than yak-count.

The amount of work for building a DBG is ofc larger than just kmer counting, but
I think that would reasonably be around 2 to 3x at the most.

I think once assembly pass for this sample could be done under 90secs!

*/
