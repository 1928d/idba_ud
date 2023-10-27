/**
 * @file scaffold.cpp
 * @brief Build scaffolds on a set of contigs with multiple levels of insert libraries.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.6
 * @date 2011-08-06
 */

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

struct IDBAOption {
    string directory;
    string read_file;
    string long_read_file;
    int mink;
    int maxk;
    int step;
    int inner_mink;
    int inner_step;
    int prefix_length;
    int min_count;
    int min_support;
    int min_contig;
    double similar;
    int max_mismatch;
    int seed_kmer_size;
    int num_threads;
    int min_pairs;
    bool is_no_local;
    bool is_use_all;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;
    string input_contig;

    IDBAOption() {
        directory = "out";
        mink = 20;
        maxk = 100;
        step = 20;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 5;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        is_no_local = false;
        is_use_all = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    string log_file() { return directory + "/log"; }
    string kmer_file() { return directory + "/kmer"; }
    string align_file(int kmer_size) { return directory + FormatString("/align-%d", kmer_size); }
    string graph_file(int kmer_size) { return directory + FormatString("/graph-%d.fa", kmer_size); }
    string contig_file(int kmer_size) { return directory + FormatString("/contig-%d.fa", kmer_size); }
    string contig_info_file(int kmer_size) { return directory + FormatString("/contig-info-%d.fa", kmer_size); }
    string local_contig_file(int kmer_size) { return directory + FormatString("/local-contig-%d.fa", kmer_size); }
    string contig_file() { return directory + "/contig.fa"; }

    // string scaffold_file() { return directory + "/scaffold.fa"; }
    string scaffold_file(int level = 0) {
        return directory + (level == 0 ? "/scaffold.fa" : FormatString("/scaffold-level-%d.fa", level + 1));
    }
};

AssemblyInfo assembly_info;
IDBAOption option;
int median = 0;
int sd = 0;
int read_length = 0;

void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file);

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file);

int main(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "min size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");
    desc.AddOption("maxk", "", option.maxk, FormatString("maximum k value"));

    try {
        desc.Parse(argc, argv);

        if (argc < 2)
            throw logic_error("not enough parameters");

    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "scaffold - build scaffolds from contigs and multiple paired-end libraries." << endl;
        cerr << "Usage: scaffold -o output_dir reads-lib-1.fa [reads-lib-2.fa] [...]" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    MakeDir(option.directory);

    if (option.num_threads == 0)
        option.num_threads = omp_get_max_threads();
    else
        omp_set_num_threads(option.num_threads);
    cout << "number of threads " << option.num_threads << endl;

    // option.read_file = argv[1];

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);
    ScaffoldGraph scaffold_graph(option.maxk, contigs);

    deque<string> read_files;
    // read_files.push_back(option.read_file);
    for (int i = 1; i < argc; ++i)
        read_files.push_back(argv[i]);

    for (int level = 0; level < (int)read_files.size(); ++level)
        AddPairs(level, scaffold_graph, read_files[level], option.align_file(option.maxk) + FormatString("-%d", level));

    for (int level = 0; level < (int)read_files.size(); ++level) {
        scaffold_graph.BuildEdges();
        scaffold_graph.FilterEdges(option.min_pairs, scaffold_graph.sd(level) * 4);
        scaffold_graph.ParseEdges();
        cout << "edgs " << scaffold_graph.num_edges(level) << endl;
        scaffold_graph.print_edges(level);
        scaffold_graph.RemoveTransitiveConnections(level);

        deque<ContigGraphPath> paths;
        scaffold_graph.Assemble(level, paths);

        deque<Sequence> contigs;
        scaffold_graph.Assemble(level, contigs);
        PrintN50(contigs);
        WriteSequence(option.scaffold_file(level), contigs, "scaffold");
        scaffold_graph.Initialize(paths);
    }

    fflush(stdout);
    return 0;
}

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file) {
    ShortReadLibrary short_read_library;
    ReadLibrary(read_file, short_read_library);
    cout << "reads " << short_read_library.size() << endl;
    AlignReads(option.contig_file(option.maxk), short_read_library, align_file);
    double mean = 0;
    double sd = 0;
    EstimateDistance(align_file, mean, sd);
    ::median = mean;
    ::sd = sd;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);
    cout << "contig" << option.contig_file(option.maxk) << "#" << contigs.size() << endl;
    deque<ContigInfo> contig_infos(contigs.size());
    vector<int> num_aligned_reads(contigs.size(), 0);
    vector<double> coverage(contigs.size());

    deque<ShortSequence> &reads = short_read_library.reads();

    cout << "align file " << align_file << endl;
    FILE *falign = OpenFile(align_file, "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size) {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0) {
#pragma omp atomic
                ++num_aligned_reads[record.ref_id];
            }
        }
    }
    fclose(falign);

    double sum_coverage = 0;
    double sum_length = 0;
#pragma omp parallel for reduction(+ : sum_coverage, sum_length)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i) {
        if ((int)contigs[i].size() > option.min_contig) {
            sum_coverage += num_aligned_reads[i];
            sum_length += contigs[i].size() - read_length + 1;
            coverage[i] = 1.0 * num_aligned_reads[i] / (contigs[i].size() - reads[0].size() + 1);
            contig_infos[i].set_kmer_count(num_aligned_reads[i]);
        }
    }
    double mean_coverage = sum_coverage / sum_length;
    cout << "expected coverage " << mean_coverage << endl;

    int num_connections = 0;
    falign = OpenFile(align_file, "rb");
    for (unsigned i = 0; i < reads.size(); i += 2) {
        deque<HashAlignerRecord> records1;
        deque<HashAlignerRecord> records2;
        ReadHashAlignerRecords(falign, records1);
        ReadHashAlignerRecords(falign, records2);

        for (unsigned j = 0; j < records1.size(); ++j) {
            for (unsigned k = 0; k < records2.size(); ++k) {
                HashAlignerRecord &r1 = records1[j];
                HashAlignerRecord &r2 = records2[k];
                r2.ReverseComplement();

                if (r1.ref_length > option.min_contig && r2.ref_length > option.min_contig &&
                    r1.ref_from - r1.query_from > r1.ref_length - median - 3 * sd &&
                    r2.ref_to + r2.query_length - r2.query_to < median + 3 * sd && r1.ref_id != r2.ref_id) {
                    int d = median - (r1.ref_length - (r1.ref_from - r1.query_from)) -
                            (r2.ref_to + r2.query_length - r2.query_to);
                    scaffold_graph.AddPair(level, (r1.ref_id * 2 + r1.is_reverse), (r2.ref_id * 2 + r2.is_reverse), d);
                    ++num_connections;
                }
            }
        }
    }

    scaffold_graph.set_library_info(level, read_length, mean_coverage, median, sd);
}

void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file) {
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    assembly_info.reads = library.reads();
    assembly_info.read_flags.resize(assembly_info.reads.size());

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);
    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    cout << "aligned " << num_aligned_reads << " reads" << endl;
}
