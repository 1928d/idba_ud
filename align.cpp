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

void Aligner(const string &contig_file, ShortReadLibrary &library, const string &align_file);
void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file);
void PAF(const HashAlignerRecord &r);

int main(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("read", "r", option.read_file, FormatString("fasta read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("long_read", "l", option.long_read_file,
                   FormatString("fasta long read file (>%d)", ShortSequence::max_size()));
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "min size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");

    try {
        desc.Parse(argc, argv);
        if (argc < 2)
            throw logic_error("not enough parameters");

    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "align - align reads against contigs." << endl;
        cerr << "Usage: align -r reads.fa -o output_dir contig.fa" << endl;
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

    string ctgs = argv[1];
    string align_file = option.align_file(option.seed_kmer_size);

    ShortReadLibrary short_read_library;
    ReadLibrary(option.read_file, short_read_library);
    cout << "reads " << short_read_library.size() << endl;
    Aligner(ctgs, short_read_library, align_file);

    deque<ShortSequence> &reads = short_read_library.reads();
    int num_aligned_reads = 0;

    FILE *falign = OpenFile(align_file, "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size) {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);
        ReadHashAlignerRecordBlock(falign, all_records);
        for (int i = 0; i < size; ++i) {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0) {
                PAF(record);
            }
        }
    }
    fclose(falign);

    cout << "aligned " << num_aligned_reads << " reads" << endl;
    fflush(stdout);
    return 0;
}

void PAF(const HashAlignerRecord &r) {
    cout << r.query_length << "\t" << r.query_from << "\t" << r.query_to << "\t" << (r.is_reverse ? '-' : '+') << "\t"
         << r.ref_id << "\t" << r.ref_length << "\t" << r.ref_from << "\t" << r.ref_to << "\t" << r.match_length
         << endl;
}

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file) {
    // int num_connections = 0;
    // falign = OpenFile(align_file, "rb");
    // for (unsigned i = 0; i < reads.size(); i += 2) {
    //     deque<HashAlignerRecord> records1;
    //     deque<HashAlignerRecord> records2;
    //     ReadHashAlignerRecords(falign, records1);
    //     ReadHashAlignerRecords(falign, records2);

    //     for (unsigned j = 0; j < records1.size(); ++j) {
    //         for (unsigned k = 0; k < records2.size(); ++k) {
    //             HashAlignerRecord &r1 = records1[j];
    //             HashAlignerRecord &r2 = records2[k];
    //             r2.ReverseComplement();

    //             if (r1.ref_length > option.min_contig && r2.ref_length > option.min_contig &&
    //                 r1.ref_from - r1.query_from > r1.ref_length - median - 3 * sd &&
    //                 r2.ref_to + r2.query_length - r2.query_to < median + 3 * sd && r1.ref_id != r2.ref_id) {
    //                 int d = median - (r1.ref_length - (r1.ref_from - r1.query_from)) -
    //                         (r2.ref_to + r2.query_length - r2.query_to);
    //                 scaffold_graph.AddPair(level, (r1.ref_id * 2 + r1.is_reverse), (r2.ref_id * 2 + r2.is_reverse),
    //                 d);
    //                 ++num_connections;
    //             }
    //         }
    //     }
    // }
}

int64_t AlignReadsEx(AssemblyInfo &assembly_info, HashAligner &hash_aligner, double similar,
                     const std::string &align_file, bool is_all) {
    deque<ShortSequence> &reads = assembly_info.reads;
    vector<bool> &read_flags = assembly_info.read_flags;

    FILE *falign = OpenFile(align_file, "wb");

    int64_t num_aligned_reads = 0;
    int buffer_size = (1 << 20) * omp_get_max_threads();
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size) {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<deque<HashAlignerRecord>> buffer_records(omp_get_max_threads());
        vector<HashAlignerRecord> all_records(size);

#pragma omp parallel for schedule(static, 1) reduction(+ : num_aligned_reads)
        for (int64_t i = 0; i < size; ++i) {
            all_records[i].match_length = 0;
            if (read_flags[offset + i] || is_all) {
                deque<HashAlignerRecord> &records = buffer_records[omp_get_thread_num()];
                Sequence seq(reads[offset + i]);
                hash_aligner.AlignRead(seq, records, seq.size() * similar, 1);

                if (records.size() == 1) {
                    ++num_aligned_reads;
                    all_records[i] = records[0];
                }
            }
        }
        WriteHashAlignerRecordBlock(falign, all_records);
    }

    fclose(falign);

    return num_aligned_reads;
}

void Aligner(const string &contig_file, ShortReadLibrary &library, const string &align_file) {
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    assembly_info.reads = library.reads();
    assembly_info.read_flags.resize(assembly_info.reads.size());

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);
    int64_t num_aligned_reads = AlignReadsEx(assembly_info, hash_aligner, option.similar, align_file, true);
    cout << "aligned " << num_aligned_reads << " reads" << endl;
}
