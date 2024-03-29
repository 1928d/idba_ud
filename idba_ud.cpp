/**
 * @file idba_ud.cpp
 * @brief An iterative de Bruijn graph assembler for sequencing data with highly
 * uneven depth.
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
#include "extras.hpp"
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
    string read_file_1, read_file_2;
    string long_read_file;
    deque<string> extra_read_files;
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
    int max_gap;
    bool is_no_bubble;
    bool is_no_local;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;
    string reference;

    IDBAOption() {
        extra_read_files.resize(4);
        directory = "out";
        mink = 20;
        maxk = 100;
        step = 20;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 3;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        max_gap = 50;
        is_no_bubble = false;
        is_no_local = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    string log_file() { return directory + "/log"; }

    string kmer_file() { return directory + "/kmer"; }

    string align_file(int kmer_size) { return directory + FormatString("/align-%d", kmer_size); }

    string graph_file(int kmer_size) { return directory + FormatString("/graph-%d.fa", kmer_size); }
    string gfa_file(int kmer_size) { return directory + FormatString("/graph-%d.gfa", kmer_size); }

    string contig_file(int kmer_size) { return directory + FormatString("/contig-%d.fa", kmer_size); }

    string local_contig_file(int kmer_size) { return directory + FormatString("/local-contig-%d.fa", kmer_size); }

    string contig_file() { return directory + "/contig.fa"; }

    string scaffold_file(int level = 0) {
        return directory + (level == 0 ? "/scaffold.fa" : FormatString("/scaffold-level-%d.fa", level + 1));
    }

    string ref_contig_file() { return directory + "/ref_contig.fa"; }
};

AssemblyInfo assembly_info;
IDBAOption option;
double median = 0;
double sd = 0;
int read_length = 0;

void BuildHashGraph(int kmer_size);
void Assemble(HashGraph &hash_graph);
void AlignReads(const string &contig_file, const string &align_file);
void CorrectReads(int kmer_size);
void LocalAssembly(int kmer_size, int new_kmer_size);
void Iterate(int kmer_size, int new_kmer_size);
// void Scaffold(int kmer_size, int min_contig);
// void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file);
void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file);
void ReadInputFastq(const std::string &read_file_1, const std::string &read_file_2, const std::string &long_read_file,
                    AssemblyInfo &assembly_info);

int main(int argc, char *argv[]) {
    OptionsDescription desc;

    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("read_1", "r", option.read_file_1,
                   FormatString("fasta/q (gz) read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("read_2", "s", option.read_file_2,
                   FormatString("fasta/q (gz) read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("read_level_2", "", option.extra_read_files[0], "paired-end reads fasta for second level scaffolds");
    desc.AddOption("read_level_3", "", option.extra_read_files[1], "paired-end reads fasta for third level scaffolds");
    desc.AddOption("read_level_4", "", option.extra_read_files[2], "paired-end reads fasta for fourth level scaffolds");
    desc.AddOption("read_level_5", "", option.extra_read_files[3], "paired-end reads fasta for fifth level scaffolds");
    desc.AddOption("long_read", "l", option.long_read_file,
                   FormatString("fasta long read file (>%d)", ShortSequence::max_size()));
    // desc.AddOption("reference", "", option.reference, "reference genome");
    desc.AddOption("mink", "", option.mink, FormatString("minimum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("maxk", "", option.maxk, FormatString("maximum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("step", "", option.step, "increment of k-mer of each iteration");
    desc.AddOption("inner_mink", "", option.inner_mink, "inner minimum k value");
    desc.AddOption("inner_step", "", option.inner_step, "inner increment of k-mer");
    desc.AddOption("prefix", "", option.prefix_length, "prefix length used to build sub k-mer table");
    desc.AddOption("min_count", "", option.min_count,
                   "minimum multiplicity for filtering k-mer when building the graph");
    desc.AddOption("min_support", "", option.min_support, "minimum supoort in each iteration");
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "minimum size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("max_mismatch", "", option.max_mismatch, "max mismatch of error correction");
    desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");
    // desc.AddOption("max_gap", "", option.max_gap, "maximum gap in
    // reference");
    desc.AddOption("no_bubble", "", option.is_no_bubble, "do not merge bubble");
    desc.AddOption("no_local", "", option.is_no_local, "do not use local assembly");
    desc.AddOption("no_coverage", "", option.is_no_coverage, "do not iterate on coverage");
    desc.AddOption("no_correct", "", option.is_no_correct, "do not do correction");
    desc.AddOption("pre_correction", "", option.is_pre_correction, "perform pre-correction before assembly");

    try {
        desc.Parse(argc, argv);

        if (option.read_file_1 == "" && option.long_read_file == "")
            throw logic_error("not enough parameters");

        if (option.maxk < option.mink)
            throw invalid_argument("mink is larger than maxk");

        if (option.maxk > (int)Kmer::max_size())
            throw invalid_argument("maxk is too large");
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "IDBA-UD - Iterative de Bruijn Graph Assembler for sequencing "
                "data "
                "with highly uneven depth."
             << endl;
        cerr << "Usage: idba_ud -r read.fa -o output_dir" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    MakeDir(option.directory);

    LogThread log_thread(option.log_file());

    string begin_file = option.directory + "/begin";
    fclose(OpenFile(begin_file, "wb"));

    if (option.num_threads == 0)
        option.num_threads = omp_get_max_threads();
    else
        omp_set_num_threads(option.num_threads);
    log() << "number of threads " << option.num_threads << endl;

    timer loading_fastq;
    ReadInputFastq(option.read_file_1, option.read_file_2, option.long_read_file, assembly_info);
    // ReadInput(option.read_file_1, option.long_read_file, assembly_info);
    deque<Sequence> extra_reads;
    for (unsigned i = 0; i < option.extra_read_files.size(); ++i) {
        if (option.extra_read_files[i] != "") {
            deque<Sequence> reads;
            ReadSequence(option.extra_read_files[i], reads);
            extra_reads.insert(extra_reads.end(), reads.begin(), reads.end());
        }
    }
    log("perf") << "loading fastq " << loading_fastq.duration() << endl;
    log() << "reads " << assembly_info.reads.size() << endl;
    log() << "long reads " << assembly_info.long_reads.size() << endl;
    log() << "extra reads " << extra_reads.size() << endl;

    assembly_info.long_reads.insert(assembly_info.long_reads.end(), extra_reads.begin(), extra_reads.end());
    assembly_info.ClearStatus();

    read_length = assembly_info.read_length();
    log() << "read_length " << read_length << endl;

    if (option.is_pre_correction) {
        int kmer_size = (option.maxk + option.mink) / 2;
        log() << "kmer " << kmer_size << endl;
        BuildHashGraph(kmer_size);
        AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));
        CorrectReads(kmer_size);
        assembly_info.ClearStatus();
    }

    int old_kmer_size = 0;
    int kmer_size = option.mink;
    while (true) {
        log() << "kmer " << kmer_size << endl;

        if (kmer_size >= (option.mink + option.maxk) / 2 || kmer_size == option.maxk)
            assembly_info.ref_contigs.clear();

        if (kmer_size == option.mink)
            BuildHashGraph(kmer_size);
        else
            Iterate(old_kmer_size, kmer_size);

        if (kmer_size < option.maxk) {
            AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));
            CorrectReads(kmer_size);
            assembly_info.ClearStatus();

            old_kmer_size = kmer_size;
            kmer_size = min(option.maxk, kmer_size + option.step);
            LocalAssembly(old_kmer_size, kmer_size);

            if (old_kmer_size == option.maxk)
                break;
        } else
            break;
    }

    kmer_size = option.maxk;

    deque<Sequence> contigs;
    deque<string> names;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    FastaWriter writer(option.contig_file());
    for (unsigned i = 0; i < contigs.size(); ++i) {
        if ((int)contigs[i].size() >= option.min_contig)
            writer.Write(contigs[i], names[i]);
    }

    // NOTE: we're skipping the scaffolding entirely!
    //
    // Scaffold(option.maxk, option.min_contig);

    string end_file = option.directory + "/end";
    fclose(OpenFile(end_file, "wb"));

    fflush(stdout);

    return 0;
}

void BuildHashGraph(int kmer_size) {
    timer kmer_cnt;
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file());

    HashGraph hash_graph(kmer_size);
    ReadKmerFile(option.kmer_file(), hash_graph);

    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);

    if (option.reference != "") {
        deque<Sequence> ref_contigs;
        ReadSequence(option.ref_contig_file(), ref_contigs);
#pragma omp parallel for
        for (int64_t i = 0; i < (int64_t)ref_contigs.size(); ++i)
            hash_graph.InsertUncountKmers(ref_contigs[i]);
        hash_graph.RefreshEdges();
    }

    log("perf") << "kmer count " << kmer_cnt.duration() << endl;
    Assemble(hash_graph);
}

void Assemble(HashGraph &hash_graph) {
    log() << "kmers " << hash_graph.num_vertices() << " " << hash_graph.num_edges() << endl;
    timer assemb;

    int kmer_size = hash_graph.kmer_size();
    double min_cover = max(1, (kmer_size == option.mink ? option.min_count : option.min_support));

    Histgram<int> hist = hash_graph.coverage_histgram();
    double expected_coverage = hist.mean();

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

    contig_graph.RemoveDeadEnd(option.min_contig);

    if (!option.is_no_bubble) {
        int bubble = contig_graph.RemoveBubble();
        log() << "merge bubble " << bubble << endl;
        contig_graph.MergeSimilarPath();
    }

    if (!option.is_no_coverage)
        contig_graph.RemoveLocalLowCoverage(min_cover, option.min_contig, 0.1);

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(option.graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage) {
        double ratio = (kmer_size < option.maxk) ? 0.5 : 0.2;
        if (ratio < 2.0 / expected_coverage)
            ratio = 2.0 / expected_coverage;
        contig_graph.IterateLocalCoverage(option.min_contig, ratio, min_cover, 1e100, 1.1);
        contig_graph.MergeSimilarPath();
    }

    deque<Sequence> multi_contigs;
    deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    auto dead = contig_graph.CountDeadEnds();
    PrintN50(multi_contigs, dead);
    // WriteSequence(option.contig_file(kmer_size), multi_contigs);
    WriteContig(option.contig_file(kmer_size), multi_contigs, multi_contig_infos, FormatString("contig-%d", kmer_size));
    // WriteContigInfo(option.contig_info_file(kmer_size), multi_contig_infos);

    std::ofstream gfa_graph(option.gfa_file(kmer_size));
    contig_graph.PrintGFASegments(gfa_graph);
    contig_graph.PrintGFAEdges(gfa_graph);

    log("perf") << "assemble graph " << assemb.duration() << endl;
}

void AlignReads(const string &contig_file, const string &align_file) {
    timer align;
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    log("perf") << "align reads " << align.duration() << endl;
    log() << "aligned " << num_aligned_reads << " reads" << endl;
}

void CorrectReads(int kmer_size) {
    if (option.is_no_correct)
        return;

    timer corr;
    deque<Sequence> contigs;
    deque<string> names;
    deque<ContigInfo> contig_infos;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    CorrectReads(assembly_info, contigs, contig_infos, option.align_file(kmer_size), option.max_mismatch);
    // WriteSequence(option.contig_file(kmer_size), contigs);
    WriteContig(option.contig_file(kmer_size), contigs, contig_infos, FormatString("contig-%d", kmer_size));
    log("perf") << "correct reads " << corr.duration() << endl;
}

void LocalAssembly(int kmer_size, int new_kmer_size) {
    // if (median == 0)
    // EstimateDistance(kmer_size);
    timer local_asm;

    EstimateDistance(option.align_file(kmer_size), median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2 * median) {
        log() << "invalid insert distance" << endl;
        deque<Sequence> local_contigs;
        WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
        return;
    }

    deque<ShortSequence> &reads = assembly_info.reads;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    LocalAssembler local_assembler;
    local_assembler.Initialize(assembly_info, contigs);
    local_assembler.set_num_threads(option.num_threads);
    local_assembler.set_mink(option.inner_mink);
    local_assembler.set_maxk(new_kmer_size);
    local_assembler.set_step(option.inner_step);
    local_assembler.set_min_contig(option.min_contig);
    local_assembler.set_insert_distance(median, sd);

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size) {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
                local_assembler.AddReadByHashAlignerRecord(record, offset + i);
        }
    }
    fclose(falign);

    deque<Sequence> local_contigs;

    if (!option.is_no_local)
        local_assembler.Assemble(local_contigs);

    int num_seed_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i) {
        if ((int)contigs[i].size() > option.min_contig)
            ++num_seed_contigs;
    }

    log() << "seed contigs " << num_seed_contigs << " local contigs " << local_contigs.size() << endl;
    WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
    log("perf") << "local assembly " << local_asm.duration() << endl;
}

void Iterate(int kmer_size, int new_kmer_size) {
    timer iter;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    deque<Sequence> local_contigs;
    ReadSequence(option.local_contig_file(kmer_size), local_contigs);

    deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(kmer_size), multi_contigs);

    uint64_t sum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
        sum += contigs[i].size();
    HashGraph hash_graph(kmer_size);
    hash_graph.reserve(sum);

    deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    contigs.clear();
    local_contigs.clear();
    multi_contigs.clear();

    IterateHashGraph(assembly_info, new_kmer_size, option.min_support, hash_graph, old_contigs);
    kmer_size = new_kmer_size;
    old_contigs.clear();

    // if (kmer_size < option.maxk)
    // if (kmer_size < read_length)
    hash_graph.RefreshEdges();
    //    else
    //        hash_graph.AddAllEdges();

    log("perf") << "iterate " << iter.duration() << endl;
    Assemble(hash_graph);
}

void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file) {
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    assembly_info.reads.swap(library.reads());
    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    log() << "aligned " << num_aligned_reads << " reads" << endl;

    assembly_info.reads.swap(library.reads());
}

///
/// Allow read_file to be both FASTA, FASTQ, with and without GZ
///
#include "kseq.h"
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

uint64_t ReadSequenceEx(const string &read_r, const string &read_s, deque<ShortSequence> &sequences) {
    gzFile fpr = gzopen(read_r.c_str(), "r");
    kseq_t *sqr = kseq_init(fpr);
    gzFile fps = 0;
    kseq_t *sqs = 0;

    bool paired = read_s != "";
    if (paired) {
        fps = gzopen(read_s.c_str(), "r");
        sqs = kseq_init(fps);
    }
    sequences.clear();
    while (kseq_read(sqr) >= 0 && (paired ? (kseq_read(sqs) >= 0) : true)) {
        Sequence r(string(sqr->seq.s));
        r.TrimN();
        ShortSequence short_r(r);
        sequences.push_back(short_r);
        if (paired) {
            Sequence s(string(sqs->seq.s));
            s.TrimN();
            ShortSequence short_s(s);
            sequences.push_back(short_s);
        }
    }

    kseq_destroy(sqr);
    gzclose(fpr);
    if (paired) {
        kseq_destroy(sqs);
        gzclose(fps);
    }
    return sequences.size();
}

void ReadInputFastq(const std::string &read_file_1, const std::string &read_file_2, const std::string &long_read_file,
                    AssemblyInfo &assembly_info) {
    if (read_file_1 != "")
        ReadSequenceEx(read_file_1, read_file_2, assembly_info.reads);
    if (long_read_file != "")
        ReadSequence(long_read_file, assembly_info.long_reads);

    assembly_info.read_flags.resize(assembly_info.reads.size());
    assembly_info.long_read_flags.resize(assembly_info.long_reads.size());

    assembly_info.ClearStatus();
}
