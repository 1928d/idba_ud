/**
 * @file hash_aligner.h
 * @brief HashAligner class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-12
 */

#ifndef __MISC_HASH_ALIGNER_H_

#define __MISC_HASH_ALIGNER_H_

#include "basic/kmer.h"
#include "container/hash_map.h"
#include "container/hash_set.h"
#include "container/managed_list.h"
#include "sequence/sequence.h"

#include <deque>
#include <list>
#include <set>

class Sequence;

/**
 * @brief It is a alignment record generated by HashAligner.
 */
struct HashAlignerRecord {
    int query_from;
    int query_to;
    int query_length;
    int ref_id;
    int ref_from;
    int ref_to;
    int ref_length;
    int match_length;
    bool is_reverse;

    bool operator<(const HashAlignerRecord &record) const {
        if (ref_id != record.ref_id)
            return ref_id < record.ref_id;
        if (is_reverse != record.is_reverse)
            return is_reverse < record.is_reverse;
        else if (ref_from - query_from != record.ref_from - record.query_from)
            return ref_from - query_from < record.ref_from - record.query_from;
        else
            return ref_from < record.ref_from;
    }

    bool Merge(const HashAlignerRecord &record) {
        if (ref_id != record.ref_id || is_reverse != record.is_reverse)
            return false;

        if (ref_from - query_from != record.ref_from - record.query_from)
            return false;

        if (ref_to < record.ref_from - 5)
            return false;

        //        if (record.query_from <= query_to)
        //            match_length += record.query_to - query_to;
        //        else
        //            match_length += record.match_length;

        query_to = record.query_to;
        ref_to = record.ref_to;

        return true;
    }

    const HashAlignerRecord &ReverseComplement() {
        int to = ref_length - ref_from;
        int from = ref_length - ref_to;
        ref_from = from;
        ref_to = to;

        to = query_length - query_from;
        from = query_length - query_to;
        query_from = from;
        query_to = to;

        is_reverse = !is_reverse;

        return *this;
    }
};

/**
 * @brief It is an aligner used to align sequence to reference by simple
 * hashing.
 */
class HashAligner {
    typedef uint64_t Position;

  public:
    typedef ManagedList<Position>::node_pool_type pool_type;
    typedef ManagedList<Position> position_list_type;
    typedef position_list_type::iterator position_iterator;

    explicit HashAligner(uint32_t kmer_size = 0, uint32_t min_length = 0, uint32_t step = 1);
    ~HashAligner() {}

    void Initialize(const std::deque<Sequence> &sequences);
    int AlignRead(const Sequence &seq, std::deque<HashAlignerRecord> &records, int min_match,
                  int max_records = (1 << 30));
    int AlignReadLocal(const Sequence &seq, std::deque<HashAlignerRecord> &records, int min_match, int max_mismatch,
                       int max_records = (1 << 30));
    int AlignSequence(const Sequence &seq, std::deque<HashAlignerRecord> &records, int min_match, double similar,
                      int max_records = (1 << 30));

    uint32_t kmer_size() const { return kmer_size_; }
    void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

    uint32_t min_length() const { return min_length_; }
    void set_min_lenghth(int min_length) { min_length_ = min_length; }

    uint32_t step() const { return step_; }
    void set_step(uint32_t step) { step_ = step; }

    void clear() {
        hash_map_.clear();
        sequences_.clear();
        // sequence_sizes_.clear();
        pool_.clear();
    }

  private:
    HashAligner(const HashAligner &);
    const HashAligner &operator=(const HashAligner &);

    void InsertSequence(const Sequence &seq, uint64_t id);
    void ExtendRecord(const Sequence &seq, HashAlignerRecord &record);
    void ExtendRecord(const Sequence &seq, HashAlignerRecord &record, int max_mismatch);
    void Match(const Sequence &seq, HashAlignerRecord &record);

    void Match(const std::vector<uint64_t> &words, HashAlignerRecord &record, int max_mismatch);
    void Convert(const Sequence &seq, std::vector<uint64_t> &words);

    uint64_t GetWord(const std::vector<uint64_t> &words, int from, int to) {
        if (from == to)
            return 0;

        uint64_t word = 0;
        if ((from >> 5) == ((to - 1) >> 5))
            word = (words[from >> 5] >> ((from & 31) << 1));
        else
            word = ((words[from >> 5] >> ((from & 31) << 1)) | (words[(from >> 5) + 1] << ((32 - (from & 31)) << 1)));

        if (to - from < 32)
            word &= (1ULL << ((to - from) << 1)) - 1;

        return word;
    }

    pool_type pool_;
    HashMap<Kmer, position_list_type> hash_map_;
    std::deque<Sequence> sequences_;
    std::deque<Sequence> reverse_sequences_;
    std::deque<std::vector<uint64_t>> sequence_words;
    std::deque<std::vector<uint64_t>> reverse_sequence_words;
    uint32_t kmer_size_;
    uint32_t min_length_;
    uint32_t step_;

    std::deque<std::vector<HashAlignerRecord>> buffer_records;
    // std::deque<HashSet<uint64_t> > buffer_tables;
    std::deque<std::set<uint64_t>> buffer_tables;
};

#endif
