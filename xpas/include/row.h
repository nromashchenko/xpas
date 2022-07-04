#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <xcl/phylo_kmer.h>
#include <array>

namespace xpas
{
    using xcl::phylo_kmer;

    struct proba_pair
    {
        phylo_kmer::score_type score;
        phylo_kmer::key_type index;
    };

    using branch_type = phylo_kmer::branch_type;
    using row_type = std::array<proba_pair, xcl::seq_traits::alphabet_size>;
}

#endif //RAPPAS_CPP_ROW_H
