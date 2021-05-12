#include "xpas/seq.h"

using namespace xpas;

#ifdef SEQ_TYPE_DNA
const seq_traits_impl<dna>::char_type seq_traits_impl<dna>::code_to_key[] = {'A', 'C', 'G', 'T'};

const seq_traits_impl<dna>::key_type seq_traits_impl<dna>::pow_sigma[] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824 };

#elif SEQ_TYPE_AA

const seq_traits_impl<aa>::char_type seq_traits_impl<aa>::code_to_key[] = {
            'R', 'H', 'K',
            'D', 'E',
            'S', 'T', 'N', 'Q',
            'C', 'G', 'P',
            'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V'
        };

const seq_traits_impl<aa>::key_type seq_traits_impl<aa>::pow_sigma[] = { 1, 20, 400, 8000, 160000, 3200000, 64000000 };

#endif
