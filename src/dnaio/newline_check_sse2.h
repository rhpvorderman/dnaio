#include <string.h>
#include <emmintrin.h>
#include <stdint.h>

static inline int contains_newline(const char *string, size_t length) {
    __m128i allcomps = _mm_setzero_si128();
    __m128i all_newlines = _mm_set1_epi8('\n');
    size_t remaining = length;
    const char *cursor = string; 
    __m128i current;
    uint64_t comp = 0;
    uint8_t tmp;
    while (remaining > 0 && ((size_t)cursor % sizeof(__m128i))) {
        tmp = ((((*cursor ^ '\n') - 1) & ~'\n') & (1 << 7));
        comp |= tmp;
        remaining -= 1;
        cursor +=1;
    }
    while(remaining > sizeof(__m128i)) {
        current = _mm_load_si128((__m128i *)cursor);
        allcomps = _mm_or_si128(allcomps, _mm_cmpeq_epi8(current, all_newlines));
        remaining -= sizeof(__m128i);
        cursor += sizeof(__m128i);
    }
    comp |= _mm_movemask_epi8(allcomps);
    while (remaining) {
        uint8_t tmp = ((((*cursor ^ '\n') - 1) & ~'\n') & (1 << 7));
        comp |= tmp;
        remaining -= 1;
        cursor +=1;
    }
    return comp != 0;
}
