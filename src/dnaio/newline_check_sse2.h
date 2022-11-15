#include <string.h>
#include <emmintrin.h>
static inline int contains_newline(const char *restrict string, size_t length) {
    __m128i allcomps = _mm_setzero_si128();
    __m128i all_newlines = _mm_set1_epi8('\n');
    size_t remaining = length;
    const char *cursor = string; 
    __m128i current;
    while(remaining > sizeof(__m128i)) {
        current = _mm_loadu_si128((__m128i *)cursor);
        allcomps = _mm_or_si128(allcomps, _mm_cmpeq_epi8(current, all_newlines));
        remaining -= sizeof(__m128i);
        cursor += sizeof(__m128i);
    }
    int comp = _mm_movemask_epi8(allcomps);
    while (remaining > 0) {
        comp |= (*cursor == '\n');
        remaining -= 1;
        cursor +=1;
    }
    return comp != 0;
}
