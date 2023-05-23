#include "emmintrin.h"
#include <stdint.h>
size_t count_newlines(const char *text, size_t text_size) {
    size_t count = 0;
    size_t bytes_to_align = (size_t)text / sizeof(__m128i);
    for (size_t i=0; i < bytes_to_align; i++) {
        if (text[i] == '\n') {
            count += 1;
        }
    }
    const char *cursor = text + bytes_to_align;
    const char *end_ptr = text + text_size;
    __m128i newlines = _mm_set1_epi8('\n');
    while (cursor < (end_ptr - sizeof(__m128i))) {
        size_t chunks_remaining = (end_ptr - cursor) / sizeof(__m128i);
        if (chunks_remaining > UINT8_MAX) {
            // This is the maximum that we can accumalate using the method below
            chunks_remaining = UINT8_MAX;
        }
        __m128i accumulator = _mm_setzero_si128();
        for (size_t i=0; i<chunks_remaining; i++) {
            __m128i data = _mm_load_si128((__m128i *)cursor);
            __m128i eq = _mm_cmpeq_epi8(data, newlines);
            // cmp_eq sets the most significant bit if a newline is found. 
            // convert this to a one.
            // 0b1000_0000 >> 7 = 0b0000_0001 
            __m128i ones = _mm_srli_si128(eq, 7);
            accumulator = _mm_adds_epu8(accumulator, ones);
            cursor += sizeof(__m128i);
        }
        static uint8_t accumulated_counts[sizeof(__m128i)];
        _mm_storeu_si128((__m128i *)&accumulated_counts, accumulator);
        for (size_t i=0; i < sizeof(__m128i); i++) {
            count += accumulated_counts[i];
        }
    }
    while (cursor < end_ptr) {
        if (*cursor == '\n') {
            count += 1;
        }
        cursor += 1;
    }
    return count;
}