#include "emmintrin.h"
#include <stdint.h>

size_t count_newlines(const char *text, size_t text_size) {
    const char *cursor = text;
    const char *end_ptr = text + text_size;
    size_t count = 0;
    // Align cursor to __m128i boundary
    while ((cursor < end_ptr) && ((size_t)cursor % sizeof(__m128i))) {
        if (*cursor == '\n') {
            count += 1;
        }
        cursor += 1;
    }
    while (cursor < (end_ptr - sizeof(__m128i))) {
        size_t chunks_remaining = (end_ptr - cursor) / sizeof(__m128i);
        if (chunks_remaining > UINT8_MAX) {
            // This is the maximum that we can accumalate using the method below
            chunks_remaining = UINT8_MAX;
        }
        __m128i accumulator = _mm_setzero_si128();
        for (size_t i=0; i<chunks_remaining; i++) {
            __m128i data = _mm_load_si128((__m128i *)cursor);
            __m128i eq = _mm_cmpeq_epi8(data, _mm_set1_epi8('\n'));
            // cmp_eq sets all bits if a newline is found.
            // convert this to a one. 255 - 254 = 1
            __m128i ones = _mm_subs_epu8(eq, _mm_set1_epi8(254));
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
