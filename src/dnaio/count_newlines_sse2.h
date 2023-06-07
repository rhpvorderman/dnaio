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
        /* Use a vector of 16 uint8 integers to count newlines. This makes it
           easy to accumulate the result of _mm_cmp_eq_epi8 which also reports 
           it results as 8-bit integers. 
           Only 255 vectors can be counted this way. If more are counted one of 
           the  accumulators might get saturated and the total count will be 
           incorrect. 
           The outer loop ensures it keeps running until all 255x16 chunks are 
           counted.
           */
        __m128i uint8x16_accumulator = _mm_setzero_si128();
        size_t chunks_remaining = (end_ptr - cursor) / sizeof(__m128i);
        if (chunks_remaining > UINT8_MAX) {
            chunks_remaining = UINT8_MAX;
        }
        for (size_t i=0; i<chunks_remaining; i++) {
            __m128i uint8x16_data = _mm_load_si128((__m128i *)cursor);
            __m128i uint8x16_eq = _mm_cmpeq_epi8(uint8x16_data, _mm_set1_epi8('\n'));
            // cmp_eq sets all bits if a newline is found.
            // convert this to a one. 255 - 254 = 1
            __m128i uint8x16_ones = _mm_subs_epu8(uint8x16_eq, _mm_set1_epi8(254));
            uint8x16_accumulator = _mm_adds_epu8(uint8x16_accumulator, uint8x16_ones);
            cursor += sizeof(__m128i);
        }
        static uint8_t accumulated_counts[sizeof(__m128i)];
        _mm_storeu_si128((__m128i *)&accumulated_counts, uint8x16_accumulator);
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
