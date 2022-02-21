#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_1BYTE 0x80

#include <stddef.h>
#include <stdint.h>
#include <stdalign.h>
#include <x86intrin.h>

static const alignas(128) uint64_t ascii_mask_16_byte[2] = {ASCII_MASK_8BYTE, ASCII_MASK_8BYTE};

static int
string_is_ascii(char * string, size_t length) {
    size_t n = length;
    char * char_ptr = string;
    typedef __m128i longword;
    longword *longword_ascii_mask = (longword *)ascii_mask_16_byte;

    // The first loop aligns the memory address. Char_ptr is cast to a size_t
    // to return the memory address. longword is 8 bytes long, and the processor
    // handles this better when its address is a multiplier of 8. This loops
    // handles the first few bytes that are not on such a multiplier boundary.
    while ((size_t)char_ptr % sizeof(longword) && n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    longword *longword_ptr = (longword *)char_ptr;
    while (n >= sizeof(longword)) {
        if (!_mm_testz_si128(*longword_ptr, *longword_ascii_mask)){
            return 0;
        }
        longword_ptr += 1;
        n -= sizeof(longword);
    }
    char_ptr = (char *)longword_ptr;
    while (n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    return 1;
}
