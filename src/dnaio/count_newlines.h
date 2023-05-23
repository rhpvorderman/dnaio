#include <stddef.h>

static inline size_t count_newlines(const char *text, size_t text_size) {
    size_t count = 0;
    const char *cursor = text;
    const char *end_ptr = text + text_size;
    while (1) {
        cursor = memchr(cursor, '\n', end_ptr - cursor);
        if (cursor == NULL) {
            return count;
        }
        count += 1;
        cursor += 1;  // Hop over newline
    }
}
