#include <string.h>

static inline int contains_newline(char *string, size_t length) {
    return memchr(string, '\n', length) != NULL;
}
