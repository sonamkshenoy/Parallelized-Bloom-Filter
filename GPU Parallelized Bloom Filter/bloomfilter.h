#include <cstdint>
#include <inttypes.h>

#define _MURMURHASH3_H_

typedef unsigned int uint32_t;
typedef unsigned char uint8_t;

#ifdef _WIN32
    typedef unsigned __int64 uint64_t;
#endif

#ifdef __linux__
    typedef __uint64_t uint64_t;
#endif

// Arguments (word to be hashed, length of the key, a seed, the array of two 64 integers storing final hash value)
__global__ void MurmurHash3_x64_128(const void* key, int len, uint32_t seed, uint64_t* hash);