#define _MURMURHASH3_H_
#include <cstdint>

typedef unsigned int uint32_t;
typedef unsigned char uint8_t;
typedef unsigned __int64 uint64_t;

// Arguments (word to be hashed, length of the key, a seed, the array of two 64 integers storing final hash value)
void MurmurHash3_x64_128(const void* key, int len, uint32_t seed, uint64_t* hash);
