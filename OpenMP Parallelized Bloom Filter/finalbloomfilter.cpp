#define _MURMURHASH3_H_
#include <cstdint>

typedef unsigned int uint32_t;
typedef unsigned char uint8_t;
typedef __uint64_t uint64_t;

// Arguments (word to be hashed, length of the key, a seed, the array of two 64 integers storing final hash value)
void MurmurHash3_x64_128(const void* key, int len, uint32_t seed, uint64_t* hash);

// # include "bloomfilter.h"
#include <stdlib.h>
#include <iostream>
#include <semaphore.h>
#include <vector>
#include <bitset>
#include <cstring>
#include <ctime>
#include <omp.h>
#include <inttypes.h>
#include <iomanip>
#include <chrono>

using namespace std;

#define BIG_CONSTANT(x) (x)
#define ROTL64(x,y) rotl64(x,y)
#define	FORCE_INLINE inline __attribute__((always_inline))

#define BIT_ARRAY_SIZE  100000
#define SEED_VALUE_1 27
#define SEED_VALUE_2 58
#define SEED_VALUE_3 99

const int MAX = 26;
sem_t semaphore;

inline uint64_t rotl64(uint64_t x, int8_t r){
  return (x << r) | (x >> (64 - r));
}

FORCE_INLINE uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}


FORCE_INLINE uint64_t getblock64 ( const uint64_t * p, int i )
{
  return p[i];
}

void MurmurHash3_x64_128(const void* key, const int len, const uint32_t seed, uint64_t* hash, uint64_t k1, uint64_t k2, uint64_t k3, uint64_t k4){

  const uint8_t* data = (const uint8_t*)key;
  const int nblocks = len/16;
  uint64_t h1 = seed;
  uint64_t h2 = seed;
  uint64_t c1;
  uint64_t c2;
  c1 = BIG_CONSTANT(0x87c37b91114253d5);
  c2 = BIG_CONSTANT(0x4cf5ad432745937f);
  const uint64_t *blocks = (const uint64_t *)(data);

  h1 ^= k1;

  h1 = ROTL64(h1,27);
  h1 += h2;
  h1 = h1*5+0x52dce729;

  // cout << "h1:        " << h1 << "\n";


  h2 ^= k2;

  h2 = ROTL64(h2,31);
  h2 += h1;
  h2 = h2*5+0x38495ab5;

  // cout << "h2:        " << h2 << "\n";

  h1 ^= k3;

  h1 = ROTL64(h1,27);
  h1 += h2;
  h1 = h1*5+0x52dce729;


  h2 ^= k4;

  h2 = ROTL64(h2,31);
  h2 += h1;
  h2 = h2*5+0x38495ab5;


  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  // uint64_t 
  k1 = 0;
  //uint64_t 
  k2 = 0;

  switch(len & 15){
    case 15: k2 ^= ((uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
    case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
             k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
    case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
    case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
    case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
    case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
    case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
    case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
    case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
             k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix64(h1);
  h2 = fmix64(h2);

  h1 += h2;
  h2 += h1;

  ((uint64_t*)hash)[0] = h1;
  ((uint64_t*)hash)[1] = h2;

}

string genRandomString(int n) 
{ 
    char alphabet[MAX] = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 
                          'h', 'i', 'j', 'k', 'l', 'm', 'n',  
                          'o', 'p', 'q', 'r', 's', 't', 'u', 
                          'v', 'w', 'x', 'y', 'z' }; 
  
    string res = ""; 

    for (int i = 0; i < n; i++)  
        res = res + alphabet[rand() % MAX]; 
      
    return res; 
} 

void insertInHashTable(bitset<BIT_ARRAY_SIZE>& HashTable, char* key, int length){
  
  // Calculate 3 hashes and insert
  uint64_t hash1[2];
  uint64_t hash2[2];
  uint64_t hash3[2];
  int bit1, bit2, bit3;




  const uint8_t* data = (const uint8_t*)key;
  const int nblocks = length/16;
  uint64_t c1;
  uint64_t c2;
  c1 = BIG_CONSTANT(0x87c37b91114253d5);
  c2 = BIG_CONSTANT(0x4cf5ad432745937f);
  const uint64_t *blocks = (const uint64_t *)(data);
  int i = 0;
  uint64_t k1, k2, k3, k4;
  k1 = getblock64(blocks,i*2+0);
  k1 *= c1;
  k1  = ROTL64(k1,31);
  k1 *= c2;

  k2 = getblock64(blocks,i*2+1);
  k2 *= c2;
  k2  = ROTL64(k2,33);
  k2 *= c1;

  k3 = getblock64(blocks,i*2+2);
  k3 *= c1;
  k3  = ROTL64(k3,31);
  k3 *= c2;

  k4 = getblock64(blocks,i*2+3);
  k4 *= c2;
  k4  = ROTL64(k4,33);
  k4 *= c1;



MurmurHash3_x64_128(key, length, SEED_VALUE_1, hash1, k1, k2, k3, k4);
bit1 = (hash1[0] % BIT_ARRAY_SIZE + hash1[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;

MurmurHash3_x64_128(key, length, SEED_VALUE_2, hash2, k1, k2, k3, k4);
bit2 = (hash2[0] % BIT_ARRAY_SIZE + hash2[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;

MurmurHash3_x64_128(key, length, SEED_VALUE_3, hash3, k1, k2, k3, k4);
bit3 = (hash3[0] % BIT_ARRAY_SIZE + hash3[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;  


  // cout << "Bits set are: " << bit1 << "," << bit2 << " and " << bit3 << "\n";


  HashTable.set(bit1);
  HashTable.set(bit2);
  HashTable.set(bit3);


  //cout << "Set bits: " << bit1 << ", " << bit2 << ", " << bit3 << "\n";
}


/*
void checkIfPresent(bitset<BIT_ARRAY_SIZE> HashTable, char* key, int length){
  
  // Calculate 3 hashes and check bit

  uint64_t hash1[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_1, hash1);
  int bit1 = (hash1[0] % BIT_ARRAY_SIZE + hash1[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;

  uint64_t hash2[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_2, hash2);
  int bit2 = (hash2[0] % BIT_ARRAY_SIZE + hash2[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;

  uint64_t hash3[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_3, hash3);
  int bit3 = (hash3[0] % BIT_ARRAY_SIZE + hash3[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
  
  if(HashTable.test(bit1) == 1 && HashTable.test(bit2) == 1 && HashTable.test(bit3) == 1){
    cout << key << " might be present" << "\n";
  }
  else{
    cout << key << " is definitely not present" << "\n";
  }
}*/

int main(){

    int lenOfWord = 32;
    string str;
    int numIterations = 100000;

    vector<string> wordsToInsert;
    for(int i = 0; i < numIterations; i++){
        str = genRandomString(lenOfWord);
        wordsToInsert.push_back(str);
    }


    char* cstr;
    bitset<BIT_ARRAY_SIZE> HashTable; 

    auto t_start = std::chrono::high_resolution_clock::now(), t_end = std::chrono::high_resolution_clock::now();
    

    #pragma omp parallel
    {
      #pragma omp single
      {
        t_start = std::chrono::high_resolution_clock::now();
      }
      
      #pragma omp for
      for(int i = 0; i < numIterations; ++i){
          str = wordsToInsert[i];
          cstr = new char[lenOfWord + 1];
          strcpy(cstr, str.c_str()); 
          insertInHashTable(HashTable, cstr, lenOfWord);
      }
     
      #pragma omp single
      {
        t_end = std::chrono::high_resolution_clock::now();
      }
    }
  
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();

  cout << "Time taken for inserting " << numIterations <<  " records in Openmp parallelized version: " << fixed << elapsed_time_ms << setprecision(9);
  cout << " ms" << endl;
}
