#include "bloomfilter.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <bitset>
#include <cstring>
#include <ctime>
#include <omp.h>
#include <inttypes.h>
#include <iomanip>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cstdio>

using namespace std;

#define BIG_CONSTANT(x) (x)
#define ROTL64(x,y) rotl64(x,y)
#define ROTL642(x,y) rotl642(x,y)

#define	FORCE_INLINE inline __attribute__((always_inline))

#define BIT_ARRAY_SIZE  1000
#define SEED_VALUE_1 27
#define SEED_VALUE_2 58
#define SEED_VALUE_3 99

const int MAX = 26;



__device__ inline uint64_t rotl64(uint64_t x, int8_t r){
  return (x << r) | (x >> (64 - r));
}

inline uint64_t rotl642(uint64_t x, int8_t r){
  return (x << r) | (x >> (64 - r));
}

__device__ FORCE_INLINE uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}

FORCE_INLINE uint64_t fmix642 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}

__device__ FORCE_INLINE uint64_t getblock64 ( const uint64_t * p, int i )
{
  return p[i];
}


FORCE_INLINE uint64_t getblock642 ( const uint64_t * p, int i )
{
  return p[i];
}

__device__ void MurmurHash3_x64_128(const void* key, const int len, const uint32_t seed, uint64_t* hash){

  
  const uint8_t* data = (const uint8_t*)key;
  
  const int nblocks = len/16;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

  //------------
  // body

  const uint64_t *blocks = (const uint64_t *)(data);

  
  for(int i = 0; i < nblocks; i++){


    uint64_t k1 = getblock64(blocks,i*2+0);
    uint64_t k2 = getblock64(blocks,i*2+1);

    k1 *= c1;
    k1  = ROTL64(k1,31);
    k1 *= c2;
    h1 ^= k1;

    h1 = ROTL64(h1,27);
    h1 += h2;
    h1 = h1*5+0x52dce729;

    k2 *= c2;
    k2  = ROTL64(k2,33);
    k2 *= c1;
    h2 ^= k2;

    h2 = ROTL64(h2,31);
    h2 += h1;
    h2 = h2*5+0x38495ab5;
  }
  

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  

  uint64_t k1 = 0;
  uint64_t k2 = 0;

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

void MurmurHash3_x64_1282(const void* key, const int len, const uint32_t seed, uint64_t* hash){

  
  const uint8_t* data = (const uint8_t*)key;
  
  const int nblocks = len/16;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

  //------------
  // body

  const uint64_t *blocks = (const uint64_t *)(data);

  
  for(int i = 0; i < nblocks; i++){


    uint64_t k1 = getblock642(blocks,i*2+0);
    uint64_t k2 = getblock642(blocks,i*2+1);

    k1 *= c1;
    k1  = ROTL642(k1,31);
    k1 *= c2;
    h1 ^= k1;

    h1 = ROTL642(h1,27);
    h1 += h2;
    h1 = h1*5+0x52dce729;

    k2 *= c2;
    k2  = ROTL642(k2,33);
    k2 *= c1;
    h2 ^= k2;

    h2 = ROTL642(h2,31);
    h2 += h1;
    h2 = h2*5+0x38495ab5;
  }
  

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch(len & 15){
    case 15: k2 ^= ((uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
    case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
             k2 *= c2; k2  = ROTL642(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
    case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
    case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
    case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
    case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
    case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
    case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
    case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
             k1 *= c1; k1  = ROTL642(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix642(h1);
  h2 = fmix642(h2);

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

__global__ void insertInHashTable(int* HashTable, char* key, int length){
  
  // Calculate 3 hashes and insert
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if(idx == 0){
  uint64_t hash1[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_1, hash1);
  int bit1 = (hash1[0] % BIT_ARRAY_SIZE + hash1[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
  printf("%d", bit1);
  HashTable[bit1] = 1;
  
  }

  if(idx == 1){
  uint64_t hash2[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_2, hash2);
  int bit2 = (hash2[0] % BIT_ARRAY_SIZE + hash2[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
  HashTable[bit2] = 1;
  //printf("bit 2: %d", bit2);
  }

  if(idx == 2){
  uint64_t hash3[2];
  MurmurHash3_x64_128(key, length, SEED_VALUE_3, hash3);
  int bit3 = (hash3[0] % BIT_ARRAY_SIZE + hash3[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
  HashTable[bit3] = 1;
  //printf("bit 3: %d", bit3);
  }
}


void checkIfPresent(int* HashTable, char* key, int length){
  
//   // Calculate 3 hashes and check bit
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;


//   if(idx == 0){
//   printf("Entering2\n");
  uint64_t hash1[2];
  MurmurHash3_x64_1282(key, length, SEED_VALUE_1, hash1);
  int bit1 = (hash1[0] % BIT_ARRAY_SIZE + hash1[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
//   printf("bit 1: %d", bit1);
// }

//   if(idx == 1){
  uint64_t hash2[2];
  MurmurHash3_x64_1282(key, length, SEED_VALUE_2, hash2);
  int bit2 = (hash2[0] % BIT_ARRAY_SIZE + hash2[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
//   printf("bit 2: %d", bit2);
// }

//   if(idx == 2){
  uint64_t hash3[2];
  MurmurHash3_x64_1282(key, length, SEED_VALUE_3, hash3);
  int bit3 = (hash3[0] % BIT_ARRAY_SIZE + hash3[1] % BIT_ARRAY_SIZE) % BIT_ARRAY_SIZE;
//   printf("bit 3: %d", bit3);
//   }

//   // __syncthreads();
  if(HashTable[bit1] == 1 && HashTable[bit2] == 1 && HashTable[bit3] == 1){
    printf("%s might be present\n", key);
  }
  else{
    printf("%s is definitely not present\n", key);
  }
}

int main(){

  time_t start, end;

  time(&start);

  int* HashTable = (int*)calloc(BIT_ARRAY_SIZE, sizeof(int));

  int* d_HashTable;
  cudaMalloc((void**)&d_HashTable, BIT_ARRAY_SIZE*sizeof(int));
  cudaMemset(d_HashTable, 0, BIT_ARRAY_SIZE*sizeof(int));
  cout << d_HashTable[0];
  //  for(int i=0; i<BIT_ARRAY_SIZE; i++){
  //   cout << d_HashTable[i];
  // }

  // copy to DEVICE
  cudaMemcpy(d_HashTable, HashTable, BIT_ARRAY_SIZE*sizeof(int), cudaMemcpyHostToDevice);

  cout << sizeof(d_HashTable) << endl;

  string str = "thefirsthashfunctionjdfsjldkjklsjlfjdslkjflsjwjadjsaijadijasjdoiajsdoiannewfjsoinsfnoiesfinosen";
  int len = str.length();
  cout << len << endl;
  char *cstr = new char[str.length() + 1];
  strcpy(cstr, str.c_str());
  
  
  insertInHashTable<<<1,3>>>(d_HashTable, cstr, len);
  cudaDeviceSynchronize();

  cudaMemcpy(HashTable, d_HashTable, BIT_ARRAY_SIZE*sizeof(int), cudaMemcpyDeviceToHost);
 
  cout << endl;
  checkIfPresent(HashTable, cstr, len);
  // /cudaDeviceSynchronize();
  cudaFree(d_HashTable);
  

  time(&end);
  double timeTaken = double(end - start);
  //cout << "Time taken for inserting " << numIterations <<  " records in unparallelized version: " << fixed << timeTaken << setprecision(9);
  //cout << "s" << endl;

}