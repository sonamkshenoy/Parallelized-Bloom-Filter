# Parallelized-Bloom-Filter

An accelerated bloom filter



# Commands to execute  
To plot comparison between openmp, cuda and serial versions
```
make
```

To plot speedups for i) Fixed number of iterations ii) Varying word size iii) Varying number of cores
```
make  
python3 benchmark_type2.py  
python3 findSpeedup_type2.py
```

To plot speedups for i) Fixed word size ii) Varying number of iterations iii) Varying number of cores
```
make  
python3 benchmark_type3.py  
python3 findSpeedup_type3.py
```


## References
* [Murmur Hash 3 implementation](https://github.com/aappleby/smhasher) by Austin Appleby.  
