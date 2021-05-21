import os
import sys
import matplotlib.pyplot as plt
import multiprocessing
import time


# Fix word size
# Vary number of iterations and number of cores
fixedWordSize = 10
numIterations = [100000, 200000, 300000]

# Find number of cpu cores in system
numOfCores = multiprocessing.cpu_count()


print("\n\nRunning for fixed word size: ", fixedWordSize, "\n")


for numIteration in numIterations:
    print("-" * 50)
    print("Number of Iterations: ", numIteration, "\n")

    basic = "./Basic.o " + str(fixedWordSize) + " " + str(numIteration) 
    os.system(basic)

    # Number of threads applicable only for parallel versions
    for numThreads in range(1, numOfCores+1):
        print("\tNumber of threads: ", numThreads, "\n")
        
        openmp = "./Openmp.o " + str(fixedWordSize) + " " + str(numIteration) + " " + str(numThreads);
        os.system(openmp)

        # cuda = "./Cuda.o " + str(fixedWordSize) + " " + str(numIteration)
        # os.system(cuda)
        #print("\n")


    # Wait for a few seconds for resources to get freed
    time.sleep(10)
