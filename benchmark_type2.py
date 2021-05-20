import os
import sys
import matplotlib.pyplot as plt
import multiprocessing


# Fix number of iterations
# Vary word size and number of cores
fixedNumIterations = 50000
wordSizes = [16, 32, 70]

# Find number of cpu cores in system
numOfCores = multiprocessing.cpu_count()


print("\n\nRunning for fixed number of Iterations: ", fixedNumIterations, "\n")


for wordSize in wordSizes:
    print("-" * 50)
    print("Word size: ", wordSize, "\n")

    basic = "./Basic.o " + str(wordSize) + " " + str(fixedNumIterations) 
    os.system(basic)

    # Number of threads applicable only for parallel versions
    for numThreads in range(1, numOfCores+1):
        print("\tNumber of threads: ", numThreads, "\n")
        
        openmp = "./Openmp.o " + str(wordSize) + " " + str(fixedNumIterations) + " " + str(numThreads);
        os.system(openmp)

        # cuda = "./Cuda.o " + str(wordSize) + " " + str(fixedNumIterations)
        # os.system(cuda)
        #print("\n")
