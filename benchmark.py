import os
import sys
import matplotlib.pyplot as plt

def buildDict(logs_type):
    type_dict = {}
    for logs in logs_type:
        arr = logs.split(":")
        wordLen = int(arr[0])
        time = float(arr[2])

        if wordLen not in type_dict:
            type_dict[wordLen] = []
        type_dict[wordLen].append(time)
    return type_dict
    


numIterations = [10, 100, 1000, 10000, 100000]
numLengths = [5, 10, 32, 64, 70]

for length in numLengths:
    for iteration in numIterations:
        #print("running: ", length, " ", iteration, "\n")
        basic = "./Basic.o " + str(length) + " " + str(iteration)
        os.system(basic)
        openmp = "./Openmp.o " + str(length) + " " + str(iteration)
        os.system(openmp)
        cuda = "./Cuda.o " + str(length) + " " + str(iteration)
        os.system(cuda)
        #print("\n")

#Read  text files and plot here
#Basic plot
#Multiple lines: each line represent word length, for each iteration how much time


f_basic = open("./Times/basic_times.txt")
f_openmp = open("./Times/openmp_times.txt")
f_cuda = open("./Times/cuda_times.txt")

basic = f_basic.read().split('\n')
openmp = f_openmp.read().split('\n')
cuda = f_cuda.read().split('\n')

f_basic.close()
f_openmp.close()
f_cuda.close()

basic = basic[:-1]
openmp = openmp[:-1]
cuda = cuda[:-1]

basic_dict = buildDict(basic)
openmp_dict = buildDict(openmp)
cuda_dict = buildDict(cuda)


x = numIterations
for wordLen in numLengths:
    y_basic = basic_dict[wordLen]
    y_openmp = openmp_dict[wordLen]
    y_cuda = cuda_dict[wordLen]

    plt.figure()
    plt.plot(x, y_basic, label = "Basic")
    #plt.plot(x, y_openmp, label = "OpenMP")
    plt.plot(x, y_cuda, label="CUDA")
    plt.legend()
    plt.title("Word Length: " + str(wordLen))
    plt.xlabel("Number of Iterations")
    plt.ylabel("Time in ms")
plt.show()



