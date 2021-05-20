f_basic = open("./Times/basic_times.txt")
f_openmp = open("./Times/openmp_times.txt")

basic = f_basic.read().split('\n')
openmp = f_openmp.read().split('\n')

print(basic)
print(openmp)

basic_sorted_by_numIterations = dict()
openmp_sorted_by_numIterations = dict()


for j in basic:
	toks = j.split(":")
	if(len(toks)>2):
		basic_sorted_by_numIterations[int(toks[1])] = float(toks[2])


for j in openmp:
	toks = j.split(":")
	if(len(toks)>3):
		if(int(toks[1]) not in openmp_sorted_by_numIterations):
			openmp_sorted_by_numIterations[int(toks[1])] = dict()

		openmp_sorted_by_numIterations[int(toks[1])][int(toks[2])] = float(toks[3])

print(basic_sorted_by_numIterations)
print(openmp_sorted_by_numIterations)


f_basic.close()
f_openmp.close()

speedups = {}
numIterations = 0;


# Find average for fixed number of iterations
# Vary word length and number of cores per each

for numIteration in basic_sorted_by_numIterations:

	time_basic = basic_sorted_by_numIterations[numIteration]


	openmp_measurements = openmp_sorted_by_numIterations[numIteration]

	for num_cores in openmp_measurements:

		time_openmp = openmp_measurements[num_cores]

		speedup = time_basic/time_openmp

		if(numIteration not in speedups):
			speedups[numIteration] = dict()

		speedups[numIteration][num_cores] = speedup


for i in speedups:
	for j in speedups[i]:
		print("Speedup for number of iterations ", i, ", and number of cores ", j, " is: ", speedups[i][j])

