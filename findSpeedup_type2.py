f_basic = open("./Times/basic_times.txt")
f_openmp = open("./Times/openmp_times.txt")

basic = f_basic.read().split('\n')
openmp = f_openmp.read().split('\n')

basic_sorted_by_wordsize = dict()
openmp_sorted_by_wordsize = dict()


for j in basic:
	toks = j.split(":")
	if(len(toks)>2):
		basic_sorted_by_wordsize[int(toks[0])] = float(toks[2])


for j in openmp:
	toks = j.split(":")
	if(len(toks)>3):
		if(int(toks[0]) not in openmp_sorted_by_wordsize):
			openmp_sorted_by_wordsize[int(toks[0])] = dict()

		openmp_sorted_by_wordsize[int(toks[0])][int(toks[2])] = float(toks[3])


f_basic.close()
f_openmp.close()

speedups = {}
fixednumIterations = basic[0].split(":")[1];


# Find average for fixed number of iterations
# Vary word length and number of cores per each

for wordsize in basic_sorted_by_wordsize:

	time_basic = basic_sorted_by_wordsize[wordsize]


	openmp_measurements = openmp_sorted_by_wordsize[wordsize]

	for num_cores in openmp_measurements:

		time_openmp = openmp_measurements[num_cores]

		speedup = time_basic/time_openmp

		if(wordsize not in speedups):
			speedups[wordsize] = dict()

		speedups[wordsize][num_cores] = speedup


print("Speedups for fixed number of iterations = ", fixednumIterations)

for i in speedups:
	for j in speedups[i]:
		print("Speedup for word of size ", i, ", and number of cores ", j, " is: ", speedups[i][j])

	
