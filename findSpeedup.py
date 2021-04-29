f_basic = open("./Times/basic_times.txt")
f_openmp = open("./Times/cuda_times.txt")

basic = f_basic.read().split('\n')
openmp = f_openmp.read().split('\n')

f_basic.close()
f_openmp.close()

speedups = {}
numIterations = 0;

# Find average for each word length
for i in zip(basic, openmp):
	tokens = i[0].split(":")
	if(len(tokens) == 3):
		tokens_o = i[1].split(":")

		wordsize = int(tokens[0])

		time = float(tokens[2])
		time_o = float(tokens_o[2])

		speedup = time/time_o

		if(wordsize not in speedups):
			speedups[wordsize] = 0
			numIterations = 0
		
		numIterations += 1

		speedups[wordsize] += speedup


for j in speedups:
	print("Speedup for word of size ", j, " is: ", speedups[j]/numIterations)


# Find average for each number of insertions/iterations

speedups_iter = {}
numwords = 0;

for i in zip(basic, openmp):
	tokens = i[0].split(":")
	if(len(tokens) == 3):
		tokens_o = i[1].split(":")

		iterationsize = int(tokens[1])

		time = float(tokens[2])
		time_o = float(tokens_o[2])

		speedup = time/time_o

		if(iterationsize not in speedups_iter):
			speedups_iter[iterationsize] = 0
		
		speedups_iter[iterationsize] += speedup

print(speedups_iter)

for j in speedups_iter:
	print("Speedup for number of insertions = ", j, " is: ", speedups_iter[j]/numIterations)



	
