all:
	rm -f *.txt
	g++ Basicversion.cpp -o Basic.o
	./Basic.o 16 50 >> basic.txt
	g++ OpenMPversion.cpp -fopenmp -o OpenMP.o
	./OpenMP.o 16 50 4 >> openmp.txt
	nvcc Cudaversion.cu -o Cuda.o
	./Cuda.o 16 50 >> cuda.txt
