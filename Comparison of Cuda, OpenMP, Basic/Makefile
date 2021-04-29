all:
	rm -f ./Times/*.txt
	rm -f Basic.o
	rm -f Openmp.o
	rm -f Cuda.o
	g++ Basicversion.cpp -o Basic.o
	g++ OpenMPversion.cpp -fopenmp -o Openmp.o
	nvcc Cudaversion.cu -o Cuda.o
	python3 benchmark.py

clean:
	rm *.o

cleanlogs:
	rm ./Times/*.txt

benchmark:
	python3 benchmark.py

cleanandbenchmark:
	rm ./Times/*.txt
	python3 benchmark.py

gen:
	g++ Basicversion.cpp -o Basic.o
	g++ OpenMPversion.cpp -o Openmp.o
	nvcc Cudaversion.cu -o Cuda.o
