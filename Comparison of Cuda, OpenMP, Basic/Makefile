all:
	g++ Basicversion.cpp -o Basic.o
	g++ OpenMPversion.cpp -o Openmp.o
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
