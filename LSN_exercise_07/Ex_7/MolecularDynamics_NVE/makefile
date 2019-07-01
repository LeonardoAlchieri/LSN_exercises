CC = g++
CFLAGS = -Wall -O3 --std=c++11

simulation.exe : main.o random.o
	$(CC) random.o main.o -o simulation.exe
main.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o main.o $(CFLAGS)
random.o : random_generator/random.cpp random_generator/random.h
	$(CC) -c random_generator/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulation.exe
