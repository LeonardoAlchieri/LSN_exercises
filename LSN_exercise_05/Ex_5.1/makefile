CC = g++
CFLAGS = -Wall -O3 --std=c++11

hidrogen_metropolis.exe : main.o random.o
	$(CC) random.o main.o -o hidrogen_metropolis
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe seed.out *.txt
