CC = g++
CFLAGS = -Wall -O3 --std=c++11

qmc1d.exe : qmc1d.o random.o
	$(CC) random.o qmc1d.o -o qmc1d.exe

qmc1d.o : qmc1d.cpp constants.h functions.h
	$(CC) -c qmc1d.cpp -o qmc1d.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o qmc1d.exe
