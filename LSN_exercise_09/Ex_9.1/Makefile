CC = g++
CFLAGS = -Wall -O3 --std=c++11

genetic.exe : main.o random.o 
	$(CC) random.o main.o -o genetic.exe
main.o : main.cpp classes.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o genetic.exe
