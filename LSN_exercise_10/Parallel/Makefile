CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

genetic.exe : main.o my_random.o 
	$(CC) my_random.o main.o -o genetic.exe
main.o : main.cpp classes.h metropolis.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
my_random.o : my_random.cpp my_random.h
	$(CC) -c my_random.cpp -o my_random.o $(CFLAGS)

clean :
	rm *.o genetic.exe
