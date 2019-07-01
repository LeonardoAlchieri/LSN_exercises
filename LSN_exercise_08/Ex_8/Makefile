CC = g++
CFLAGS = -Wall -O3 --std=c++11

VCM.exe : VCM.o random.o
	$(CC) random.o VCM.o -o VCM.exe
VCM.o : VCM.cpp
	$(CC) -c VCM.cpp -o VCM.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VCM.exe
