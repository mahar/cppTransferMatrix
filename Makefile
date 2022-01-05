all: calculation

calculation: calculation.o tmm.o
	g++ calculation.o tmm.o -o calculation

tmm.o: tmm.cpp tmm.h
	g++ -c tmm.cpp 

calculation.o: calculation.cpp
	g++ -c calculation.cpp

clean:
	rm *.o calculation