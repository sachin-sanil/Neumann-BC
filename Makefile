
all: main solver grid grid_int mgsolve

mgsolve: main.o Grid_int.o Grid.o Solver.o
	g++ -O3 -Wall -Winline -Wshadow main.o  Grid_int.o Grid.o Solver.o -o mgsolve

grid_int: Grid_int.h Grid_int.cpp
	g++ -c -O3 -Wall -Winline -Wshadow -std=c++11 Grid_int.cpp -o Grid_int.o

grid: Grid.h Grid.cpp
	g++ -c -O3 -Wall -Winline -Wshadow -std=c++11 Grid.cpp -o Grid.o

solver: Solver.h Solver.cpp
	g++ -c -O3 -Wall -Winline -Wshadow -std=c++11 Solver.cpp -o Solver.o

main: main.cpp Grid_int.cpp Grid_int.o
	g++ -c -O3 -Wall -Winline -Wshadow -std=c++11 main.cpp -o main.o

clean:
	rm -rf solution.txt
	rm -rf *.o
	rm mgsolve

