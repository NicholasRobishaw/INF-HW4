# Makefile

my_program: main.o queries_nw.o queries_bl.o
	g++ -std=c++17 main.o Queries_NW.o Queries_BL.o -o my_program

main.o: main.cpp main.h
	g++ -std=c++17 -c -Wall -Wextra main.cpp

queries_nw.o: Queries_NW.cpp Queries_NW.h
	g++ -std=c++17 -c -Wall -Wextra Queries_NW.cpp

queries_bl.o: Queries_BL.cpp Queries_BL.h
	g++ -std=c++17 -c -Wall -Wextra Queries_BL.cpp

clean:
	rm -f *.o my_program
