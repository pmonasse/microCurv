default : extractLines

clean :
	rm -f io_png.o levelLine.o main.o

io_png.o : io_png.c io_png.h
	gcc -c -g io_png.c
levelLine.o : levelLine.cpp levelLine.h
	g++ -c -g levelLine.cpp
main.o : main.cpp levelLine.h io_png.h
	g++ -c -g main.cpp

extractLines : main.o levelLine.o io_png.o
	g++ -o extractLines -lpng main.o levelLine.o io_png.o
