default : extractLines

clean :
	rm -f io_png.o levelLine.o main.o

io_png.o : io_png.c io_png.h
	gcc -c -g io_png.c
levelLine.o : levelLine.cpp levelLine.h
	g++ -c -g levelLine.cpp
fill_curve.o : fill_curve.cpp levelLine.h
	g++ -c -g fill_curve.cpp
lltree.o : lltree.cpp lltree.h levelLine.h
	g++ -c -g lltree.cpp
main.o : main.cpp fill_curve.h lltree.h levelLine.h cmdLine.h io_png.h
	g++ -c -g main.cpp

extractLines : main.o fill_curve.o llfill.o levelLine.o io_png.o
	g++ -o extractLines main.o fill_curve.o lltree.o levelLine.o io_png.o -lpng
