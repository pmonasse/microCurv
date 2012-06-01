default : extractLines smoothLines

clean :
	rm -f *.o extractLines smoothLines

io_png.o : io_png.c io_png.h
	gcc -c -g io_png.c
levelLine.o : levelLine.cpp levelLine.h
	g++ -c -g levelLine.cpp
fill_curve.o : fill_curve.cpp levelLine.h
	g++ -c -g fill_curve.cpp
lltree.o : lltree.cpp lltree.h levelLine.h
	g++ -c -g lltree.cpp
gass.o : gass.cpp gass.h
	g++ -c -g gass.cpp
main.o : main.cpp gass.h fill_curve.h lltree.h levelLine.h cmdLine.h io_png.h
	g++ -c -g main.cpp

main_gass.o : main_gass.cpp gass.h levelLine.h cmdLine.h
	g++ -c -g main_gass.cpp

extractLines : main.o gass.o fill_curve.o lltree.o levelLine.o io_png.o
	g++ -o extractLines main.o gass.o fill_curve.o lltree.o levelLine.o io_png.o -lpng

smoothLines : main_gass.o gass.o levelLine.o
	g++ -o smoothLines main_gass.o gass.o levelLine.o
