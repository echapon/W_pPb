all: wresults

wresults: wresults.o statchannel.o
	g++ wresults.o statchannel.o -o wresults `root-config --libs`

wresults.o: wresults.cpp
	g++ -c wresults.cpp -o wresults.o `root-config --cflags`

statchannel.o: statchannel.C statchannel.h lumierror.h
	g++ -c statchannel.C -o statchannel.o `root-config --cflags`

clear:
	rm *.o wresults
