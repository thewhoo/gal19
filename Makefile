CCPP = g++
CPPFLAGS = -std=c++17 -Wall -pedantic -Wno-sign-compare -fopenmp

gal:  gal.o
	$(CCPP) $(CPPFLAGS) -o gal gal.o -lboost_graph

clean:
	rm -f *.o gal
