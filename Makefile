CCPP = g++
CPPFLAGS = -std=c++17 -pedantic -Wall -fopenmp

gal:  gal.o
	$(CCPP) $(CPPFLAGS) -o gal gal.o -lboost_graph -fopenmp

clean:
	rm -f *.o
