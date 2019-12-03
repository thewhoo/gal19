CCPP = g++
CPPFLAGS = -std=c++11 -pedantic -Wall

gal:  gal.o
	$(CCPP) $(CPPFLAGS) -o gal gal.o -lboost_graph

clean:
	rm -f *.o
