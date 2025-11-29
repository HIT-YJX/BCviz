GCC=g++ -std=c++11 -march=native -mavx -O3
OBJS=main.o BCVIZ.o Bigraph.o MB.o

BCviz: $(OBJS)
	$(GCC) -o BCviz $(OBJS)
clean:
	rm $(OBJS); rm BCviz

MB.o: MB.cpp MB.h
	$(GCC) -c MB.cpp -lpthread -g -DDISCRETE
MB.h:
BCVIZ.o: BCVIZ.cpp BCVIZ.h
	$(GCC) -c BCVIZ.cpp -lpthread -g -DDISCRETE
BCVIZ.h: Bigraph.h MB.h
Bigraph.o: Bigraph.cpp Bigraph.h
	$(GCC) -c Bigraph.cpp -lpthread -g -DDISCRETE
Bigraph.h:
main.o: main.cpp BCVIZ.h
	$(GCC) -c main.cpp -lpthread -g -DDISCRETE
