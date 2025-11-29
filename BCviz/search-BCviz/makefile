GCC=g++ -O3
PTHREADS = -lpthread
OBJS=main.o ooMBE.o bgraph.o
EXECUTABLE = MBS

all: $(EXECUTABLE)

MBS: $(OBJS)
	$(GCC) -o MBS $(OBJS) $(PTHREADS)
main.o: main.cpp bgraph.h ooMBE.h
	$(GCC) -c main.cpp -g
ooMBE.o: ooMBE.cpp ooMBE.h bgraph.h
	$(GCC) -c ooMBE.cpp -g
ooMBE.h: bgraph.h
bgraph.o: bgraph.cpp bgraph.h
	$(GCC) -c bgraph.cpp -g
bgraph.h:

cleanall:
	rm -f core* *.o $(EXECUTABLE)

clean:
	rm -f core* *.o