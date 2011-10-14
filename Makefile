OSXARGS = -I/sw/include -L/sw/lib
CFLAGS = -O2
all : wallop scallop

wallop.o : wallop.c 
	gcc -c wallop.c $(CFLAGS) $(OSXARGS) `pkg-config --cflags --libs gtk+-2.0`

op.o : op.c 
	gcc -c op.c $(CFLAGS) $(OSXARGS) `pkg-config --cflags --libs gtk+-2.0`

#io.o : io.c 
#	gcc -c io.c $(CFLAGS) $(OSXARGS) `pkg-config --cflags --libs gtk+-2.0`

fatgraph.o : fatgraph.c 
	gcc -c fatgraph.c $(CFLAGS) $(OSXARGS) `pkg-config --cflags --libs gtk+-2.0`

vector2d.o : vector2d.c 
	gcc -c vector2d.c $(CFLAGS) $(OSXARGS) `pkg-config --cflags --libs gtk+-2.0`

wallop : wallop.o op.o vector2d.o fatgraph.o
	gcc $(CFLAGS) $(OSXARGS) -o wallop wallop.o op.o vector2d.o fatgraph.o `pkg-config --cflags --libs gtk+-2.0`

.PHONY : scallop
scallop :
	cd scallop; make

clean : 
	rm wallop
	rm *.o
	cd scallop; make clean

