OSXARGS = -I/sw/include -L/sw/lib
all : wallop.c op.c io.c vector2d.c
	gcc -O2 $(OSXARGS) -o wallop wallop.c op.c vector2d.c fatgraph.c `pkg-config --cflags --libs gtk+-2.0`

clean : 
	rm wallop
