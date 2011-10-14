#ifndef IO_DEF
#define IO_DEF

#include "fatgraph.h"


fatgraph* read_fatgraph_from_file(char* filename, double screen_width, double screen_height);
int write_fatgraph_to_file(fatgraph* fg, char* filename);

#endif
