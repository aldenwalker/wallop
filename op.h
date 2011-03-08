#ifndef OP_H
#define OP_H


#include "vector2d.h"
#include "fatgraph.h"

typedef struct {
  vector2d loc;
  int num_edges;
  int* dest_verts;
  int* edges;
  double* theta;
} op_vert;



void fatgraph_optimize_drawing(fatgraph* fg, 
                               double screen_width,
                               double screen_height, 
                               int optimize_positions_also);






#endif
