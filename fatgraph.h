#ifndef FG_DEF
#define FG_DEF

#include "vector2d.h"

/******************************************************************************/
/* fatgraph type                                                              */
/******************************************************************************/

typedef struct {
  char* name;
  int* edges;  //in order!
  int* edges_initial;
  int num_edges;
  vector2d* bezier;
  vector2d loc;
} vert;

typedef struct {
  char* name;
  char* label_forward;
  char* label_backward;
  int start;
  int end;
} edge;

typedef struct {
  vert* verts;
  edge* edges;
  int num_verts;
  int num_edges;
} fatgraph;

void fatgraph_print(fatgraph* f);
void fatgraph_free(fatgraph* f);
void fatgraph_copy(fatgraph* F, fatgraph* G);
void fatgraph_add_vertex(fatgraph* fg, double x, double y);
void fatgraph_add_edge(fatgraph* fg, 
                       int s, 
                       int index_in_s, 
                       vector2d sbezier,
                       int e,
                       int index_in_e,
                       vector2d ebezier);
fatgraph* fatgraph_create();
int fatgraph_close_to_edge(fatgraph* fg, double x, double y);
double fatgraph_remotest(fatgraph* fg, int edge);
int fatgraph_chi(fatgraph* fg);
void fatgraph_next_edge(fatgraph* fg, int* edge, int* dir);
char* fatgraph_follow_edge(fatgraph* fg, int dir, int edge, int* edges_remaining);
char** fatgraph_boundaries(fatgraph* fg, int* numBoundaries);
fatgraph* read_fatgraph_from_file_new(char* filename, double screen_width,
                                                      double screen_height);
fatgraph* read_fatgraph_from_file(char* filename, double screen_width, double screen_height);
int write_fatgraph_to_file(fatgraph* fg, char* filename);
int draw_fatgraph_to_file(fatgraph* fg, char* filename, double width, double height);

#endif
