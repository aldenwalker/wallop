#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fatgraph.h"

/******************************************************************************/
/* fatgraph functions                                                         */
/******************************************************************************/
void fatgraph_print(fatgraph* f) {
  int i,j;
  for (i=0; i<f->num_verts; i++) {
    printf("vert %d: %s\n", i, f->verts[i].name);
    printf("location: %f, %f\n", f->verts[i].loc.x, f->verts[i].loc.y);
    printf("edges:\n");
    for (j=0; j<f->verts[i].num_edges; j++) {
      printf("%d ", f->verts[i].edges[j]);
    } printf("\n");
    for (j=0; j<f->verts[i].num_edges; j++) {
      printf("%d ", f->verts[i].edges_initial[j]);
    } printf("\n");    
    for (j=0; j<f->verts[i].num_edges; j++) {
      printf("(%f,%f) ", f->verts[i].bezier[j].x, f->verts[i].bezier[j].y);
    } printf("\n");
  }
  printf("edges:\n");
  for (i=0; i<f->num_edges; i++) {
    printf("%d: %s %s %s %d %d\n", i, f->edges[i].name, 
                                   f->edges[i].label_forward,
                                   f->edges[i].label_backward, 
                                   f->edges[i].start,
                                   f->edges[i].end);
  }
} 


    
void fatgraph_free(fatgraph* f) {
  int i;
  for (i=0; i<f->num_verts; i++) {
    free(f->verts[i].edges);
    free(f->verts[i].edges_initial);
    free(f->verts[i].bezier);
    free(f->verts[i].name);
  }
  for (i=0; i<f->num_edges; i++) {
    free(f->edges[i].name);
    free(f->edges[i].label_forward);
    free(f->edges[i].label_backward);
  }
  free(f->verts);
  free(f->edges);
}

void fatgraph_copy(fatgraph* F, fatgraph* G) {
  fatgraph_free(F);
  int i,j;
  F->num_verts = G->num_verts;
  F->num_edges = G->num_edges;
  F->verts = (vert*)malloc((F->num_verts)*sizeof(vert));
  F->edges = (edge*)malloc((F->num_edges)*sizeof(edge));
  for (i=0; i<G->num_verts; i++) {
    F->verts[i].num_edges = G->verts[i].num_edges;
    F->verts[i].name = (char*)malloc((strlen(G->verts[i].name)+1)*sizeof(char));
    strcpy(F->verts[i].name, G->verts[i].name);
    F->verts[i].edges = (int*)malloc((F->verts[i].num_edges)*sizeof(int));
    F->verts[i].bezier = (vector2d*)malloc((F->verts[i].num_edges)*sizeof(vector2d));
    F->verts[i].edges_initial = (int*)malloc((F->verts[i].num_edges)*sizeof(int));
    F->verts[i].loc.x = G->verts[i].loc.x;
    F->verts[i].loc.y = G->verts[i].loc.y;
    for (j=0; j<F->verts[i].num_edges; j++) {
      F->verts[i].edges[j] = G->verts[i].edges[j];
      F->verts[i].edges_initial[j] = G->verts[i].edges_initial[j];
      F->verts[i].bezier[j] = G->verts[i].bezier[j];
    }
  }
  for (i=0; i<G->num_edges; i++) {
    F->edges[i].name = (char*)malloc((strlen(G->edges[i].name)+1)*sizeof(char));
    F->edges[i].label_forward = (char*)malloc((strlen(G->edges[i].label_forward)+1)*sizeof(char));
    F->edges[i].label_backward = (char*)malloc((strlen(G->edges[i].label_backward)+1)*sizeof(char));
    strcpy(F->edges[i].name, G->edges[i].name);
    strcpy(F->edges[i].label_forward, G->edges[i].label_forward);
    strcpy(F->edges[i].label_backward, G->edges[i].label_backward);
    F->edges[i].start = G->edges[i].start;
    F->edges[i].end = G->edges[i].end;
  }
}


void fatgraph_add_vertex(fatgraph* fg, double x, double y) {
  fg->verts = (vert*)realloc((void*)(fg->verts), (fg->num_verts+1)*sizeof(vert));
  fg->verts[fg->num_verts].name = (char*)malloc(10*sizeof(char));
  sprintf(fg->verts[fg->num_verts].name, "VERT%d", fg->num_verts);
  fg->verts[fg->num_verts].loc.x = x;
  fg->verts[fg->num_verts].loc.y = y;
  fg->verts[fg->num_verts].num_edges = 0;
  fg->verts[fg->num_verts].edges = NULL;
  fg->verts[fg->num_verts].edges_initial = NULL;
  fg->verts[fg->num_verts].bezier = NULL;
  fg->num_verts ++;
}

void fatgraph_add_edge(fatgraph* fg, 
                       int s, 
                       int index_in_s, 
                       vector2d sbezier,
                       int e,
                       int index_in_e,
                       vector2d ebezier){
  //printf("adding edge from %d to %d with bezier %f,%f and %f,%f\n", s,e,sbezier.x,sbezier.y,ebezier.x,ebezier.y);
  fg->edges = (edge*)realloc((void*)(fg->edges), (fg->num_edges+1)*sizeof(edge));
  fg->edges[fg->num_edges].name = (char*)malloc(10*sizeof(char));
  sprintf(fg->edges[fg->num_edges].name, "EDGE%d", fg->num_edges);
  fg->edges[fg->num_edges].start = s;
  fg->edges[fg->num_edges].end = e;
  fg->edges[fg->num_edges].label_forward = (char*)malloc(2*sizeof(char));
  strcpy(fg->edges[fg->num_edges].label_forward, "x");
  fg->edges[fg->num_edges].label_backward = (char*)malloc(2*sizeof(char));
  strcpy(fg->edges[fg->num_edges].label_backward, "X");
  
  
  insert_space_int(&(fg->verts[s].edges), fg->verts[s].num_edges, index_in_s);
  insert_space_int(&(fg->verts[s].edges_initial), fg->verts[s].num_edges, index_in_s);
  insert_space_vector2d(&(fg->verts[s].bezier), fg->verts[s].num_edges, index_in_s);
  fg->verts[s].num_edges ++;
  fg->verts[s].edges[index_in_s] = fg->num_edges;
  fg->verts[s].edges_initial[index_in_s] = 1;
  fg->verts[s].bezier[index_in_s] = sbezier;
  
  insert_space_int(&(fg->verts[e].edges), fg->verts[e].num_edges, index_in_e);
  insert_space_int(&(fg->verts[e].edges_initial), fg->verts[e].num_edges, index_in_e);
  insert_space_vector2d(&(fg->verts[e].bezier), fg->verts[e].num_edges, index_in_e);
  fg->verts[e].num_edges ++;
  fg->verts[e].edges[index_in_e] = fg->num_edges;
  fg->verts[e].edges_initial[index_in_e] = 0;
  fg->verts[e].bezier[index_in_e] = ebezier;

  fg->num_edges ++;
}

fatgraph* fatgraph_create() {
  fatgraph* fg = (fatgraph*)malloc(sizeof(fatgraph));
  fg->verts = NULL;
  fg->edges = NULL;
  fg->num_verts = 0;
  fg->num_edges = 0;
  return fg;
}
  


int fatgraph_close_to_edge(fatgraph* fg, double x, double y) {
  int i,j;
  double closest = 1e10;
  int closestEdge = -1;
  double dist;
  vector2d bezier1, bezier2;
  vector2d control1, control2;
  for (i=0; i<fg->num_edges; i++) {
    for (j=0; j<fg->verts[fg->edges[i].start].num_edges; j++) {
      if (i == fg->verts[fg->edges[i].start].edges[j] 
          && fg->verts[fg->edges[i].start].edges_initial[j] == 1) {
        bezier1 = fg->verts[fg->edges[i].start].bezier[j];
        break;
      }
    }
    for (j=0; j<fg->verts[fg->edges[i].end].num_edges; j++) {
      if (i == fg->verts[fg->edges[i].end].edges[j]
          && fg->verts[fg->edges[i].end].edges_initial[j] == 0) {
        bezier2 = fg->verts[fg->edges[i].end].bezier[j];
        break;
      }
    }
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2);
    //printf("I'm checking edge %d between %d and %d\n", i, fg->edges[i].start, fg->edges[i].end);
    dist = distance_to_interior_bezier(fg->verts[fg->edges[i].start].loc,
                                        fg->verts[fg->edges[i].end].loc,
                                        control1,
                                        control2,
                                        x,y);
    //printf("distance to bezier: %f\n", dist);
    if (dist < closest) {
      closestEdge = i;
      closest = dist;
    }
  }
  if (closest < 100) {
    return closestEdge;
  } else {
    return -1;
  }
}


//returns the distance along the edge which is farthest from the other edges
double fatgraph_remotest(fatgraph* fg, int edge) {
  //go through 100 points on edge and find the one which is farthest from 
  //all the other edges
  vector2d currentPoint;
  double currentDist;
  double minDistToThisPoint;
  double pointMaxDistAway;
  double maxDist;
  double i;
  int j;
  int k;
  vector2d* bezier1 = (vector2d*)malloc((fg->num_edges)*sizeof(vector2d));
  vector2d* bezier2 = (vector2d*)malloc((fg->num_edges)*sizeof(vector2d));
  
  for (j=0; j<fg->num_edges; j++) {
    //find the beziers
    for (k=0; k<fg->verts[fg->edges[j].start].num_edges; k++) {
      if (j == fg->verts[fg->edges[j].start].edges[k] 
            && fg->verts[fg->edges[j].start].edges_initial[k] == 1) {
          bezier1[j] = fg->verts[fg->edges[j].start].bezier[k];
          break;
      }
    }
    for (k=0; k<fg->verts[fg->edges[j].end].num_edges; k++) {
      if (j == fg->verts[fg->edges[j].end].edges[k]
          && fg->verts[fg->edges[j].end].edges_initial[k] == 0) {
          bezier2[j] = fg->verts[fg->edges[j].end].bezier[k];
          break;
      }
    }
    bezier1[j] = vector2d_add(fg->verts[fg->edges[j].start].loc, bezier1[j]);
    bezier2[j] = vector2d_add(fg->verts[fg->edges[j].end].loc, bezier2[j]);  
  }
  
  //printf("I'm finding the remotest point on edge %d\n", edge);
  
  //now find the best point
  maxDist = 0;
  
  for (i=0; i<1; i+=0.01) {
    minDistToThisPoint = 1e8;
    currentPoint = point_on_bezier(fg->verts[fg->edges[edge].start].loc,
                                   fg->verts[fg->edges[edge].end].loc,
                                   bezier1[edge],
                                   bezier2[edge],
                                   i);
    //distances to edges
    for (j=0; j<fg->num_edges; j++) {
      if (j==edge) {
        continue;
      }
      currentDist = distance_to_bezier_fast(fg->verts[fg->edges[j].start].loc,
                                            fg->verts[fg->edges[j].end].loc,
                                            bezier1[j],
                                            bezier2[j],
                                            currentPoint.x,
                                            currentPoint.y);
      if (currentDist < minDistToThisPoint) {
        minDistToThisPoint = currentDist;
      }
    }
    //distances to vertices
    for (j=0; j<fg->num_verts; j++) {
      currentDist = sqrt( (currentPoint.x - fg->verts[j].loc.x)*
                          (currentPoint.x - fg->verts[j].loc.x) +
                          (currentPoint.y - fg->verts[j].loc.y)*
                          (currentPoint.y - fg->verts[j].loc.y) );
      if (currentDist < minDistToThisPoint) {
        minDistToThisPoint = currentDist;
      }     
    }
    if (minDistToThisPoint > maxDist) {
      maxDist = minDistToThisPoint;
      pointMaxDistAway = i;
      //printf("New remotest point (distance %f) at %f\n", maxDist, i);
    }
  }    
  free(bezier1);
  free(bezier2);         
  return pointMaxDistAway;
}


/********/
/* fatgraph optimization functions are in op.h */
/********/

int fatgraph_chi(fatgraph* fg) {
  return fg->num_verts - fg->num_edges;
}


void fatgraph_next_edge(fatgraph* fg, int* edge, int* dir) {
  int nextVert;
  int indInVert;
  int nextEdgeInVert;
  if (*dir == 1) {
    nextVert = fg->edges[*edge].end;
    for (indInVert = 0; indInVert<fg->verts[nextVert].num_edges; indInVert++) {
      if (fg->verts[nextVert].edges[indInVert] == *edge  &&
          fg->verts[nextVert].edges_initial[indInVert] == 0) {
        break;
      }
    }
    if (indInVert == fg->verts[nextVert].num_edges) {
      printf("Couldn't find an edge?\n");
      return;
    }
    nextEdgeInVert =  (indInVert+1)%(fg->verts[nextVert].num_edges);
    //printf("The next edge after %d in direction %d is", *edge, *dir);
    *edge = fg->verts[nextVert].edges[nextEdgeInVert];
    *dir = (fg->verts[nextVert].edges_initial[nextEdgeInVert] == 1 ? 1 : -1);
    //printf(" %d in direction %d\n", *edge, *dir);
  } else { 
    nextVert = fg->edges[*edge].start;
    for (indInVert = 0; indInVert<fg->verts[nextVert].num_edges; indInVert++) {
      if (fg->verts[nextVert].edges[indInVert] == *edge  &&
          fg->verts[nextVert].edges_initial[indInVert] == 1) {
        break;
      }
    }
    if (indInVert == fg->verts[nextVert].num_edges) {
      printf("Couldn't find an edge?\n");
      return;
    }
    nextEdgeInVert =  (indInVert+1)%(fg->verts[nextVert].num_edges);
    //printf("The next edge after %d in direction %d is", *edge, *dir);
    *edge = fg->verts[nextVert].edges[nextEdgeInVert];
    *dir = (fg->verts[nextVert].edges_initial[nextEdgeInVert] == 1 ? 1 : -1);
    //printf(" %d in direction %d\n", *edge, *dir);
  }
}

char* fatgraph_follow_edge(fatgraph* fg, int dir, int edge, int* edges_remaining) {
  char* curBound = NULL;
  int current_edge = edge;
  int current_dir = dir;   //1 = forwards, -1 = backwards
  //printf("Following edge: %d in direction %d\n", current_edge, current_dir);
  if (dir == 1) {
    curBound = (char*)malloc((strlen(fg->edges[edge].label_forward)+1)*sizeof(char));
    strcpy(curBound, fg->edges[edge].label_forward);
    edges_remaining[edge] ^= 1;
  } else {
    curBound = (char*)malloc((strlen(fg->edges[edge].label_backward)+1)*sizeof(char));
    strcpy(curBound, fg->edges[edge].label_backward);
    edges_remaining[edge] ^= 2;
  }
  fatgraph_next_edge(fg, &current_edge, &current_dir);
  while (current_edge != edge || current_dir != dir) {
    //printf("Following edge: %d in direction %d\n", current_edge, current_dir);
    if (current_dir == 1) {
      curBound = (char*)realloc((void*)curBound,
                                (strlen(curBound) + strlen(fg->edges[current_edge].label_forward) + 1)*sizeof(char));
      strcat(curBound, fg->edges[current_edge].label_forward);
      edges_remaining[current_edge] ^= 1;
    } else {
      curBound = (char*)realloc((void*)curBound,
                                (strlen(curBound) + strlen(fg->edges[current_edge].label_backward) + 1)*sizeof(char));
      strcat(curBound, fg->edges[current_edge].label_backward);
      edges_remaining[current_edge] ^= 2;
    }
    fatgraph_next_edge(fg, &current_edge, &current_dir);
  }
  //printf("Broke with boundary: %s\n", curBound);
  return curBound; 
}




char** fatgraph_boundaries(fatgraph* fg, int* numBoundaries) {
  char** boundaries = NULL;
  *numBoundaries = 0;
  int* edges_remaining = (int*)malloc((fg->num_edges)*sizeof(int));
  int i;//,j;
  for (i=0; i<fg->num_edges; i++) {
    edges_remaining[i] = 3; //i.e. 11  (edges_remaining[i] & 1) is forward edge, (edges_remaining[i] >> 1)&1 is backward
  }
  while (1) {
    for (i=0; i<fg->num_edges; i++) {
      if (edges_remaining[i] != 0) {
        break;
      }
    }
    if (i == fg->num_edges) {
      break; //we're done
    }
    //printf("Trying to find another edge to start on; edges_remaining:\n");
    //for (j=0; j<fg->num_edges; j++) {
    //  printf("%d ", edges_remaining[j]);
    //} printf("\n");
    //printf("Trying forward edge %d\n", i);
    if ((edges_remaining[i]&1) == 1) { //if we need to do forward
      //printf("Doing forward\n");
      //follow the forward edge
      boundaries = (char**)realloc((void*)boundaries, (*numBoundaries+1)*sizeof(char*));
      boundaries[*numBoundaries] = fatgraph_follow_edge(fg, 1, i, edges_remaining);
      (*numBoundaries)++;
    }
    //printf("Trying backward edge %d\n", i);
    if (((edges_remaining[i]>>1)&1) == 1) { //if we need to do backward
      //printf("Doing backward\n");
      //follow the forward edge
      boundaries = (char**)realloc((void*)boundaries, (*numBoundaries+1)*sizeof(char*));
      boundaries[*numBoundaries] = fatgraph_follow_edge(fg, -1, i, edges_remaining);
      (*numBoundaries)++;
    }
  }
  return boundaries;
}




int write_fatgraph_to_file(fatgraph* fg, char* filename) {
  int i,j;
  FILE* ofile = fopen(filename, "w");
  if (ofile == NULL) {
    return 1;
  }
  fprintf(ofile, "#output produced from fatgraph program\n");
  fprintf(ofile, "vertices %d\n", fg->num_verts);
  for (i=0; i<fg->num_verts; i++) {
    if (fg->verts[i].name == NULL) {
      fprintf(ofile, "VERT%d %d\n", i, fg->verts[i].num_edges);
    } else {
      fprintf(ofile, "%s %d\n", fg->verts[i].name, fg->verts[i].num_edges);
    }
    for (j=0; j<fg->verts[i].num_edges; j++) {
      if (fg->edges[fg->verts[i].edges[j]].name == NULL) {
        fprintf(ofile, "EDGE%d ", fg->verts[i].edges[j]);
      } else {
        fprintf(ofile, "%s ", fg->edges[fg->verts[i].edges[j]].name);
      }
    }
    fprintf(ofile, "\n");
    for (j=0; j<fg->verts[i].num_edges; j++) {
      fprintf(ofile, "%d ", fg->verts[i].edges_initial[j]);
    }
    if (fg->verts[i].num_edges > 0) {
      fprintf(ofile, "\nbezier ");
      for (j=0; j<fg->verts[i].num_edges; j++) {
        fprintf(ofile, "%f %f ", fg->verts[i].bezier[j].x, fg->verts[i].bezier[j].y);
      }
    }
    fprintf(ofile, "\nloc ");
    fprintf(ofile, "%f %f ", fg->verts[i].loc.x, fg->verts[i].loc.y);
    fprintf(ofile, "\n");
  }
  fprintf(ofile, "edges %d\n", fg->num_edges);
  for (i=0; i<fg->num_edges; i++) {
    if (fg->edges[i].name == NULL) {
      fprintf(ofile, "EDGE%d ", i);
    } else {
      fprintf(ofile, "%s ", fg->edges[i].name);
    }
    fprintf(ofile, "%s %s %s %s\n", fg->edges[i].label_forward, 
                                    fg->edges[i].label_backward,
                                    fg->verts[fg->edges[i].start].name,
                                    fg->verts[fg->edges[i].end].name);
  }
  fclose(ofile);
  return 0;
}
  

fatgraph* read_fatgraph_from_file(char* filename, double screen_width,
                                                  double screen_height) {
  fatgraph* fg = NULL;
  FILE* ifile = fopen(filename, "r");
  if (ifile == NULL) {
    return NULL;
  }
  
  fg = (fatgraph*)malloc(sizeof(fatgraph));
  double th;
  int i,j,k;
  char tempc;
  char temps[1000];
  char temp_edge_name[50];
  char temp_label_forward[50];
  char temp_label_backward[50];
  char temp_vert_start[50];
  char temp_vert_end[50];
  
  
  while (fgetc(ifile) == '#') {  
    do {
      tempc = fgetc(ifile);
    } while (tempc != EOF && tempc != '\n');
  }
  fseek(ifile, -1, SEEK_CUR);
  
  
  fscanf(ifile, "vertices %d\n", &(fg->num_verts));
  fg->verts = (vert*)malloc((fg->num_verts)*sizeof(vert));
  for (i=0; i<fg->num_verts; i++) {
    fscanf(ifile, "%s %d\n", temps, &(fg->verts[i].num_edges));
    printf("new vertex %s with %d edges\n", temps, fg->verts[i].num_edges); 
    fg->verts[i].name = (char*)malloc((strlen(temps)+1)*sizeof(char));
    strcpy(fg->verts[i].name, temps);  
    if (fg->verts[i].num_edges > 0 ) {
      fg->verts[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
      fg->verts[i].edges_initial = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
      fg->verts[i].bezier = (vector2d*)malloc((fg->verts[i].num_edges)*sizeof(vector2d));
    } else {
      fg->verts[i].edges = NULL;
      fg->verts[i].edges_initial = NULL;
      fg->verts[i].bezier = NULL;
    }
    
    if (fg->verts[i].num_edges > 0) {
      //skip the edge line for now
      fgets(temps, 1000, ifile);
      
      //read in the initial edges (even though we don't know what they are!)
      for (j=0; j<fg->verts[i].num_edges; j++) {
        fscanf(ifile, "%d", &(fg->verts[i].edges_initial[j]));
      }
    }
    
    //read the bezier and location, if we've got them
    fscanf(ifile, " %s", temps);
    printf("next word: %s\n", temps);
    if (strcmp(temps, "bezier") == 0) {
      printf("It's got bezier stuff\n");
      for (j=0; j<fg->verts[i].num_edges; j++) {
        fscanf(ifile, "%lf %lf", &(fg->verts[i].bezier[j].x),&(fg->verts[i].bezier[j].y));
        printf("read bezier: %lf, %lf\n", fg->verts[i].bezier[j].x, fg->verts[i].bezier[j].y);
      }
      fscanf(ifile, "%s", temps);
    } else {
      //there's no bezier, so deal with that
      printf("There's no bezier\n");
      if (fg->verts[i].num_edges > 0) {
        printf("But there are %d edges, so I'm adding that in\n", fg->verts[i].num_edges);
        th = 0;
        for (j=0; j<fg->verts[i].num_edges; j++) {
          th += 6.4/(double)fg->verts[i].num_edges;
          fg->verts[i].bezier[j].x = cos(th);
          fg->verts[i].bezier[j].y = sin(th);
        }
      }

    }
    if (strcmp(temps, "loc") == 0) {
      fscanf(ifile, "%lf %lf", &(fg->verts[i].loc.x),&(fg->verts[i].loc.y));
      printf("read loc: %lf, %lf\n", fg->verts[i].loc.x, fg->verts[i].loc.y);
    } else {
      //there's no loc, so deal with that
      fg->verts[i].loc.x = (rand()/(double)RAND_MAX)*screen_width;
      fg->verts[i].loc.y = (rand()/(double)RAND_MAX)*screen_height;
      fseek(ifile, -strlen(temps)-1, SEEK_CUR);
    } 
  }
  
  printf("I've read in the vertices\n"); fflush(stdout);
  //fgets(temps, 1000, ifile);
  //fgets(temps, 1000, ifile);
  //printf("The next thing is: %s\n", temps);
  //fseek(ifile, -strlen(temps)-1, SEEK_CUR);
  
  
  fscanf(ifile, " edges %d\n", &(fg->num_edges));
  printf("There should be %d edges\n", fg->num_edges); fflush(stdout);
  fg->edges = (edge*)malloc((fg->num_edges)*sizeof(edge));
  for (i=0; i<fg->num_edges; i++) {
    fscanf(ifile, " %s %s %s %s %s\n", temp_edge_name, 
                                      temp_label_forward,
                                      temp_label_backward,
                                      temp_vert_start,
                                      temp_vert_end);
    printf("Read edge: %s, %s, %s, %s, %s\n", temp_edge_name, 
                                      temp_label_forward,
                                      temp_label_backward,
                                      temp_vert_start,
                                      temp_vert_end);
    fg->edges[i].name = (char*)malloc((strlen(temp_edge_name)+1)*sizeof(char));
    fg->edges[i].label_forward = (char*)malloc((strlen(temp_label_forward)+1)*sizeof(char));
    fg->edges[i].label_backward = (char*)malloc((strlen(temp_label_backward)+1)*sizeof(char));
    strcpy(fg->edges[i].name, temp_edge_name);
    strcpy(fg->edges[i].label_forward, temp_label_forward);
    strcpy(fg->edges[i].label_backward, temp_label_backward);
    for (j=0; j<fg->num_verts; j++) {
      if (strcmp(temp_vert_start, fg->verts[j].name)==0) {
        fg->edges[i].start = j;
        break;
      }
    }
    for (j=0; j<fg->num_verts; j++) {
      if (strcmp(temp_vert_end, fg->verts[j].name)==0) {
        fg->edges[i].end = j;
        break;
      }
    }
  }
  
  //reset to the beginning
  fseek(ifile, 0, SEEK_SET);
  
  if (fgetc(ifile) == '#') {
    do {
      tempc = fgetc(ifile);
    } while (tempc != EOF && tempc != '\n');
  }
  
  fscanf(ifile, "vertices %*d\n");
  for (i=0; i<fg->num_verts; i++) {
    fscanf(ifile, " %s ", temps);
    while (strcmp(temps, fg->verts[i].name) != 0) {
      fscanf(ifile, " %s ", temps);
    }
    fscanf(ifile, "%*d"); //get rid of the number of edges
    for (j=0; j<fg->verts[i].num_edges; j++) {
      fscanf(ifile, " %s ", temp_edge_name);
      for (k=0; k<fg->num_edges; k++) {
        if (strcmp(temp_edge_name, fg->edges[k].name)==0) {
          break;
        }
      }
      fg->verts[i].edges[j] = k;
    } 
  }
  
  return fg;
} 


/*****************************************************************************/
/* draw a fatgraph as an eps file                                            */
/*****************************************************************************/
int draw_fatgraph_to_file(fatgraph* fg, char* filename, double width, double height) {
  int i,j;
  double minx, miny, maxx, maxy, x, y;
  vector2d bezier1, bezier2, control1, control2;
  
  FILE* ofile = fopen(filename, "w");
  if (ofile == NULL) {
    return 1;
  }
  //////////////////////////
  //Fast drawing (just the vertices and edges)
  ///////////////////////////
  if (fg->num_verts == 0) {
    fprintf(ofile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(ofile, "%%%%BoundingBox: 0 0 %f %f\n\n", width, height);      
    fprintf(ofile, "1 setlinejoin\n");
    fprintf(ofile, "2 setlinewidth\n");    
    fclose(ofile);
    return 0;
  }
  /*
  //find a good bounding box
  minx = maxx = fg->verts[0].loc.x;
  miny = maxy = fg->verts[0].loc.y;
  for (i=0; i<fg->num_verts; i++) {
    x = fg->verts[i].loc.x;
    y = fg->verts[i].loc.y;
    if (x < minx) {
      minx = x;
    } else if (x > maxx) {
      maxx = x;
    }
    if (y < miny) {
      miny = y;
    } else if (y > maxy) {
      maxy = y;
    }
  }
  */
    
  fprintf(ofile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(ofile, "%%%%BoundingBox: %f %f %f %f\n\n", 0.0, 0.0, width, height);      
  fprintf(ofile, "1 setlinejoin\n");
  fprintf(ofile, "5 setlinewidth\n");
  fprintf(ofile, "0 0 0 setrgbcolor\n");
  for (i=0; i<fg->num_edges; i++) {
    for (j=0; j<fg->verts[fg->edges[i].start].num_edges; j++) {
      //printf("Looking at edge %d from vertex %d--it's edge %d\n", j, fg->edges[i].start, 
      //                                                               fg->verts[fg->edges[i].start].edges[j] );
      //printf("Vertex %d says that edge %d (%d) is initial: %d here\n", fg->edges[i].start,
      //                                                                 j,
      //                                                                 fg->verts[fg->edges[i].start].edges[j],
      //                                                                 fg->verts[fg->edges[i].start].edges_initial[j]);
      if (i == fg->verts[fg->edges[i].start].edges[j] 
          && fg->verts[fg->edges[i].start].edges_initial[j] == 1) {
        bezier1 = fg->verts[fg->edges[i].start].bezier[j];
        break;
      }
    }
    for (j=0; j<fg->verts[fg->edges[i].end].num_edges; j++) {
      //printf("Looking at edge %d from vertex %d--it's edge %d\n", j, fg->edges[i].end, 
      //                                                               fg->verts[fg->edges[i].end].edges[j] );
      //printf("Vertex %d says that edge %d (%d) is initial: %d here\n", fg->edges[i].end,
      //                                                                 j,
      //                                                                 fg->verts[fg->edges[i].end].edges[j],
      //                                                                 fg->verts[fg->edges[i].end].edges_initial[j]);
      if (i == fg->verts[fg->edges[i].end].edges[j]
          && fg->verts[fg->edges[i].end].edges_initial[j] == 0) {
        bezier2 = fg->verts[fg->edges[i].end].bezier[j];
        break;
      }
    }
    //set the bezier controls
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2);
    fprintf(ofile, "1 1 1 setrgbcolor\n");
    fprintf(ofile, "16 setlinewidth\n");
    fprintf(ofile, "%f %f moveto\n", fg->verts[fg->edges[i].start].loc.x, 
                                     fg->verts[fg->edges[i].start].loc.y);
    fprintf(ofile, "%f %f %f %f %f %f curveto\n", control1.x, control1.y,
                                                  control2.x, control2.y,
                                                  fg->verts[fg->edges[i].end].loc.x,
                                                  fg->verts[fg->edges[i].end].loc.y);
    fprintf(ofile, "stroke\n");
    fprintf(ofile, "0 0 0 setrgbcolor\n");
    fprintf(ofile, "8 setlinewidth\n");
    fprintf(ofile, "%f %f moveto\n", fg->verts[fg->edges[i].start].loc.x, 
                                     fg->verts[fg->edges[i].start].loc.y);
    fprintf(ofile, "%f %f %f %f %f %f curveto\n", control1.x, control1.y,
                                                  control2.x, control2.y,
                                                  fg->verts[fg->edges[i].end].loc.x,
                                                  fg->verts[fg->edges[i].end].loc.y);
    fprintf(ofile, "stroke\n");
  }
  
  //draw the vertices
  fprintf(ofile, "11 setlinewidth\n");
  for (i=0; i<fg->num_verts; i++) {
    fprintf(ofile, "%f %f moveto\n", fg->verts[i].loc.x, fg->verts[i].loc.y);
    fprintf(ofile, "%f %f %f %f %f arc\n", fg->verts[i].loc.x, 
                                    fg->verts[i].loc.y,
                                    4.0,
                                    0.0, 360.0);
    fprintf(ofile, "stroke\n");
  }
  fclose(ofile);
  return 0;
  
}






