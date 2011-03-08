#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "io.h"



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

  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  




