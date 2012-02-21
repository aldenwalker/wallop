#include <gtk/gtk.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vector2d.h"
#include "fatgraph.h"
#include "op.h"
//#include "io.h"
#include "wallop.h"

 
#define SQUARE(x) (x)*(x)


/******************************************************************************/
/* globals                                                                    */
/******************************************************************************/
static GdkPixmap* pixmap = NULL;
static cairo_t* cr = NULL;
static fatgraph* fg = NULL;
static struct { 
          enum { FG_NO_ACTION, 
                 FG_DRAWING_EDGE1,
                 FG_DRAWING_EDGE2,
                 FG_DRAGGING_VERT } action;
          enum { FG_NO_FLAG,
                 FG_DONE_EDGE1,
                 FG_EDGE_HIGHLIGHT } flag;
          double edge1_down_x; 
          double edge1_down_y;
          double edge1_up_x;
          double edge1_up_y;
          double edge2_down_x;
          double edge2_down_y;
          double edge2_up_x;
          double edge2_up_y;
          int dragging_vert; 
          int highlightedEdge;
       } dstatus;

typedef struct {
  GtkEntry* entry1;
  GtkEntry* entry2;
  GtkWidget* widget;
} entryWidgetPair;


typedef struct {
  GtkEntry* entry1;
  GtkEntryBuffer* text1;
  GtkEntry* entry2;
  GtkEntryBuffer* text2;
  GtkTextView* textView;
} displayBoxes;

typedef struct {
  GtkEntry* filenameEntry;
  GtkEntry* forwardLabelEntry;
  GtkEntry* backwardLabelEntry;
  GtkEntry* scallopEntry;
  GtkTextView* dataTextView;
  GtkWidget* widget;
} widgetList;
  


/******************************************************************************/
/* various handy functions                                                    */
/******************************************************************************/
int gcd(int a, int b) {
  int t1 = (b > a ? b : a); //max
  int t2 = (b > a ? a : b); //min
  int t3;
  if (a==0 || b==0) {
    return 1;
  }
  while (t2 != 0) {
    t3 = t1 % t2;
    t1 = t2;
    t2 = t3;
  }
  return t1;
}


int lcm(int a, int b) {
  int g = gcd(a,b);
  return (a*b)/g;
}

int list_gcd(int* L, int len) {
  int g;
  int i;
  if (len == 0) {
    return 1;
  } else if (len == 1) {
    return L[0];
  } else if (len == 2) {
    return gcd(L[0],L[1]);
  } else {
    g = gcd(L[0], L[1]);
    for (i=2; i<len; i++) {
      g = gcd(g, L[i]);
    }
    return g;
  }
} 



void insert_space_int(int** L, int len, int index_to_clear) {
  int i;
  (*L) = (int*)realloc((void*)(*L), (len+1)*sizeof(int));
  for (i=len; i>index_to_clear; i--) {
    (*L)[i] = (*L)[i-1];
  }
}
void insert_space_vector2d(vector2d** L, int len, int index_to_clear) {
  int i;
  (*L) = (vector2d*)realloc((void*)(*L), (len+1)*sizeof(vector2d));
  for (i=len; i>index_to_clear; i--) {
    (*L)[i] = (*L)[i-1];
  }
}

int in_circular_order(double a, double b, double c) {
  if ( (a < b && b < c) || (b < c && c < a) || (c < a && a < b) ) {
    return 1;
  } else {
    return 0;
  }
}




  











void fg_cairo_draw_arrow(cairo_t* ca, vector2d loc, vector2d tangent, double size) {
  vector2d normal;
  normal.x = tangent.y;
  normal.y = -tangent.x;  
  cairo_move_to(ca, loc.x - size*0.5*tangent.x + size*0.5*normal.x, 
                    loc.y - size*0.5*tangent.y + size*0.5*normal.y);
  cairo_line_to(ca, loc.x, loc.y);
  cairo_line_to(ca, loc.x - size*0.5*tangent.x - size*0.5*normal.x, 
                    loc.y - size*0.5*tangent.y - size*0.5*normal.y);
  cairo_stroke(ca);
}
  
void show_text_flipped(cairo_t* ca, double x, double y, char* text) {
  double dx, dy;
  dx = x;
  dy = y;
  cairo_user_to_device(ca, &dx, &dy);
  cairo_save(ca);
  cairo_identity_matrix(ca);
  cairo_move_to(ca, dx, dy);
  cairo_show_text(ca, text);
  cairo_restore(ca);
}
  







void fatgraph_draw(fatgraph* fg, 
                   cairo_t* ca, 
                   int width, 
                   int height, 
                   int redraw,
                   int fast) {
  int i,j;
  
  double outside_rectangle_width=20;
  double rectangle_width = 15.0;
  double inside_rectangle_width = 11.0;
  
  vector2d* bezier1 = (vector2d*)malloc((fg->num_edges)*sizeof(vector2d));
  vector2d* bezier2 = (vector2d*)malloc((fg->num_edges)*sizeof(vector2d));
  
  vector2d control1, control2;
  vector2d midpoint, tangent, normal, font_offset, pt_to_center;
  vector2d trunc_start, trunc_end;
  vector2d trunc_start_bezier, trunc_end_bezier;
  vector2d arrowCenter;
  cairo_text_extents_t extents;
  double normalNorm, tangentNorm;
  double font_diagonal;
  double* edge_arrow_position = (double*)malloc((fg->num_edges)*sizeof(double));

  cairo_set_source_rgb(ca, 1,1,1);
  cairo_set_line_width(ca, 0.8*rectangle_width);
  cairo_select_font_face(ca, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(ca, 1*rectangle_width);
  
  if (redraw==1) {
    cairo_rectangle(ca, 0,0, width,height);
    cairo_fill(ca);
  }
  
  //////////////////////////
  //Fast drawing (just the vertices and edges)
  ///////////////////////////
  if (fast == 1) {
    cairo_set_source_rgb(ca, 0,0,0);
    cairo_set_line_width(ca, rectangle_width/3);
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
          bezier1[i] = fg->verts[fg->edges[i].start].bezier[j];
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
          bezier2[i] = fg->verts[fg->edges[i].end].bezier[j];
          break;
        }
      }
      //set the bezier controls
      control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1[i]);
      control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2[i]);
      cairo_move_to(ca, fg->verts[fg->edges[i].start].loc.x, 
                        fg->verts[fg->edges[i].start].loc.y);
      cairo_curve_to(ca, control1.x,
                         control1.y,
                         control2.x,
                         control2.y,
                         fg->verts[fg->edges[i].end].loc.x,
                         fg->verts[fg->edges[i].end].loc.y);
      cairo_stroke(ca);
    }
    
    //draw the vertices
    cairo_set_line_width(ca, inside_rectangle_width/4);
    for (i=0; i<fg->num_verts; i++) {
      cairo_move_to(ca, fg->verts[i].loc.x, fg->verts[i].loc.y);
      cairo_arc(ca, fg->verts[i].loc.x, fg->verts[i].loc.y, inside_rectangle_width/6, 0, 2*3.14159);
      cairo_stroke(ca);
    }
    
    return;
  }
   
  
  
  
  //////////////////////////
  //detailed drawing
  //////////////////////////
  //draw the arcs
  cairo_set_source_rgb(ca, 0,0,0);
  cairo_set_line_width(ca, rectangle_width);
  for (i=0; i<fg->num_edges; i++) {  //this needs to be faster
    //printf("drawing edge %d, firom %d to %d\n", i, fg->edges[i].start, fg->edges[i].end); fflush(stdout);
    for (j=0; j<fg->verts[fg->edges[i].start].num_edges; j++) {
      //printf("Looking at edge %d from vertex %d--it's edge %d\n", j, fg->edges[i].start, 
      //                                                               fg->verts[fg->edges[i].start].edges[j] );
      //printf("Vertex %d says that edge %d (%d) is initial: %d here\n", fg->edges[i].start,
      //                                                                 j,
      //                                                                 fg->verts[fg->edges[i].start].edges[j],
      //                                                                 fg->verts[fg->edges[i].start].edges_initial[j]);
      if (i == fg->verts[fg->edges[i].start].edges[j] 
          && fg->verts[fg->edges[i].start].edges_initial[j] == 1) {
        bezier1[i] = fg->verts[fg->edges[i].start].bezier[j];
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
        bezier2[i] = fg->verts[fg->edges[i].end].bezier[j];
        break;
      }
    }
    //printf("set bezier vectors: (%f, %f) and (%f, %f)\n", bezier1[i].x,
     //                                                     bezier1[i].y,
     //                                                     bezier2[i].x,
    //                                                      bezier2[i].y);
    
    //printf("Set control vectors\n");
  
    //printf("Drawing a bezier with points: %f,%f, %f,%f, %f,%f, %f,%f\n", fg->verts[fg->edges[i].start].loc.x,
    //                                                                     fg->verts[fg->edges[i].start].loc.y,
    //                                                                     control1.x,
    //                                                                     control1.y,
    //                                                                     control2.x, 
    //                                                                     control2.y,
    //                                                                     fg->verts[fg->edges[i].end].loc.x,
     //                                                                    fg->verts[fg->edges[i].end].loc.y);
    
    
    //set the bezier controls
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1[i]);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2[i]);
    
    
    
    //DRAW the middle of all the arcs
    
    //draw the outside (white background) arcs
    sub_bezier(fg->verts[fg->edges[i].start].loc, 
               fg->verts[fg->edges[i].end].loc,
               control1, 
               control2,
               0.05, 0.95,
               &(trunc_start),
               &(trunc_end),
               &(trunc_start_bezier),
               &(trunc_end_bezier));
    cairo_set_line_width(ca, outside_rectangle_width);
    cairo_set_source_rgb(ca, 1,1,1);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    cairo_stroke(ca);
    
    
    //draw the middle arcs
    if (dstatus.flag == FG_EDGE_HIGHLIGHT && dstatus.highlightedEdge == i) {
      cairo_set_source_rgb(ca, 1,0.1,0.1);
    } else {
      cairo_set_source_rgb(ca, 0,0,0);
    }
    cairo_set_line_width(ca, rectangle_width);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    //cairo_line_to(ca, fg->verts[fg->edges[i].end].x, 
    //                 fg->verts[fg->edges[i].end].y);
    cairo_stroke(ca);
    if (dstatus.flag == FG_EDGE_HIGHLIGHT && dstatus.highlightedEdge == i) {
      cairo_set_source_rgb(ca, 0,0,0);
    }
    
  
    //draw the inside arcs
    cairo_set_source_rgb(ca, 1,1,1);
    cairo_set_line_width(ca, inside_rectangle_width);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    //cairo_line_to(ca, fg->verts[fg->edges[i].end].x, 
    //                 fg->verts[fg->edges[i].end].y);
    cairo_stroke(ca);
  
  }
  
  
  //now go back and fill in the ends
  //first, the backgrounds
  for (i=0; i<fg->num_edges; i++) {  //this needs to be faster
    //set the bezier controls
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1[i]);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2[i]);
    //start
    sub_bezier(fg->verts[fg->edges[i].start].loc, 
               fg->verts[fg->edges[i].end].loc,
               control1, 
               control2,
               0, 0.05,
               &(trunc_start),
               &(trunc_end),
               &(trunc_start_bezier),
               &(trunc_end_bezier));
    cairo_set_line_width(ca, rectangle_width);
    cairo_set_source_rgb(ca, 0,0,0);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    cairo_stroke(ca);
    //end
    sub_bezier(fg->verts[fg->edges[i].start].loc, 
               fg->verts[fg->edges[i].end].loc,
               control1, 
               control2,
               0.95, 1,
               &(trunc_start),
               &(trunc_end),
               &(trunc_start_bezier),
               &(trunc_end_bezier));
    cairo_set_line_width(ca, rectangle_width);
    cairo_set_source_rgb(ca, 0,0,0);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    cairo_stroke(ca);
  } 
  
    //now, the white part, the backgrounds
  for (i=0; i<fg->num_edges; i++) {  
    //set the bezier controls
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1[i]);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2[i]);
    //start
    sub_bezier(fg->verts[fg->edges[i].start].loc, 
               fg->verts[fg->edges[i].end].loc,
               control1, 
               control2,
               0, 0.052,
               &(trunc_start),
               &(trunc_end),
               &(trunc_start_bezier),
               &(trunc_end_bezier));
    cairo_set_line_width(ca, inside_rectangle_width);
    cairo_set_source_rgb(ca, 1,1,1);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    cairo_stroke(ca);
    //end
    sub_bezier(fg->verts[fg->edges[i].start].loc, 
               fg->verts[fg->edges[i].end].loc,
               control1, 
               control2,
               0.948, 1,
               &(trunc_start),
               &(trunc_end),
               &(trunc_start_bezier),
               &(trunc_end_bezier));
    cairo_set_line_width(ca, inside_rectangle_width);
    cairo_set_source_rgb(ca, 1,1,1);
    cairo_move_to(ca, trunc_start.x, trunc_start.y);
    cairo_curve_to(ca, trunc_start_bezier.x,
                       trunc_start_bezier.y,
                       trunc_end_bezier.x,
                       trunc_end_bezier.y,
                       trunc_end.x,
                       trunc_end.y);
    cairo_stroke(ca);
  } 
  
  
      
      
      
  //decorate the arcs with letters and arrows
  cairo_set_source_rgb(ca, 0,0,0);
  for (i=0; i<fg->num_edges; i++) {
  
    //find the best positions for the arrows
    edge_arrow_position[i] = fatgraph_remotest(fg, i);  
    
    //letter 1
    control1 = vector2d_add(fg->verts[fg->edges[i].start].loc, bezier1[i]);
    control2 = vector2d_add(fg->verts[fg->edges[i].end].loc, bezier2[i]);
    point_and_tangent_to_bezier(fg->verts[fg->edges[i].start].loc,
                                fg->verts[fg->edges[i].end].loc,
                                control1,
                                control2,
                                edge_arrow_position[i],
                                &midpoint,
                                &tangent);
    normal.x = tangent.y; 
    normal.y = -tangent.x;
    //printf("tangent from %f,%f to %f,%f is: %f,%f\n", fg->verts[fg->edges[i].start].loc.x,
    //                                                  fg->verts[fg->edges[i].start].loc.y,
    //                                                  fg->verts[fg->edges[i].end].loc.x,
    //                                                  fg->verts[fg->edges[i].end].loc.y,
    //                                                  tangent.x,
    //                                                  tangent.y);
    //printf("normal is: %f,%f\n", normal.x, normal.y);
    //printf("0.35 of the way is: %f,%f\n", midpoint.x, midpoint.y);
    normalNorm = sqrt((normal.x * normal.x) + (normal.y * normal.y));
    
    //cairo_move_to(ca, midpoint.x, midpoint.y);
    
    //cairo_arc(ca, midpoint.x, midpoint.y, 3, 0, 2*3.14159);
    //cairo_stroke(ca);
   
    cairo_move_to(ca, midpoint.x, midpoint.y);
    cairo_text_extents(ca, fg->edges[i].label_forward, &extents);
    font_diagonal = sqrt(extents.width*extents.width + 
                         extents.height*extents.height);
                         
    //font_offset is where we want the center of the glyph to be
    font_offset.x = normal.x * (rectangle_width+font_diagonal/6)/normalNorm;
    font_offset.y = normal.y * (rectangle_width+font_diagonal/6)/normalNorm;
    //pt_to_center is a vector recording where the center appears if we
    //were to draw it at midpoint (up and to the right)
    pt_to_center.x = extents.width/2;
    pt_to_center.y = extents.height/2; 
    //printf("midpoint: %f, %f\n", midpoint.x, midpoint.y);
    //printf("extents: x_bearing: %f  y_bearing: %f, width: %f, height: %f\n", extents.x_bearing, extents.y_bearing, extents.width, extents.height);
    //printf("pt_to_center: %f, %f\n", pt_to_center.x, pt_to_center.y);   
    
    cairo_move_to(ca, midpoint.x + font_offset.x - pt_to_center.x, 
                      midpoint.y + font_offset.y - pt_to_center.y);
    //cairo_show_text(ca, fg->edges[i].label_forward);
    show_text_flipped(ca, midpoint.x + font_offset.x - pt_to_center.x, 
                          midpoint.y + font_offset.y - pt_to_center.y,
                          fg->edges[i].label_forward);
    
    
    
    //cairo_move_to(ca, midpoint.x, midpoint.y);
    //cairo_line_to(ca, midpoint.x + font_offset.x, midpoint.y + font_offset.y);
    //cairo_stroke(ca);
    //cairo_move_to(ca, midpoint.x, midpoint.y);
    
    
    //letter 2
    point_and_tangent_to_bezier(fg->verts[fg->edges[i].start].loc,
                                fg->verts[fg->edges[i].end].loc,
                                control1,
                                control2,
                                edge_arrow_position[i],
                                &midpoint,
                                &tangent);
    normal.x = tangent.y; 
    normal.y = -tangent.x;
    normalNorm = sqrt((normal.x * normal.x) + (normal.y * normal.y));   
    cairo_move_to(ca, midpoint.x, midpoint.y);
    cairo_text_extents(ca, fg->edges[i].label_backward, &extents);
    font_diagonal = sqrt(extents.width*extents.width + 
                         extents.height*extents.height);
    font_offset.x = normal.x * (rectangle_width+font_diagonal/6)/normalNorm;
    font_offset.y = normal.y * (rectangle_width+font_diagonal/6)/normalNorm;    
    pt_to_center.x = extents.width/2;
    pt_to_center.y = extents.height/2; 
    cairo_move_to(ca, midpoint.x - font_offset.x - pt_to_center.x, 
                      midpoint.y - font_offset.y - pt_to_center.y);
    //cairo_show_text(ca, fg->edges[i].label_backward);
    show_text_flipped(ca, midpoint.x - font_offset.x - pt_to_center.x, 
                          midpoint.y - font_offset.y - pt_to_center.y,
                          fg->edges[i].label_backward);
    
    //cairo_move_to(ca, midpoint.x, midpoint.y);
    //cairo_line_to(ca, midpoint.x - font_offset.x, midpoint.y - font_offset.y);
    //cairo_stroke(ca);
    //cairo_move_to(ca, midpoint.x, midpoint.y);
    
    //decorate the arcs with arrows
    //forward arrow
    point_and_tangent_to_bezier(fg->verts[fg->edges[i].start].loc,
                                fg->verts[fg->edges[i].end].loc,
                                control1,
                                control2,
                                edge_arrow_position[i],
                                &midpoint,
                                &tangent);
    tangentNorm = sqrt((tangent.x * tangent.x) + (tangent.y * tangent.y));
    tangent.x /= tangentNorm;
    tangent.y /= tangentNorm;
    normal.x = tangent.y;
    normal.y = -tangent.x;
    cairo_set_line_width(ca, (rectangle_width - inside_rectangle_width)/2);
    arrowCenter.x = midpoint.x + normal.x * (inside_rectangle_width/2);
    arrowCenter.y = midpoint.y + normal.y * (inside_rectangle_width/2);
    fg_cairo_draw_arrow(ca, arrowCenter, tangent, 9);
    
    //backward arrow
    point_and_tangent_to_bezier(fg->verts[fg->edges[i].start].loc,
                                fg->verts[fg->edges[i].end].loc,
                                control1,
                                control2,
                                edge_arrow_position[i],
                                &midpoint,
                                &tangent);
    tangentNorm = sqrt((tangent.x * tangent.x) + (tangent.y * tangent.y));
    tangent.x /= -tangentNorm;
    tangent.y /= -tangentNorm;
    normal.x = tangent.y;
    normal.y = -tangent.x;
    cairo_set_line_width(ca, (rectangle_width - inside_rectangle_width)/2);
    arrowCenter.x = midpoint.x + normal.x * (inside_rectangle_width/2);
    arrowCenter.y = midpoint.y + normal.y * (inside_rectangle_width/2);
    fg_cairo_draw_arrow(ca, arrowCenter, tangent, 9);   
    
  }
  
  //draw the vertices
  /*
  cairo_set_source_rgb(ca, 1,1,1);
  cairo_set_line_width(ca, inside_rectangle_width/2);
  for (i=0; i<fg->num_verts; i++) {
    cairo_move_to(ca, fg->verts[i].loc.x, fg->verts[i].loc.y);
    cairo_arc(ca, fg->verts[i].loc.x, fg->verts[i].loc.y, inside_rectangle_width/3, 0, 2*3.14159);
    cairo_stroke(ca);
  }
  */
  
  cairo_set_source_rgb(ca, 0,0,0);
  cairo_set_line_width(ca, inside_rectangle_width/4);
  for (i=0; i<fg->num_verts; i++) {
    cairo_move_to(ca, fg->verts[i].loc.x, fg->verts[i].loc.y);
    cairo_arc(ca, fg->verts[i].loc.x, fg->verts[i].loc.y, inside_rectangle_width/6, 0, 2*3.14159);
    cairo_stroke(ca);
  }
  
}


    
    

//do two different starting points in two cyclic words read off the same word?
//the words can be the same (pointer)!
//the words can be different, but they must be the same length!
int cyclic_words_are_same(char* w1, int p1, char* w2, int p2, int len) {
  int i;
  for (i=0; i<len; i++) {
    if (w1[(p1+i)%len] != w2[(p2+i)%len]) {
      return 0;
    }
  }
  return 1;
}
  

//returns -1 if w[p1] < w[p2], etc, where w[p1] means w starting at p1
int cyclic_positions_order(char* w, int p1, int p2) {
  int wL = strlen(w);
  int i;
  for (i=0; i<wL; i++) {
    if (w[(p1+i)%wL] == w[(p2+i)%wL]) continue;
    if (w[(p1+i)%wL] < w[(p2+i)%wL]) {
      return -1;
    } else if (w[(p1+i)%wL] > w[(p2+i)%wL]) {
      return 1;
    }
  }
  return 0;
}

void rotate_word(char* w, int shift) {
  int wL = strlen(w);
  char* tempWord = (char*)malloc((wL+1)*sizeof(char));
  int i;
  for (i=0; i<wL; i++) {
    tempWord[i] = w[(shift+i)%wL];
  }
  for (i=0; i<wL; i++) {
    w[i] = tempWord[i];
  }
  free(tempWord);
}


void rotate_to_alpha_max(char* w) {
  //and we rotate the word until it is alphabetically maximal
  //first, find the shift
  int maximalShift = 0;
  int shift;
  int wL = strlen(w);
  for (shift=1; shift<wL; shift++) {
    if (cyclic_positions_order(w, maximalShift, shift)==-1) { //shift is better
      maximalShift = shift;
    }
  }
  //now rotate it
  rotate_word(w, maximalShift);
}


//write a cyclic word in standard form:
//w^n ->(n,w), where n, w are such that n in maximal
//and w is rotated to be alphabetically maximal (and thus unique)
void standard_cyclic_form(char* w, int* n) {
  //find the smallest period -- this will be the smallest integer such that 
  //a shift of this gives the same word
  int period;
  int wL = strlen(w);
  for (period = 1; period < wL; period++) { //not all of these are necessary
    if (1==cyclic_words_are_same(w, 0, w, period, wL)) {
      break;
    }
  }
  //ok now the word is just the first part
  w[period] = '\0';
  *n = wL/period; //period divides the length
  wL = strlen(w); //=period-1
  
  rotate_to_alpha_max(w);
  
}
  
 

//combine the boundaries
//this starts at the first boundary and gets the standard cyclic form and searches
//through the remaining boundaries, picking up any copies it finds.  It deletes
//this from the list
//numBoundaries is updated at the end to reflect the new length
//this DOES free the memory for the removed words
void collect_boundaries(char** boundaries, int* numBoundaries, int* n) {
  int i,j,k;
  for (i=0; i<*numBoundaries; i++) {
    standard_cyclic_form(boundaries[i], &(n[i]));
  }
  for (i=0; i<*numBoundaries; i++) { //note numBoundaries will change
    j=i+1;
    while (j < *numBoundaries) {
      if (0==strcmp(boundaries[i], boundaries[j])) { 
        n[i] += n[j];
        //remove boundary j
        free(boundaries[j]);
        for (k=j; k<*numBoundaries-1; k++) {
          n[k] = n[k+1];
          boundaries[k] = boundaries[k+1];
        }
        (*numBoundaries)--;
      } else {
        j++;
      }
    }
  }
  
  
}
    




void fatgraph_update_data_fields(fatgraph* fg, widgetList* widgets) {
  GtkTextBuffer* dataView = gtk_text_view_get_buffer(widgets->dataTextView);
  char* newDataText = (char*)malloc(2000*sizeof(char));
  char tempString[20];
  char** boundaries;
  int* n = NULL;
  int numBoundaries;
  int i;
  int chi = fatgraph_chi(fg);
  strcpy(newDataText, "chi: ");
  sprintf(tempString, "%d\n", chi);
  strcat(newDataText, tempString);
  strcat(newDataText, "Boundaries:\n");
  boundaries = fatgraph_boundaries(fg, &numBoundaries);
  
  if (numBoundaries > 0) {
    for (i=0; i<numBoundaries; i++) {
      rotate_to_alpha_max(boundaries[i]);
      strcat(newDataText, "\t");
      strcat(newDataText, boundaries[i]);
      strcat(newDataText, "\n");
    }
    //now collect the boundaries
    n = (int*)malloc(numBoundaries*sizeof(int));
    int g,g2;
    collect_boundaries(boundaries, &numBoundaries, n);
    g = list_gcd(n, numBoundaries);
    strcat(newDataText, "Exhibits:\n\tscl(");
    if (n[0] != g) {
      sprintf(tempString, "%d", n[0]/g);
      strcat(newDataText, tempString);
    }
    strcat(newDataText, boundaries[0]);
    for (i=1; i<numBoundaries; i++) {
      strcat(newDataText, "\n\t     + ");
      if (n[i] != g) {
        sprintf(tempString, "%d", n[i]/g);
        strcat(newDataText, tempString);
      }
      strcat(newDataText, boundaries[i]);
      free(boundaries[i]);
    }
    strcat(newDataText, ") <= ");
    g2 = gcd(-chi, 2*g);
    if (2*g == g2) {
      sprintf(tempString, "%d\n", (-chi)/g2);
    } else {
      sprintf(tempString, "%d/%d\n", (-chi)/g2, (2*g)/g2);
    }
    strcat(newDataText, tempString);
  }

  
  
  
  free(boundaries);
  free(n);
  gtk_text_buffer_set_text(dataView, newDataText, -1);
  free(newDataText);
}
  




/******************************************************************************/
/* event handling                                                             */
/******************************************************************************/
static gboolean load_button_press(GtkWidget* widget, 
                                  GdkEventButton* event, 
                                  widgetList* widgets) {
  fatgraph* newFG;
  if (strlen((char*)gtk_entry_get_text(widgets->filenameEntry)) == 0) {
    //printf("The filename box was empty, so I didn't do anything\n");
    return TRUE;
  }
  if (event->type == GDK_LEAVE_NOTIFY) {
    //make the button pop up
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_NORMAL);
  } else if (event->type == GDK_BUTTON_PRESS) {
    //make the button go down
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_HALF);
  } else if (event->type == GDK_BUTTON_RELEASE) {
    //actually load the file 
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_NORMAL);    
    //load the fatgraph
    newFG = read_fatgraph_from_file((char*)gtk_entry_get_text(widgets->filenameEntry),
                                    400,
                                    400);
    if (newFG == NULL) {
      printf("That file doesn't seem to exist\n");
    } else {
      printf("Read graph\n");
      //fatgraph_print(newFG);
      fatgraph_copy(fg, newFG);
      //printf("Copied into fg:\n");
      //fatgraph_print(fg);
      fatgraph_draw(fg, cr, widgets->widget->allocation.width, widgets->widget->allocation.height, 1,0);
      gtk_widget_queue_draw_area(widgets->widget, 0, 0, widgets->widget->allocation.width,
                                                        widgets->widget->allocation.height);
      fatgraph_update_data_fields(fg, widgets);
    }
  }
  return TRUE;
}

static gboolean save_button_press(GtkWidget* widget, 
                                  GdkEventButton* event, 
                                  widgetList* widgets) {
  if (strlen((char*)gtk_entry_get_text(widgets->filenameEntry)) == 0) {
    //printf("The filename box was empty, so I didn't do anything\n");
    return TRUE;
  }
  if (event->type == GDK_LEAVE_NOTIFY) {    //make the button go down
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_HALF);
  } else if (event->type == GDK_BUTTON_RELEASE) {
    //actually load the file 
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_NORMAL);    
    //load the fatgraph
    int result = write_fatgraph_to_file(fg, (char*)gtk_entry_get_text(widgets->filenameEntry));
    if (result == 1) {
      printf("I couldn't open that file\n");
    }
    printf("Saved graph\n");
  }
  return TRUE;
}

static gboolean save_button_eps_press(GtkWidget* widget, 
                                      GdkEventButton* event, 
                                      widgetList* widgets) {
  if (strlen((char*)gtk_entry_get_text(widgets->filenameEntry)) == 0) {
    //printf("The filename box was empty, so I didn't do anything\n");
    return TRUE;
  }
  if (event->type == GDK_LEAVE_NOTIFY) {    //make the button go down
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_HALF);
  } else if (event->type == GDK_BUTTON_RELEASE) {
    //actually load the file 
    gtk_button_set_relief(GTK_BUTTON(widget), GTK_RELIEF_NORMAL);    
    //load the fatgraph
    int result = draw_fatgraph_to_file(fg, (char*)gtk_entry_get_text(widgets->filenameEntry),
                                           widgets->widget->allocation.width,
                                           widgets->widget->allocation.height);
    if (result == 1) {
      printf("I couldn't open that file\n");
    }
    printf("Drew graph\n");
  }
  return TRUE;
}



static gboolean opArcs_button_press(GtkWidget* widget, 
                                  GdkEventButton* event, 
                                  widgetList* widgets) {
  if (event->type == GDK_BUTTON_RELEASE) {
    fatgraph_optimize_drawing(fg, 
                              widgets->widget->allocation.width,
                              widgets->widget->allocation.height, 
                              0);
    fatgraph_draw(fg, cr, widgets->widget->allocation.width, widgets->widget->allocation.height, 1,0);
    gtk_widget_queue_draw_area(widgets->widget, 0, 0, widgets->widget->allocation.width,
                                                   widgets->widget->allocation.height);
  }
  return TRUE;
}



static gboolean opEverything_button_press(GtkWidget* widget, 
                                  GdkEventButton* event, 
                                  widgetList* widgets) {
  if (event->type == GDK_BUTTON_RELEASE) {
    fatgraph_optimize_drawing(fg, 
                              widgets->widget->allocation.width,
                              widgets->widget->allocation.height, 
                              1);
    fatgraph_draw(fg, cr, widgets->widget->allocation.width, widgets->widget->allocation.height, 1,0);
    gtk_widget_queue_draw_area(widgets->widget, 0, 0, widgets->widget->allocation.width,
                                                   widgets->widget->allocation.height);
  }
  return TRUE;
}

static gboolean relabelEdge_button_press(GtkWidget* widget, 
                                          GdkEventButton* event, 
                                          widgetList* widgets) {
  if (event->type == GDK_BUTTON_RELEASE) {
    if (dstatus.flag == FG_EDGE_HIGHLIGHT) {
      fg->edges[dstatus.highlightedEdge].label_forward = 
        (char*)realloc((void*)(fg->edges[dstatus.highlightedEdge].label_forward),
                       (strlen(gtk_entry_get_text(widgets->forwardLabelEntry))+1)*sizeof(char));
      fg->edges[dstatus.highlightedEdge].label_backward = 
        (char*)realloc((void*)(fg->edges[dstatus.highlightedEdge].label_backward),
                       (strlen(gtk_entry_get_text(widgets->backwardLabelEntry))+1)*sizeof(char));
      strcpy(fg->edges[dstatus.highlightedEdge].label_forward,
             gtk_entry_get_text(widgets->forwardLabelEntry));
      strcpy(fg->edges[dstatus.highlightedEdge].label_backward,
             gtk_entry_get_text(widgets->backwardLabelEntry));
      //unselect the edge
      dstatus.flag = FG_NO_FLAG;
      gtk_entry_set_text(widgets->forwardLabelEntry, "");
      gtk_entry_set_text(widgets->backwardLabelEntry, "");
      fatgraph_draw(fg, cr, widgets->widget->allocation.width, widgets->widget->allocation.height, 1,0);
      gtk_widget_queue_draw_area(widgets->widget, 0, 0, widgets->widget->allocation.width,
                                                     widgets->widget->allocation.height); 
      //printf("I relabeled an edge with %s and %s; now the graph is:\n", gtk_entry_get_text(widgets->forwardLabelEntry),
      //                                                                  gtk_entry_get_text(widgets->backwardLabelEntry));
      //fatgraph_print(fg);
      fatgraph_update_data_fields(fg, widgets);     
    }
  }
  return TRUE;
}
                                  
                                  
static gboolean scallopRun_button_press(GtkWidget* widget, 
                                        GdkEventButton* event, 
                                        widgetList* widgets) {        
  char command[300];
  fatgraph* newFG;
  if (event->type == GDK_BUTTON_RELEASE) {
    if (strlen(gtk_entry_get_text(widgets->scallopEntry)) > 0) {
      //run scallop
      strcpy(command, "./scallop/scallop -s fatgraphSurface ");
      strcat(command, gtk_entry_get_text(widgets->scallopEntry));
      if (0 != system(command)) {
        printf("Error running scallop\n");
        return TRUE;
      }
      printf("Ran successfully\n");
      
      //load the file
      newFG = read_fatgraph_from_file("fatgraphSurface.fg",
                                      400,
                                      400);
      fatgraph_copy(fg, newFG);
      //printf("Copied into fg:\n");
      //fatgraph_print(fg);
      fatgraph_draw(fg, cr, widgets->widget->allocation.width, widgets->widget->allocation.height, 1,0);
      gtk_widget_queue_draw_area(widgets->widget, 0, 0, widgets->widget->allocation.width,
                                                     widgets->widget->allocation.height); 
      fatgraph_update_data_fields(fg, widgets); 
    }
      
  }
  return TRUE;
}


static gint configure_event(GtkWidget* widget, GdkEventConfigure* event) {
  if (pixmap) 
    g_object_unref(pixmap);
  //printf("configure event has been called\n"); fflush(stdout);
  pixmap = gdk_pixmap_new(widget->window,
                          widget->allocation.width,
                          widget->allocation.height,
                          -1);
  cr = gdk_cairo_create(pixmap);
  //we need to apply a tranformation matrix to get into normal orientation
  cairo_matrix_t trans;
  trans.xx = 1;
  trans.xy = 0;
  trans.x0 = 0;
  trans.yx = 0;
  trans.yy = -1;
  trans.y0 = widget->allocation.height;
  cairo_set_matrix(cr, &trans); 
  
  
  
  cairo_set_source_rgb(cr, 1,1,1);
  cairo_rectangle(cr, 0,0,widget->allocation.width,widget->allocation.height);
  cairo_fill(cr);
  cairo_set_line_width(cr, 6.0);
  cairo_set_source_rgb(cr, 0,0,0);
  fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 0,0);
  return TRUE;
}

static gboolean expose_event_callback(GtkWidget* widget, GdkEventExpose* event) {
  //printf("About to try to draw the widget\n"); fflush(stdout);
  //printf("The value of pixmap: %x\n", (int)pixmap); fflush(stdout);
  gdk_draw_drawable(widget->window,
                    widget->style->fg_gc[gtk_widget_get_state (widget)],
                    pixmap,
                    event->area.x, event->area.y,
                    event->area.x, event->area.y,
                    event->area.width, event->area.height);
  //printf("I drew it\n"); fflush(stdout);
  return FALSE;
}
  
  
 
/*
static void draw_circle(GtkWidget* widget, gdouble x, gdouble y) {
  cairo_move_to(cr, x,y);
  cairo_arc(cr, x, y, 3, 0, 2*3.141519);
  cairo_stroke(cr);
  gtk_widget_queue_draw_area(widget, x-30, y-30, 60, 60);
  //g_print("I'm drawing a line\n");
}
*/



  
static gboolean drawing_button_press(GtkWidget* widget, 
                                     GdkEventButton* event,
                                     widgetList* widgets) {
  int i,v1,v2,j,destination;
  vector2d sbezier;
  vector2d ebezier;
  double our_theta, first_theta, current_theta;
  int index_in_s, index_in_e;
  double k,m1,m2;
  int edge_to_label;
  double eventFlipX = event->x;
  double eventFlipY = widgets->widget->allocation.height - event->y;
  
  //when a button is pressed:  if it's button 1 down on a vertex, we're dragging
  //if it's button 1 down not on a vertex, then make a new vertex
  //if it's button 3 down, we're starting an edge
  //if it's button 1 up and we're dragging, then stop dragging
  
  if (pixmap == NULL) { //what?
    return TRUE;
  }
  
  //g_print("button: %d, event->type = %d\n", event->button, event->type);
  
  if (event->button == 1) { //left button
    if (event->type == GDK_BUTTON_PRESS) {
    
      //find out if we're over a vertex
      for (i=0; i<fg->num_verts; i++) {
        if ((eventFlipX - fg->verts[i].loc.x)*(eventFlipX - fg->verts[i].loc.x) + 
            (eventFlipY - fg->verts[i].loc.y)*(eventFlipY - fg->verts[i].loc.y) < 40) {
          break;
        }
      }
      if (i < fg->num_verts) { //we're dragging
        dstatus.action = FG_DRAGGING_VERT;
        dstatus.dragging_vert = i;
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,1);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                             widget->allocation.height);
              
      
      } else if (0 <= (edge_to_label = fatgraph_close_to_edge(fg, eventFlipX, eventFlipY))) { 
        //handle this 
        if (dstatus.flag == FG_EDGE_HIGHLIGHT) {
          if(dstatus.highlightedEdge == edge_to_label) {
            dstatus.flag = FG_NO_FLAG;
            //unset the text boxes
            gtk_entry_set_text(widgets->forwardLabelEntry, "");
            gtk_entry_set_text(widgets->backwardLabelEntry, "");
            fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,0);
            gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
            return TRUE;
          } else {
            //unset the other edge
            //I don't need to do anything here?
          }
        }
        dstatus.flag = FG_EDGE_HIGHLIGHT;
        dstatus.highlightedEdge = edge_to_label;
        gtk_entry_set_text(widgets->forwardLabelEntry, fg->edges[edge_to_label].label_forward);
        gtk_entry_set_text(widgets->backwardLabelEntry, fg->edges[edge_to_label].label_backward);
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,0);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
        
        
        
      
      } else {
        //we're not by a vertex or an edge -- make a new vertex
        //g_print("making a new vertex\n");
        fatgraph_add_vertex(fg, eventFlipX, eventFlipY);
        printf("Making a new vertex at %f, %f\n", eventFlipX, eventFlipY);
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 0,0);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
        fatgraph_update_data_fields(fg, widgets);
        
      }
    } else { //it's button 1 up
      if (dstatus.action == FG_NO_ACTION) {
        //do nothing
      } else if (dstatus.action == FG_DRAGGING_VERT) {
        //stop dragging the vertex
        fg->verts[dstatus.dragging_vert].loc.x = eventFlipX;
        fg->verts[dstatus.dragging_vert].loc.y = eventFlipY;
        dstatus.action = FG_NO_ACTION;
        
        //optimize the graph!
        fatgraph_optimize_drawing(fg, 
                                  widget->allocation.width,
                                  widget->allocation.height, 
                                  0);
        
        
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1, 0);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
      } else {
        //reset to no action
        dstatus.action = FG_NO_ACTION;
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 0,0);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
      }
    }
  
  
  
 
  } else if (event->button == 3) { //it's button 3
    
    
    if (event->type == GDK_BUTTON_PRESS) {
      if (dstatus.action == FG_NO_ACTION) {
        dstatus.edge1_down_x = eventFlipX;
        dstatus.edge1_down_y = eventFlipY;
        dstatus.action = FG_DRAWING_EDGE1;
      } else if (dstatus.action == FG_DRAWING_EDGE1) {
        dstatus.edge2_down_x = eventFlipX;
        dstatus.edge2_down_y = eventFlipY;
        dstatus.action = FG_DRAWING_EDGE2;
      } else {
        //what?
      }
    
    
    } else if (event->type == GDK_BUTTON_RELEASE) {
      if (dstatus.action == FG_DRAWING_EDGE1) { //first edge release
        dstatus.edge1_up_x = eventFlipX;
        dstatus.edge1_up_y = eventFlipY;
        dstatus.flag = FG_DONE_EDGE1;

      } else if (dstatus.action == FG_DRAWING_EDGE2) { //finish the edge
        dstatus.action = FG_NO_ACTION;
        dstatus.flag = FG_NO_FLAG;
        
        dstatus.edge2_up_x = eventFlipX;
        dstatus.edge2_up_y = eventFlipY;
        
        //find the vertices
        m1 = m2 = 100000;
        //g_print("previous click: %f, %f\ncurrent click: %f,%f\n", tx, ty, eventFlipX, eventFlipY);
        for (i=0; i<fg->num_verts; i++) {
          k = SQUARE(dstatus.edge1_down_x - fg->verts[i].loc.x) +
              SQUARE(dstatus.edge1_down_y - fg->verts[i].loc.y);
          if ( k < m1 ) {
            v1 = i;
            m1 = k;
          }
          //g_print("k v2 for %d (at %f,%f): %f\n", i, fg->verts[i].x, fg->verts[i].y, k);
          //and for v1
          k = SQUARE(dstatus.edge2_down_x - fg->verts[i].loc.x) + 
              SQUARE(dstatus.edge2_down_y - fg->verts[i].loc.y);
          if (k < m2){
            v2 = i;
            m2 = k;
          }
          //g_print("k v1 for %d: %f\n", i, k);
        }
        
        sbezier.x = dstatus.edge1_up_x - dstatus.edge1_down_x;
        sbezier.y = dstatus.edge1_up_y - dstatus.edge1_down_y;
        ebezier.x = dstatus.edge2_up_x - dstatus.edge2_down_x;
        ebezier.y = dstatus.edge2_up_y - dstatus.edge2_down_y;
        
        //to figure out the indices, we go through the other beziers and stop when
        //the cross product has a sign change        
        if (fg->verts[v1].num_edges < 2) {
          index_in_s = 0;
        } else {
          first_theta = atan2(fg->verts[v1].bezier[0].y, 
                              fg->verts[v1].bezier[0].x);
          our_theta = atan2(sbezier.y, sbezier.x);
          for (i=1; i<fg->verts[v1].num_edges; i++) {
            //find if our theta is between first_theta and this theta
            current_theta = atan2(fg->verts[v1].bezier[i].y, 
                                  fg->verts[v1].bezier[i].x);
            if (in_circular_order(first_theta, our_theta, current_theta)==1) {
              break;
            }
          }
          index_in_s = i % (fg->verts[v1].num_edges);
        }
        if (fg->verts[v2].num_edges < 2) {
          index_in_e = 0;
        } else {
          first_theta = atan2(fg->verts[v2].bezier[0].y, 
                              fg->verts[v2].bezier[0].x);
          our_theta = atan2(ebezier.y, ebezier.x);
          for (i=1; i<fg->verts[v2].num_edges; i++) {
            //find if our theta is between first_theta and this theta
            current_theta = atan2(fg->verts[v2].bezier[i].y, 
                                  fg->verts[v2].bezier[i].x);
            if (in_circular_order(first_theta, our_theta, current_theta)==1) {
              break;
            }
          }
          index_in_e = i % (fg->verts[v2].num_edges);
        }
        fatgraph_add_edge(fg, v1, index_in_s, sbezier,
                              v2, index_in_e, ebezier);
        //printf("Added edge.  Current destinations:\n");
        for (i=0; i<fg->num_verts; i++) {
          for (j=0; j<fg->verts[i].num_edges; j++) {
            destination = ( fg->verts[i].edges_initial[j] == 1 ?
                            fg->edges[fg->verts[i].edges[j]].end :
                            fg->edges[fg->verts[i].edges[j]].start );
            //printf("vert %d, edge %d (%d): destination %d\n", i, j, fg->verts[i].edges[j], destination);
          }
        }
        fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,0);
        gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                                 widget->allocation.height);
        fatgraph_update_data_fields(fg, widgets);
      }
    }
  }
 
      
  return TRUE;
}

static gboolean drawing_motion_notify(GtkWidget* widget, GdkEventMotion* event) {
  int x,y;
  GdkModifierType state;
  
  if (event->is_hint) {
    gdk_window_get_pointer(event->window, &x, &y, &state);
  } else {
    x = event->x;
    y = event->y;
    state = event->state;
  }
  
  double eventFlipY = widget->allocation.height - event->y; 
  
  if (state && 
      GDK_BUTTON1_MASK && 
      pixmap != NULL && 
      dstatus.action == FG_DRAGGING_VERT) {
    fg->verts[dstatus.dragging_vert].loc.x = event->x;
    fg->verts[dstatus.dragging_vert].loc.y = eventFlipY;
    fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,1);
    gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                             widget->allocation.height);
  
  } else if (state && GDK_BUTTON3_MASK) {
    if (dstatus.action == FG_DRAWING_EDGE1 && dstatus.flag != FG_DONE_EDGE1) {
      fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,1);
      cairo_set_source_rgb(cr, 1,0,0);
      cairo_set_line_width(cr, 2.0);
      cairo_move_to(cr, dstatus.edge1_down_x, dstatus.edge1_down_y);
      cairo_line_to(cr, event->x, eventFlipY);
      cairo_stroke(cr);
      cairo_set_source_rgb(cr, 0,0,0);
      cairo_set_line_width(cr, 6.0);      
      gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                             widget->allocation.height);
    } else if (dstatus.action == FG_DRAWING_EDGE2 || 
                (dstatus.action == FG_DRAWING_EDGE1 && dstatus.flag == FG_DONE_EDGE1)) {
      fatgraph_draw(fg, cr, widget->allocation.width, widget->allocation.height, 1,1);
      cairo_set_source_rgb(cr, 1,0,0);
      cairo_set_line_width(cr, 2.0);
      cairo_move_to(cr, dstatus.edge1_down_x, dstatus.edge1_down_y);
      cairo_line_to(cr, dstatus.edge1_up_x, dstatus.edge1_up_y);
      cairo_stroke(cr);
      cairo_move_to(cr, dstatus.edge2_down_x, dstatus.edge2_down_y);
      cairo_line_to(cr, event->x, eventFlipY);
      cairo_stroke(cr);
      cairo_set_source_rgb(cr, 0,0,0);
      cairo_set_line_width(cr, 6.0);      
      gtk_widget_queue_draw_area(widget, 0, 0, widget->allocation.width,
                                             widget->allocation.height);
    }
  }
  
  
  
  return TRUE;
}

  
  

static gboolean delete_event(GtkWidget* widget, GdkEvent* event, gpointer data) {
  return FALSE;
}
static void destroy(GtkWidget* widget, gpointer data) {
  gtk_main_quit();
}


int main(int argc, char* argv[]) {
  
  fg = fatgraph_create();
  dstatus.action = FG_NO_ACTION;
  
  GtkWidget* window;
  GtkWidget* drawingArea;
  GtkWidget* vBox;
  GtkWidget* toolbarBox;
  GtkWidget* drawingAndDataBox;
  GtkWidget* sidebarVBox;
  GtkWidget* relabelingHBox;
  GtkWidget* forwardLabelEntry;
  GtkEntryBuffer* forwardLabelText;
  GtkWidget* backwardLabelEntry;
  GtkEntryBuffer* backwardLabelText;
  GtkWidget* relabelEdgeButton;
  GtkWidget* scallopEntry;
  GtkEntryBuffer* scallopText;
  GtkWidget* scallopRunButton;
  GtkWidget* dataBox;
  GtkWidget* opArcsButton;
  GtkWidget* opEverythingButton;
  GtkWidget* filenameEntry;
  GtkEntryBuffer* filenameText;
  GtkWidget* loadGraphButton;
  GtkWidget* saveGraphButton;
  GtkWidget* saveGraphEPSButton;
  
  
  gtk_init(&argc, &argv);
  
  
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  
  //this is the drawing area
  drawingArea = gtk_drawing_area_new();
  gtk_widget_set_size_request(drawingArea, 400, 400);
  
  //this box contains everything
  vBox = gtk_vbox_new(FALSE,0);
  
  //this box contains the drawing area and the data window
  drawingAndDataBox = gtk_hbox_new(FALSE, 0);
  
  //this box contains the toolbar buttons
  toolbarBox = gtk_hbox_new(FALSE, 0);
  
  //this is the side bar area
  sidebarVBox = gtk_vbox_new(FALSE, 0);
  relabelingHBox = gtk_hbox_new(FALSE,0);
  relabelEdgeButton = gtk_button_new_with_label("Relabel edge");
  forwardLabelText = gtk_entry_buffer_new(NULL, -1);
  forwardLabelEntry = gtk_entry_new_with_buffer(forwardLabelText);
  gtk_widget_set_size_request(forwardLabelEntry, 100, -1);
  backwardLabelText = gtk_entry_buffer_new(NULL, -1);
  backwardLabelEntry = gtk_entry_new_with_buffer(backwardLabelText);
  gtk_widget_set_size_request(backwardLabelEntry, 100, -1);
  scallopText = gtk_entry_buffer_new(NULL, -1);
  scallopEntry =  gtk_entry_new_with_buffer(scallopText);
  scallopRunButton = gtk_button_new_with_label("Load from scallop");
  
  dataBox = gtk_text_view_new();
  gtk_widget_set_size_request(dataBox, 100, 200);
  
  //these are the toolbar things
  opArcsButton = gtk_button_new_with_label("Smooth arcs");
  opEverythingButton = gtk_button_new_with_label("Redraw graph");
  filenameText = gtk_entry_buffer_new(NULL, -1);
  filenameEntry = gtk_entry_new_with_buffer(filenameText);
  loadGraphButton = gtk_button_new_with_label("Load");
  saveGraphButton = gtk_button_new_with_label("Save");
  saveGraphEPSButton = gtk_button_new_with_label("Save eps");
  
  //printf("I made all the stuff\n"); fflush(stdout);
  
  //assemble the sidebar
  gtk_box_pack_start(GTK_BOX(relabelingHBox), forwardLabelEntry, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(relabelingHBox), backwardLabelEntry, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(sidebarVBox), relabelingHBox, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sidebarVBox), relabelEdgeButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sidebarVBox), scallopEntry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sidebarVBox), scallopRunButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sidebarVBox), dataBox, TRUE, TRUE, 0);
  gtk_text_view_set_editable(GTK_TEXT_VIEW(dataBox), FALSE);
  
  
  //put the drawing area and data window into the big box
  gtk_box_pack_start(GTK_BOX(drawingAndDataBox), drawingArea, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(drawingAndDataBox), sidebarVBox, FALSE, FALSE, 0);
  
  //put the toolbar crap in
  gtk_box_pack_start(GTK_BOX(toolbarBox), opArcsButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(toolbarBox), opEverythingButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(toolbarBox), filenameEntry, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(toolbarBox), loadGraphButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(toolbarBox), saveGraphButton, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(toolbarBox), saveGraphEPSButton, FALSE, FALSE, 0);
  
  
  //assemble the big box
  gtk_box_pack_start(GTK_BOX(vBox), toolbarBox, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vBox), drawingAndDataBox, TRUE, TRUE, 0);
  
  //put it in the window
  gtk_container_add(GTK_CONTAINER(window), vBox);
  
 
  //printf("I packed everything\n"); fflush(stdout);
  
  g_signal_connect(GTK_OBJECT(window), 
                   "delete_event", 
                   G_CALLBACK(delete_event), 
                   NULL);
  g_signal_connect(GTK_OBJECT(window),
                   "destroy",
                   G_CALLBACK(destroy),
                   NULL);
  
  //button events
  widgetList widgets;
  widgets.filenameEntry = GTK_ENTRY(filenameEntry);
  widgets.forwardLabelEntry = GTK_ENTRY(forwardLabelEntry);
  widgets.backwardLabelEntry = GTK_ENTRY(backwardLabelEntry);
  widgets.scallopEntry = GTK_ENTRY(scallopEntry);
  widgets.dataTextView = GTK_TEXT_VIEW(dataBox);
  widgets.widget = drawingArea;
  
  
  
  //load
  gtk_signal_connect(GTK_OBJECT(loadGraphButton),
                     "button_press_event",
                     G_CALLBACK(load_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(loadGraphButton),
                     "button_release_event",
                     G_CALLBACK(load_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(loadGraphButton),
                     "leave_notify_event",
                     G_CALLBACK(load_button_press),
                     &widgets);
  gtk_widget_add_events(loadGraphButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);
  //save
  gtk_signal_connect(GTK_OBJECT(saveGraphButton),
                     "button_press_event",
                     G_CALLBACK(save_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(saveGraphButton),
                     "button_release_event",
                     G_CALLBACK(save_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(saveGraphButton),
                     "leave_notify_event",
                     G_CALLBACK(save_button_press),
                     &widgets);
  gtk_widget_add_events(saveGraphButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);
  
  //save eps
  gtk_signal_connect(GTK_OBJECT(saveGraphEPSButton),
                     "button_press_event",
                     G_CALLBACK(save_button_eps_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(saveGraphEPSButton),
                     "button_release_event",
                     G_CALLBACK(save_button_eps_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(saveGraphEPSButton),
                     "leave_notify_event",
                     G_CALLBACK(save_button_eps_press),
                     &widgets);
  gtk_widget_add_events(saveGraphEPSButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);
  
  //smooth arcs
  gtk_signal_connect(GTK_OBJECT(opArcsButton),
                     "button_press_event",
                     G_CALLBACK(opArcs_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(opArcsButton),
                     "button_release_event",
                     G_CALLBACK(opArcs_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(opArcsButton),
                     "leave_notify_event",
                     G_CALLBACK(opArcs_button_press),
                     &widgets);
  gtk_widget_add_events(opArcsButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);  
  //optimize everything
  gtk_signal_connect(GTK_OBJECT(opEverythingButton),
                     "button_press_event",
                     G_CALLBACK(opEverything_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(opEverythingButton),
                     "button_release_event",
                     G_CALLBACK(opEverything_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(opEverythingButton),
                     "leave_notify_event",
                     G_CALLBACK(opEverything_button_press),
                     &widgets);
  gtk_widget_add_events(opEverythingButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);  
  
  
  
  //relabel edge
  gtk_signal_connect(GTK_OBJECT(relabelEdgeButton),
                     "button_press_event",
                     G_CALLBACK(relabelEdge_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(relabelEdgeButton),
                     "button_release_event",
                     G_CALLBACK(relabelEdge_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(relabelEdgeButton),
                     "leave_notify_event",
                     G_CALLBACK(relabelEdge_button_press),
                     &widgets);
  gtk_widget_add_events(relabelEdgeButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);   
  
  //run scallop
  gtk_signal_connect(GTK_OBJECT(scallopRunButton),
                     "button_press_event",
                     G_CALLBACK(scallopRun_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(scallopRunButton),
                     "button_release_event",
                     G_CALLBACK(scallopRun_button_press),
                     &widgets);
  gtk_signal_connect(GTK_OBJECT(scallopRunButton),
                     "leave_notify_event",
                     G_CALLBACK(scallopRun_button_press),
                     &widgets);
  gtk_widget_add_events(scallopRunButton, GDK_BUTTON_RELEASE_MASK 
                                       | GDK_BUTTON_PRESS_MASK
                                       | GDK_LEAVE_NOTIFY_MASK);    
  
  
  
  
  
  
  
  //drawing stuff
  
  g_signal_connect(GTK_OBJECT(drawingArea), 
                   "expose_event", 
                   G_CALLBACK(expose_event_callback), 
                   NULL);
  g_signal_connect(GTK_OBJECT(drawingArea),
                   "configure_event",
                   G_CALLBACK(configure_event),
                   NULL);
  g_signal_connect(GTK_OBJECT(drawingArea),
                   "button_press_event",
                   G_CALLBACK(drawing_button_press),
                   &widgets);
  g_signal_connect(GTK_OBJECT(drawingArea),
                   "button_release_event",
                   G_CALLBACK(drawing_button_press),
                   &widgets);
  g_signal_connect(GTK_OBJECT(drawingArea),
                   "motion_notify_event",
                   G_CALLBACK(drawing_motion_notify),
                   NULL);
  gtk_widget_add_events(drawingArea, GDK_CONFIGURE
                                    | GDK_BUTTON_PRESS_MASK 
                                    | GDK_BUTTON_RELEASE_MASK
                                    | GDK_POINTER_MOTION_MASK 
                                    | GDK_POINTER_MOTION_HINT_MASK);
  
  //show everything
  gtk_widget_show_all(window);  
 
  
  gtk_main();
  
  return 0;
}
  
  
  
