#ifndef _VECTOR2D_H
#define _VECTOR2D_H

typedef struct {
  double x;
  double y;
} vector2d;

vector2d vector2d_add(vector2d v, vector2d w);
vector2d vector2d_sub(vector2d v, vector2d w);
double vector2d_dot(vector2d a, vector2d b);
vector2d linear_interpolate(vector2d v, vector2d w, double t);
vector2d point_on_bezier(vector2d p1, 
                         vector2d p2, 
                         vector2d h1, 
                         vector2d h2,
                         double t);
void point_and_tangent_to_bezier(vector2d p1, 
                                 vector2d p2, 
                                 vector2d h1, 
                                 vector2d h2,
                                 double t,
                                 vector2d* pt,
                                 vector2d* tangent);
void sub_bezier(vector2d start, vector2d end, 
                vector2d bezier1, vector2d bezier2, double t1, double t2,
                vector2d* trunc_start, vector2d* trunc_end,
                vector2d* trunc_start_bezier, vector2d* trunc_end_bezier);
double distance_to_interior_bezier(vector2d p1, vector2d p2,
                                   vector2d control1, vector2d control2,
                                   double x, double y);
double t_val_for_min(double y1, double y2);
double distance_to_segment(vector2d p1, vector2d p2, double x, double y);
double distance_to_bezier_fast(vector2d p1, vector2d p2,
                               vector2d control1, vector2d control2,
                               double x, double y);


#endif
