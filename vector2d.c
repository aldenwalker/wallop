#include <math.h>

#include "vector2d.h"

vector2d vector2d_add(vector2d v, vector2d w) {
  vector2d ans;
  ans.x = v.x + w.x;
  ans.y = v.y + w.y;
  return ans;
}

vector2d vector2d_sub(vector2d v, vector2d w) {
  vector2d ans;
  ans.x = v.x - w.x;
  ans.y = v.y - w.y;
  return ans;
}


double vector2d_dot(vector2d a, vector2d b) {
  return a.x*b.x + a.y*b.y;
}



//compute (1-t)v + tw
vector2d linear_interpolate(vector2d v, vector2d w, double t) {
  vector2d ans;
  ans.x = (1-t)*v.x + t*w.x;
  ans.y = (1-t)*v.y + t*w.y;
  return ans;
}

vector2d point_on_bezier(vector2d p1, 
                         vector2d p2, 
                         vector2d h1, 
                         vector2d h2,
                         double t) {
  vector2d ans;
  vector2d point_on_p1h1 = linear_interpolate(p1,h1,t);
  vector2d point_on_h1h2 = linear_interpolate(h1,h2,t);
  vector2d point_on_h2p2 = linear_interpolate(h2,p2,t);
  vector2d point_on_p1h1_h1h2 = linear_interpolate( point_on_p1h1, 
                                                    point_on_h1h2, t );
  vector2d point_on_h1h2_h2p2 = linear_interpolate( point_on_h1h2, 
                                                    point_on_h2p2, t );
  ans = linear_interpolate(point_on_p1h1_h1h2, point_on_h1h2_h2p2, t);
  return ans;
}
 
void point_and_tangent_to_bezier(vector2d p1, 
                                 vector2d p2, 
                                 vector2d h1, 
                                 vector2d h2,
                                 double t,
                                 vector2d* pt,
                                 vector2d* tangent) {
  vector2d point_on_p1h1 = linear_interpolate(p1,h1,t);
  vector2d point_on_h1h2 = linear_interpolate(h1,h2,t);
  vector2d point_on_h2p2 = linear_interpolate(h2,p2,t);
  vector2d point_on_p1h1_h1h2 = linear_interpolate( point_on_p1h1, 
                                                    point_on_h1h2, t );
  vector2d point_on_h1h2_h2p2 = linear_interpolate( point_on_h1h2, 
                                                    point_on_h2p2, t );
  *pt = linear_interpolate(point_on_p1h1_h1h2, point_on_h1h2_h2p2, t);
  (*tangent).x = point_on_h1h2_h2p2.x - point_on_p1h1_h1h2.x;
  (*tangent).y = point_on_h1h2_h2p2.y - point_on_p1h1_h1h2.y;
}


void sub_bezier(vector2d start, vector2d end, 
                vector2d bezier1, vector2d bezier2, double t1, double t2,
                vector2d* trunc_start, vector2d* trunc_end,
                vector2d* trunc_start_bezier, vector2d* trunc_end_bezier){
  vector2d point_on_p1h1 = linear_interpolate(start, bezier1, t1);
  vector2d point_on_h1h2 = linear_interpolate(bezier1, bezier2, t1);
  vector2d point_on_h2p2 = linear_interpolate(bezier2, end, t1);
  vector2d point_on_p1h1_h1h2 = linear_interpolate( point_on_p1h1, 
                                                    point_on_h1h2, t1 );
  vector2d point_on_h1h2_h2p2 = linear_interpolate( point_on_h1h2, 
                                                    point_on_h2p2, t1 );
  vector2d t1ToEndp1 = linear_interpolate(point_on_p1h1_h1h2, 
                                          point_on_h1h2_h2p2,
                                          t1);
  vector2d t1ToEndh1 = point_on_h1h2_h2p2;
  vector2d t1ToEndh2 = point_on_h2p2;
  vector2d t1ToEndp2 = end;
  
  //repeat for t1ToEnd
  double t2new = (t2-t1)/(1-t1);  //note t1 shouldn't be 1-- that would be dumb
  point_on_p1h1 = linear_interpolate(t1ToEndp1, t1ToEndh1, t2new);
  point_on_h1h2 = linear_interpolate(t1ToEndh1, t1ToEndh2, t2new);
  point_on_h2p2 = linear_interpolate(t1ToEndh2, t1ToEndp2, t2new);
  point_on_p1h1_h1h2 = linear_interpolate( point_on_p1h1, 
                                                    point_on_h1h2, t2new );
  point_on_h1h2_h2p2 = linear_interpolate( point_on_h1h2, 
                                                    point_on_h2p2, t2new );
  *trunc_end = linear_interpolate(point_on_p1h1_h1h2, 
                                          point_on_h1h2_h2p2,
                                          t2new);
  *trunc_start = t1ToEndp1;
  *trunc_start_bezier = point_on_p1h1;
  *trunc_end_bezier = point_on_p1h1_h1h2;
}
  
  
  
  
  
    



double distance_to_interior_bezier(vector2d p1, vector2d p2,
                                   vector2d control1, vector2d control2,
                                   double x, double y){
  //this is totally crappy -- but why not, right?
  double t;
  vector2d bezPt;
  double dist=1e10; 
  double current_dist;
  for (t=0.1; t<=0.9; t += 0.005) {
    bezPt = point_on_bezier(p1,p2,control1,control2,t);
    current_dist = sqrt((bezPt.x-x)*(bezPt.x-x) + (bezPt.y-y)*(bezPt.y-y));
    //printf("Current distance: %f\n", current_dist);
    if (current_dist < dist) {
      dist = current_dist;
      if (dist < 20) {
        return dist;
      }
    }
  }
  return dist;          
}



double t_val_for_min(double y1, double y2) {
  if (fabs(y2-y1) < 0.0001) {
    return 0;
  } else {
    if (y2*y1 > 0) { //same sign
      return (fabs(y2) < fabs(y1) ? 1 : 0);
    }
    return fabs(y1/(y2-y1));
  }
}


double distance_to_segment(vector2d p1, vector2d p2, double x, double y) {
  double initial_deriv = (p1.x-x)*(p2.x-p1.x) + (p1.y-y)*(p2.y-p1.y);
  double final_deriv = (p2.x-y)*(p2.x-p1.x) + (p2.y-y)*(p2.y-p1.y);
  double t = t_val_for_min(initial_deriv, final_deriv);
  return sqrt( (t*p1.x + (1-t)*p2.x - x)*(t*p1.x + (1-t)*p2.x - x) +
               (t*p1.y + (1-t)*p2.y - y)*(t*p1.y + (1-t)*p2.y - y) );
}


double distance_to_bezier_fast(vector2d p1, vector2d p2,
                               vector2d control1, vector2d control2,
                               double x, double y) {
  //approximate the bezier by five points
  vector2d m1,m2,m3;
  m1 = point_on_bezier(p1,p2,control1, control2, 0.25);
  m2 = point_on_bezier(p1,p2,control1, control2, 0.5);
  m3 = point_on_bezier(p1,p2,control1, control2, 0.75);
  double d, d2;
  d = distance_to_segment(p1, m1, x, y);
  d2 = distance_to_segment(m1, m2, x, y);
  d = (d > d2 ? d2 : d);  
  d2 = distance_to_segment(m2, m3, x, y);
  d = (d > d2 ? d2 : d);  
  d2 = distance_to_segment(m3, p2, x, y);
  d = (d > d2 ? d2 : d);  
  return d;
}
