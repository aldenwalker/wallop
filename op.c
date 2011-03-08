#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "op.h"

#define PI 3.14159265

//my angles are in [-pi,pi]

double vector2d_to_theta(vector2d v) {
  return atan2(v.y, v.x); 
}

vector2d theta_to_vector2d(double theta) {
  vector2d ans;
  ans.x = cos(theta);
  ans.y = sin(theta);
  return ans;
}


void fatgraph_op_point_on_line(op_vert* dest, op_vert* ov, op_vert* dir, int len, double t) {
  int i,j;
  for (i=0; i<len; i++) {
    dest[i].loc.x = ov[i].loc.x + (t*(dir[i].loc.x));
    dest[i].loc.y = ov[i].loc.y + (t*(dir[i].loc.y));
    for (j=0; j<ov[i].num_edges; j++) {
      dest[i].theta[j] = ov[i].theta[j] + (t*(dir[i].theta[j]));
    }
  }  
}

void fatgraph_op_print(op_vert* ov, int len) {
  int i,j;
  for (i=0; i<len; i++) {
    printf("vert %d: %f, %f\n", i, ov[i].loc.x, ov[i].loc.y);
    for (j=0; j<ov[i].num_edges; j++) {
      printf("\t%f to edge: %d\n", ov[i].theta[j], ov[i].edges[j]);
    }
  }
}



int connected_verts(op_vert* ov, int i, int j) {
  int k;
  for (k=0; k<ov[i].num_edges; k++) {
    if (ov[i].dest_verts[k] == j) {
      //printf("I concluded that %d and %d were connected\n", i,j);
      return 1;
    }
  }
  //printf("I concluded that %d and %d were NOT connected\n", i,j);
  return 0;
}



void fatgraph_op_gradient(op_vert* ov, 
                          op_vert* grad, 
                          int len,
                          int op_positions,
                          double target_length,
                          double alpha) {
  int i,j,k;
  vector2d vi;
  vector2d vj;
  vector2d vjvi;
  vector2d bezj,bezk;
  double bezjbezk;
  int destvert;
  double normvjvi;
  double thetaj, thetak;
  int edge_num;
  
  //compute the gradient in each of the vertex indices
  //note that these formulas work for both vertices < i and > i
  for (i=0; i<len; i++) {
    grad[i].loc.x = 0;
    grad[i].loc.y = 0;
  }
  if (op_positions == 0) {
    goto skipPositions;
  }
  for (i=0; i<len; i++) {
    vi = ov[i].loc;
    //location of vj affects both 
    //difference parts:
    for (j=0; j<len; j++) {
      if (i==j) continue;
      vj = ov[j].loc;
      vjvi = vector2d_sub(vj,vi);
      normvjvi = sqrt(vector2d_dot(vjvi,vjvi));
      //printf("About to run connected_verts\n");
      if (connected_verts(ov, i,j)==1) {
        grad[i].loc.x = (1/(alpha*alpha)) * 
                        ( 2*(vj.x - vi.x) + 
                          2*target_length*(1/(2*normvjvi))*2*(vj.x-vi.x));
        grad[i].loc.y = (1/(alpha*alpha)) *
                        ( 2*(vj.y - vi.y) + 
                          2*target_length*(1/(2*normvjvi))*2*(vj.y-vi.y));
      } else {
        grad[i].loc.x = -target_length*(vj.x-vi.x)*(1/normvjvi)*
                                                     (1/normvjvi)*
                                                     (1/normvjvi);
        grad[i].loc.y = -target_length*(vj.y-vi.y)*(1/normvjvi)*
                                                     (1/normvjvi)*
                                                     (1/normvjvi);
      }
    }
    //and angle parts:
    //note a vertex location affects its and its neighbors functions
    for (j=0; j<ov[i].num_edges; j++) {
      destvert = ov[i].dest_verts[j];
      if (destvert == i) continue;
      vj = ov[destvert].loc;
      edge_num = ov[i].edges[j];
      for (k=0; k<ov[destvert].num_edges; k++) {
        if (ov[destvert].edges[k] == edge_num) {
          break;
        }
      }
      vjvi = vector2d_sub(vj,vi);
      normvjvi = sqrt(vector2d_dot(vjvi,vjvi));
      thetaj = ov[i].theta[j];
      thetak = ov[destvert].theta[k];
      bezj.x = cos(thetaj);
      bezj.y = sin(thetaj);
      bezk.x = cos(thetak);
      bezk.y = sin(thetak);
      //affecting it's own functions:
      grad[i].loc.x += (-cos(thetaj)/normvjvi) + 
                       (vector2d_dot(vjvi,bezj)*(1/normvjvi)*
                                                (1/normvjvi)*
                                                (1/normvjvi));
      grad[i].loc.y += (-sin(thetaj)/normvjvi) + 
                       (vector2d_dot(vjvi,bezj)*(1/normvjvi)*
                                                (1/normvjvi)*
                                                (1/normvjvi));
      //affecting the other end:
      grad[i].loc.x += (cos(thetak)/normvjvi) + 
                       (vector2d_dot(vjvi,bezj)*(1/normvjvi)*
                                                (1/normvjvi)*
                                                (1/normvjvi));
      grad[i].loc.y += (sin(thetak)/normvjvi) + 
                       (vector2d_dot(vjvi,bezj)*(1/normvjvi)*
                                                (1/normvjvi)*
                                                (1/normvjvi));
    }
  }
  
  
  //LABEL
  //printf("I did NOT skip positions\n");
  skipPositions:
  
  //now for the theta parts:
  for (i=0; i<len; i++) {
    vi = ov[i].loc;
    for (j=0; j<ov[i].num_edges; j++) {
      destvert = ov[i].dest_verts[j];
      thetaj = ov[i].theta[j];
      bezj.x = cos(thetaj);
      bezj.y = sin(thetaj);
      if (destvert != i) {
        //affects pointing towards target:
        vj = ov[destvert].loc;
        vjvi = vector2d_sub(vj,vi);
        normvjvi = sqrt(vector2d_dot(vjvi,vjvi));
        grad[i].theta[j] = (1/normvjvi)*(vjvi.x*(-sin(thetaj)) + vjvi.y*cos(thetaj));
      } else {
        grad[i].theta[j] = 0;
      }
    
      //and affects interactions between edges
      for (k=0; k<ov[i].num_edges; k++) {
        if (k==j) continue;
        thetak = ov[i].theta[k];
        bezk.x = cos(thetak);
        bezk.y = sin(thetak);
        bezjbezk = vector2d_dot(bezj,bezk);
        if (ov[i].edges[j] == ov[i].edges[k]) { //loop edge
          //printf("Loop edge!\n");
          grad[i].theta[j] += (-6)*bezjbezk*(1/(1-bezjbezk))*
                              (2 + bezjbezk*(1/(1-bezjbezk)))*
                              (bezk.x*(-sin(thetaj)) + bezk.y*cos(thetaj));
        } else { //normal pair
          grad[i].theta[j] += (-1/((1-bezjbezk)*(1-bezjbezk))) *
                              (bezk.x * (-sin(thetaj)) + bezk.y * cos(thetaj));
        }
      }
    
    }
  }
  
  
  
}
        
                                    
      
                                                    
    
  
  
  
  
/*               
void fatgraph_op_gradient(op_vert* ov, 
                          op_vert* grad, 
                          int len,
                          int op_positions,
                          double target_length,
                          double alpha_scale) {
  //in each slot, put that partial derivative
  int i,j,k;
  double destx,desty,th,vertx,verty;
  double bjbk,bjbkm1;
  double diffNorm;
  for (i=0; i<len; i++) {
    vertx = ov[i].loc.x;
    verty = ov[i].loc.y;
    grad[i].loc.x = 0;
    grad[i].loc.y = 0;
    for (j=0; j<ov[i].num_edges; j++) {
      if (ov[i].dest_verts[j] == i) {
        continue;
      }
      destx = ov[ov[i].dest_verts[j]].loc.x;
      desty = ov[ov[i].dest_verts[j]].loc.y;
      //printf("Destination: %d at %f, %f\n", ov[i].dest_verts[j], destx,desty);
      
      if (op_positions==1) {
        diffNorm = (destx-vertx)*(destx-vertx) + (desty-verty)*(desty-verty);
        grad[i].loc.x += (-1/(alpha_scale*alpha_scale)) * 
                           (2*(vertx - destx) 
                            - 2*target_length*((vertx-destx)/sqrt(diffNorm)) );      
        grad[i].loc.y += (-1/(alpha_scale*alpha_scale)) * 
                           (2*(verty - desty) 
                            - 2*target_length*((verty-desty)/sqrt(diffNorm)) ); 

      }
    }
    
    for (j=0; j<ov[i].num_edges; j++) {
      destx = ov[ov[i].dest_verts[j]].loc.x;
      desty = ov[ov[i].dest_verts[j]].loc.y;     
      th = ov[i].theta[j];
      if (ov[i].dest_verts[j] == i) {
        grad[i].theta[j] = 0;
      } else {
        grad[i].theta[j] = 2*((destx-vertx)*(-sin(th)) + (desty-verty)*cos(th)) /
                            sqrt((destx-vertx)*(destx-vertx) + (desty-verty)*(desty-verty));
      }
      //printf("Dot with theta: %f\n", grad[i].theta[j]);
      for (k=0; k<ov[i].num_edges; k++) {
        if (k==j) continue;
        bjbk = cos(ov[i].theta[k])*cos(th) + sin(ov[i].theta[k])*sin(th);
        bjbkm1 = cos(ov[i].theta[k])*cos(th-1) + sin(ov[i].theta[k])*sin(th-1);
        //printf("Other theta: %f - dot is %f\n", ov[i].theta[k], bjbk);
        if (ov[i].dest_verts[k] == i) {
          grad[i].theta[j] -= (2/((1-bjbkm1)*(1-bjbkm1))) * 
                             (cos(ov[i].theta[k])*(-sin(th-1)) + sin(ov[i].theta[k])*cos(th-1));
        } else {
          grad[i].theta[j] -= ((ov[i].num_edges==2 ? 10 : 2)/((1-bjbk)*(1-bjbk))) * 
                             (cos(ov[i].theta[k])*(-sin(th)) + sin(ov[i].theta[k])*cos(th));
        }
      }
    }
  }
}
*/




/*
The objective function is:
maximize!

sum_i,j  if joined: -1/alpha^2 ( ||v_i-v_j||^2 - 2D||v_i-v_j|| +(D-alpha)(D+alpha)
         if not:    1-D/||v_i-v_j||
         
         
sum_i,bez_j=dest_j  (v_j-v_i).bez_j/||v_j-v_i||

sum_v_i,bez_j,bez_k  
        if not dest_j=dest_k=v_i:  -1/(1-b_j.b_k)
                           if so:  -2(b_j.b_k)^2/(1-b_j.b_k)


*/
double fatgraph_op_objective_value(op_vert* ov, 
                                   int len, 
                                   int op_positions,
                                   double target_length, 
                                   double alpha) {
  int i,j,k;
  vector2d vi;
  vector2d vj;
  vector2d vjvi;
  vector2d bezj,bezk;
  int destvert;
  double normvjvi;
  double bjbkdot;
  double val = 0;
  //go through and compute the vert distance parts
  if (op_positions==1) {
    for (i=0; i<len; i++) {
      vi.x = ov[i].loc.x;
      vi.y = ov[i].loc.y;
      for (j=i+1; j<len; j++) {
        vj.x = ov[j].loc.y;
        vj.y = ov[i].loc.y;
        normvjvi = sqrt((vi.x-vj.x)*(vi.x-vj.x) + (vi.y-vj.y)*(vi.y-vj.y));
        //printf("About to run connected_verts\n");
        if (connected_verts(ov, i,j)==1) {
          val += (-1/(alpha*alpha)) *
                 (normvjvi*normvjvi - 2*target_length*normvjvi + (target_length + alpha)*
                                                          (target_length - alpha)); 
        } else { //the verts aren't connected
          val += 1 - target_length/normvjvi;
        }
      }
    }
  }
  
  //make the edges want to point towards their destination
  for (i=0; i<len; i++) {
    vi.x = ov[i].loc.x;
    vi.y = ov[i].loc.y;
    for (j=0; j<ov[i].num_edges; j++) {
      destvert = ov[i].dest_verts[j]; 
      if (destvert == i) {
        continue;
      }
      vj.x = ov[destvert].loc.x;
      vj.y = ov[destvert].loc.y;
      vjvi.x = vj.x - vi.x;
      vjvi.y = vj.y - vi.y;
      normvjvi = sqrt(vector2d_dot(vjvi,vjvi));
      bezj.x = cos(ov[i].theta[j]);
      bezj.y = sin(ov[i].theta[j]);
      val += (1/normvjvi)*vector2d_dot(vjvi, bezj);
    }
  }
  
  //make the angles want to stay apart
  for (i=0; i<len; i++) {
    vi.x = ov[i].loc.x;
    vi.y = ov[i].loc.y;
    for (j=0; j<ov[i].num_edges; j++) {
      bezj.x = cos(ov[i].theta[j]);
      bezj.y = sin(ov[i].theta[j]);
      for (k=j+1; k<ov[i].num_edges; k++) {
        bezk.x = cos(ov[i].theta[k]);
        bezk.y = sin(ov[i].theta[k]);
        

        bjbkdot = vector2d_dot(bezj,bezk);
        if (ov[i].edges[j] == ov[i].edges[k]) {        
          //if they are the same edge, bias it towards not pi
          val +=  (-6*bjbkdot*bjbkdot)/(1-bjbkdot);
        } else { 
          //if not, do normal
          val += (-2)/(1-bjbkdot);
        }
      }
    }
  }
  
  
  return val;
          
}
        
  
  
                                   
/* old version                                
double fatgraph_op_objective_value(op_vert* ov, 
                                   int len, 
                                   int op_positions,
                                   double target_length, 
                                   double alpha_scale) {
  int i,j,k;
  double destx,desty,th,vertx,verty;
  double diffNorm;
  double val = 0;
  for (i=0; i<len; i++) {
    vertx = ov[i].loc.x;
    verty = ov[i].loc.y;
    for (j=0; j<ov[i].num_edges; j++) {
      destx = ov[ov[i].dest_verts[j]].loc.x;
      desty = ov[ov[i].dest_verts[j]].loc.y;
      
      th = ov[i].theta[j];
      
      if (ov[i].dest_verts[j] != i) {
        if (op_positions == 1) {
          diffNorm = (destx-vertx)*(destx-vertx) + (desty-verty)*(desty-verty);
          val += (-1/(alpha_scale*alpha_scale)) * 
                  (diffNorm 
                   - 2*target_length*sqrt(diffNorm) 
                   + (target_length-alpha_scale)*(target_length+alpha_scale));  
        }    

        val += 2 * ((destx-vertx)*(cos(th)) + (desty-verty)*sin(th)) /
                sqrt((destx-vertx)*(destx-vertx) + (desty-verty)*(desty-verty));
      }
      //printf("val from vertex %d -> %d: %f\n", i,  ov[i].dest_verts[j], val ); fflush(stdout);
      for (k=j+1; k<ov[i].num_edges; k++) {
        //printf("The bottom is: cos(%f)*cos(%f) + sin(%f)*sin(%f) = %f\n", ov[i].theta[k], th,ov[i].theta[k], th, cos(ov[i].theta[k])*cos(th) + sin(ov[i].theta[k])*sin(th)); fflush(stdout);
        if (ov[i].dest_verts[k] == i) {
          val -= 2/(1-(cos(ov[i].theta[k])*cos(th-1) + sin(ov[i].theta[k])*sin(th-1)));
        } else {
          val -= (ov[i].num_edges==2 ? 10 : 2)/(1-(cos(ov[i].theta[k])*cos(th) + sin(ov[i].theta[k])*sin(th)));
        }
      }
    }
    //printf("val after vertex %d: %f\n", i, val); fflush(stdout);
  }
  return val;
} 
*/



void fatgraph_op_line_max(op_vert* ov, 
                          op_vert* dir, 
                          op_vert* temp, 
                          int len,
                          int op_positions,
                          double target_length,
                          double alpha_scale) {
  double f0, f1, f2, fnew;
  double t0, t1, t2, tnew;
  double tt0, tt1, tt2;
  double ff0, ff1, ff2;
  double tol = 0.0001;
  int num_tries;
  
  //first we need to find the maximum value that t can be --
  //note we can't cross vectors, and no point in spinning around
  double max_t=1e10;
  double diff0, diff1, diff2;
  int i,j,k;
  for (i=0; i<len; i++) {
    for (j=0; j<dir[i].num_edges; j++) {
      if (dir[i].theta[j] > 0) {
        max_t = (2*PI/dir[i].theta[j] < max_t ? 2*PI/dir[i].theta[j] : max_t);
      } else if (dir[i].theta[j] < 0) {
        max_t = (-2*PI/dir[i].theta[j] < max_t ? -2*PI/dir[i].theta[j] : max_t);
      }
    }
  }
  if (max_t == 1e10) {
    printf("No edges?\n");
    return;
  }
  for (i=0; i<len; i++) {
    for (j=0; j<ov[i].num_edges; j++) {
      for (k=j+1; k<ov[i].num_edges; k++) {
        diff0 = ((ov[i].theta[j] - ov[i].theta[k])) / (dir[i].theta[j] - dir[i].theta[k]);
        diff1 = (2*PI - (ov[i].theta[j] - ov[i].theta[k])) / (dir[i].theta[j] - dir[i].theta[k]);
        diff2 = (-2*PI - (ov[i].theta[j] - ov[i].theta[k])) / (dir[i].theta[j] - dir[i].theta[k]);
        //printf("I get %f, with t mult of %f, so I get %f, %f, %f\n", (ov[i].theta[j] - ov[i].theta[k]),
        //                                                             (dir[i].theta[j] - dir[i].theta[k]),
        //                                                             diff0, diff1, diff2);
        max_t = (diff0 < max_t && diff0 > 0 ? diff0 : max_t);
        max_t = (diff1 < max_t && diff1 > 0 ? diff1 : max_t);
        max_t = (diff2 < max_t && diff2 > 0 ? diff2 : max_t);
      }
    }
  }
    
  //printf("max_t: %f\n", max_t); fflush(stdout);  
 // for (i=0; i<32; i++) {
  //  fatgraph_op_point_on_line(temp, ov, dir, len, (double)i*max_t/32.0);
  //  f0 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
 //   printf("%f -> %f\n", i*max_t/32.0, f0); fflush(stdout);
  //}
  
  t0 = 0;
  fatgraph_op_point_on_line(temp, ov, dir, len, 0);
  //fatgraph_op_print(temp, len);
  f0 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
  
  t1 = max_t/8.0;
  fatgraph_op_point_on_line(temp, ov, dir, len, t1);
  f1 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
  
  while (f1 < f0){
    t1 /= 2.0;
    fatgraph_op_point_on_line(temp, ov, dir, len, t1);
    f1 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
  }  
  
  t2 = (max_t-t1)/2.0;
  fatgraph_op_point_on_line(temp, ov, dir, len, t2);
  f2 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
  
  while (f2 > f1) {
    t2 += max_t/32.0;
    fatgraph_op_point_on_line(temp, ov, dir, len, t2);
    f2 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
  }
  
  //printf("Initial bracket: [%f,%f,%f] -> [%f,%f,%f]\n", t0,t1,t2,f0,f1,f2);
  
  num_tries = 0;
  while (fabs(t2-t0) > tol) {
    num_tries ++;
    if (num_tries > 20) {
      break;
    }
    //guess the parabolic maximum
    tnew = t1 - 0.5*( ( (t1-t0)*(t1-t0)*(f1-f2) - (t1-t2)*(t1-t2)*(f1-f0)) /
                      ( (t1-t0)*(f1-f2)         - (t1-t2)*(f1-f0) )          );
    //figure out if this is a reasonable guess
    if (t0 < tnew && tnew < t2) {
      //use this guess
      fatgraph_op_point_on_line(temp, ov, dir, len, tnew);
      fnew = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
      //printf("I'm using a parabola guess, giving %f at %f\n", fnew, tnew);
      
      //generally, this parabolic guess is extremely good -- let's try a new
      //bracket of size 1/10 of what we had before, centered at this point
      tt0 = tnew - (0.001*(t2-t0));
      tt1 = tnew;
      tt2 = tnew + (0.001*(t2-t0));
      fatgraph_op_point_on_line(temp, ov, dir, len, tt0);
      ff0 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
      ff1 = fnew;
      fatgraph_op_point_on_line(temp, ov, dir, len, tt2);
      ff2 = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
      //printf("Trying [%f,%f,%f] at [%f,%f,%f]\n", ff0,ff1,ff2,tt0,tt1,tt2);
      if (ff0 < ff1 && ff1 > ff2) {
        //printf("My super squish worked!\n");
        f0 = ff0;
        t0 = tt0;
        f1 = ff1;
        t1 = tt1;
        f2 = ff2;
        t2 = tt2;
        goto cutTheCrap;
      }
      
      if (tnew < t1) {
        if (fnew > f1) {
          t2 = t1;
          f2 = f1;
          t1 = tnew;
          f1 = fnew;
        } else {
          t0 = tnew;
          f0 = fnew;
        }
      } else { 
        if (fnew > f1) {
          t0 = t1;
          f0 = f1;
          t1 = tnew;
          f1 = fnew;
        } else {
          t2 = tnew;
          f2 = fnew;
        }
      }
        
    } else {
      //bisection
      //printf ("eh, bisection\n");
      if (fabs(t1-t0) > fabs(t2-t1)) {    
        tnew = t0 + (t1-t0)/2.0;
        fatgraph_op_point_on_line(temp, ov, dir, len, tnew);
        fnew = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
        if (fnew > f1) {
          t2 = t1;
          f2 = f1;
          t1 = tnew;
          f1 = fnew;
        } else {
          t0 = tnew;
          f0 = fnew;
        }
      } else {
        tnew = t1 + (t2-t1)/2.0;
        fatgraph_op_point_on_line(temp, ov, dir, len, tnew);
        fnew = fatgraph_op_objective_value(temp, len, op_positions, target_length, alpha_scale);
        if (fnew > f1) {
          t0 = t1;
          f0 = f1;
          t1 = tnew;
          f1 = fnew;
        } else {
          t2 = tnew;
          f2 = fnew;
        }
      }
    }
    cutTheCrap:;
    //printf("Current bracket: [%f,%f,%f] -> [%f,%f,%f]\n", t0,t1,t2,f0,f1,f2);
  }
  
  fatgraph_op_point_on_line(temp, ov, dir, len, t1);
  fatgraph_op_point_on_line(ov, temp, dir, len, 0);//copy temp to ov
    
}






void fatgraph_optimize_drawing(fatgraph* fg, 
                               double screen_width,
                               double screen_height, 
                               int optimize_positions_also) {
  int i,j;
  int num_verts = fg->num_verts;
  int total_dimension;
  int current_iteration;
  
  double objVal, prevObjVal;
  double gradNorm;
  double oldGradNorm;
  double beta;
  double tol = 0.001;
  double target_length, alpha_scale;
  
  printf ("Entered optimization function\n"); fflush(stdout);
  //for (i=0; i<16; i++) {
  //  th = ((double)i)*2*PI/16.0;
  //  printf("Theta: %f ", th);
  //  temp = theta_to_vector2d(th);
  //  printf(" -> (%f,%f) ", temp.x, temp.y);
  //  printf(" -> %f\n", vector2d_to_theta(temp));
  //}
    
  
  
  
  
  //make the fatgraph more amenable to optimization, and make a temp one
  //AND swap the y axis
  op_vert* ov = (op_vert*)malloc((fg->num_verts)*sizeof(op_vert));
  op_vert* grad = (op_vert*)malloc((fg->num_verts)*sizeof(op_vert));
  op_vert* oldGrad = (op_vert*)malloc((fg->num_verts)*sizeof(op_vert));
  op_vert* dir = (op_vert*)malloc((fg->num_verts)*sizeof(op_vert));
  op_vert* temp = (op_vert*)malloc((fg->num_verts)*sizeof(op_vert));
  total_dimension = 0;
  for (i=0; i<fg->num_verts; i++) {
    total_dimension += 2 + fg->verts[i].num_edges;
    ov[i].loc.x = fg->verts[i].loc.x;
    ov[i].loc.y = fg->verts[i].loc.y;
    ov[i].num_edges = fg->verts[i].num_edges;
    ov[i].dest_verts = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    ov[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    ov[i].theta = (double*)malloc((fg->verts[i].num_edges)*sizeof(double));
    grad[i].num_edges = fg->verts[i].num_edges;
    grad[i].dest_verts = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    grad[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    grad[i].theta = (double*)malloc((fg->verts[i].num_edges)*sizeof(double));
    oldGrad[i].num_edges = fg->verts[i].num_edges;
    oldGrad[i].dest_verts = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    oldGrad[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    oldGrad[i].theta = (double*)malloc((fg->verts[i].num_edges)*sizeof(double));
    dir[i].num_edges = fg->verts[i].num_edges;
    dir[i].dest_verts = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    dir[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    dir[i].theta = (double*)malloc((fg->verts[i].num_edges)*sizeof(double));
    temp[i].num_edges = fg->verts[i].num_edges;
    temp[i].dest_verts = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    temp[i].edges = (int*)malloc((fg->verts[i].num_edges)*sizeof(int));
    temp[i].theta = (double*)malloc((fg->verts[i].num_edges)*sizeof(double));
    for (j=0; j<fg->verts[i].num_edges; j++) {
      ov[i].dest_verts[j] = ( fg->edges[fg->verts[i].edges[j]].start == i ?
                              fg->edges[fg->verts[i].edges[j]].end :
                              fg->edges[fg->verts[i].edges[j]].start );
      ov[i].edges[j] = fg->verts[i].edges[j];
      grad[i].dest_verts[j] = ov[i].dest_verts[j];
      grad[i].edges[j] = ov[i].edges[j];
      oldGrad[i].dest_verts[j] = ov[i].dest_verts[j];
      oldGrad[i].edges[j] = ov[i].edges[j];
      dir[i].dest_verts[j] = ov[i].dest_verts[j];
      dir[i].edges[j] = ov[i].edges[j];
      temp[i].dest_verts[j] = ov[i].dest_verts[j];
      temp[i].edges[j] = ov[i].edges[j];
      ov[i].theta[j] = vector2d_to_theta(fg->verts[i].bezier[j]);
      //printf("converted (%f,%f) to %f\n", fg->verts[i].bezier[j].x, 
      //                                    fg->verts[i].bezier[j].y,
      //                                    ov[i].theta[j]);
    }
  }
  printf("created data structures\n"); fflush(stdout);
  
  
  //we need to find the appropriate target length for the edges
  if (optimize_positions_also == 0) {
    target_length = alpha_scale = 0;
  } else {
    target_length = ((screen_width + screen_height)/2.0)/num_verts;
    alpha_scale = target_length/3.0;
  }  
  
  
  
  objVal = 1e10;
  current_iteration = 0;
  printf("Total dimension: %d\n", total_dimension);
  
  //main loop:
  do {
    prevObjVal = objVal;
    
    //printf("Current iteration: %d\n", current_iteration);
    printf("Current point:\n");
    fatgraph_op_print(ov, num_verts);
    
    
    //get gradient vector
    fatgraph_op_gradient(ov, 
                         grad, 
                         num_verts, 
                         optimize_positions_also,
                         target_length,
                         alpha_scale);
    printf("Gradient:\n");
    fatgraph_op_print(grad, num_verts);
    
    //conjugate gradient part
    if (current_iteration % total_dimension == 0) {
      //regular steepest descent
      //printf("Regular steepest descent iteration\n");
      for (i=0; i<num_verts; i++) {
        dir[i].loc = grad[i].loc;
        for (j=0; j<grad[i].num_edges; j++) {
          dir[i].theta[j] = grad[i].theta[j];
        }
      }
      //calculate the grad norm, as we will need it in the future
      gradNorm = 0;
      for (i=0; i<num_verts; i++) {
        gradNorm += (grad[i].loc.x * grad[i].loc.x) + 
                    (grad[i].loc.y * grad[i].loc.y);
        for (j=0; j<grad[i].num_edges; j++) {
          gradNorm += grad[i].theta[j] * grad[i].theta[j];
        }
      }
    } else {
      //conjugate stuff
      //printf("Conjugate iteration\n");
      gradNorm = 0;
      for (i=0; i<num_verts; i++) {
        gradNorm += (grad[i].loc.x * grad[i].loc.x) + 
                    (grad[i].loc.y * grad[i].loc.y);
        for (j=0; j<grad[i].num_edges; j++) {
          gradNorm += grad[i].theta[j] * grad[i].theta[j];
        }
      }
      //printf("Grad norm: %f\n", gradNorm);
      //printf("Old grad norm: %f\n", oldGradNorm);
      beta = gradNorm / oldGradNorm;
      //printf("beta: %f\n", beta);
      
      //new direction vector
      for (i=0; i<num_verts; i++) {
        dir[i].loc.x = grad[i].loc.x + beta*dir[i].loc.x; 
        dir[i].loc.y = grad[i].loc.y + beta*dir[i].loc.y;
        for (j=0; j<grad[i].num_edges; j++) {
          dir[i].theta[j] = grad[i].theta[j] + beta*dir[i].theta[j];
        }
      }
      
    }
    
    printf("Direction we will take:\n");
    fatgraph_op_print(dir, num_verts);
    
      
    //follow gradient for line minimization
    fatgraph_op_line_max(ov, 
                         dir, 
                         temp, 
                         num_verts,
                         optimize_positions_also,
                         target_length,
                         alpha_scale);
    
    //get the new objective value
    objVal = fatgraph_op_objective_value(ov, 
                                         num_verts, 
                                         optimize_positions_also, 
                                         target_length, 
                                         alpha_scale);
    printf("new objective value: %f\n", objVal);
    printf("old objective value: %f\n", prevObjVal);
    
    //copy the gradient into oldGrad
    oldGradNorm = gradNorm;
    for (i=0; i<num_verts; i++) {
      oldGrad[i].loc = grad[i].loc;
      for (j=0; j<grad[i].num_edges; j++) {
        oldGrad[i].theta[j] = grad[i].theta[j];
      }
    }
    
    current_iteration++;
    
  } while (fabs(objVal - prevObjVal) > tol);
  
  printf("Finished main loop since %f = fabs(%f - %f) <= %f\n", fabs(objVal - prevObjVal), objVal, prevObjVal, tol);
  
  
  
   
  //turn thetas back into vectors and scale
  double distance_to_target;
  double max_x, min_x, max_y, min_y;
  double scale_factor, scale_factor_x, scale_factor_y;
  double offset_x, offset_y;
  
  //we need to make sure that all of our vertices are actually within the 
  //correct range -- find the max x and y differences, then scale to make that 
  //amount fit, then translate so it's in the right place
  if (optimize_positions_also == 1) {
    max_x = min_x = ov[0].loc.x;
    max_y = min_y = ov[0].loc.y;
    for (i=1; i<num_verts; i++) {
      if (ov[i].loc.x < min_x) {
        min_x = ov[i].loc.x;
      } else if (ov[i].loc.x > max_x) {
        max_x = ov[i].loc.x;
      }
      if (ov[i].loc.y < min_y) {
        min_y = ov[i].loc.y;
      } else if (ov[i].loc.y > max_y) {
        max_y = ov[i].loc.y;
      }
    }
    //scale everything
    scale_factor_x = screen_width/(1.3*(max_x-min_x));
    scale_factor_y = screen_height/(1.3*(max_y-min_y));
    scale_factor = (scale_factor_x > scale_factor_y ? 
                    scale_factor_y : 
                    scale_factor_x);
    for (i=0; i<num_verts; i++) {
      ov[i].loc.x *= scale_factor;
      ov[i].loc.y *= scale_factor;
    }
    //printf("min_x: %f\nmax_x: %f\nmin_y: %f\nmax_y: %f\nscale_factor: %f\n", min_x, max_x, min_y, max_y, scale_factor);
    //now we need to move everything so it's in the right window
    offset_y = (min_y*scale_factor) - (screen_height*(3.0/26.0));
    offset_x = (min_x*scale_factor) - (screen_width*(3.0/26.0));
    for (i=0; i<num_verts; i++) {
      ov[i].loc.x -= offset_x;
      ov[i].loc.y -= offset_y;
      //printf("ov vertex %d at: %f, %f\n", i, ov[i].loc.x , ov[i].loc.y);
    }
    //printf("offset_x: %f\noffset_y: %f\n", offset_x, offset_y);
  }
    
  
  
  
  //now build the fatgraph back from the ov
  for (i=0; i<num_verts; i++) {
    for (j=0; j<fg->verts[i].num_edges; j++) {
      fg->verts[i].bezier[j] = theta_to_vector2d(ov[i].theta[j]);
      if (ov[i].dest_verts[j] != i) {
        distance_to_target = sqrt( (ov[ov[i].dest_verts[j]].loc.x - ov[i].loc.x) *
                             (ov[ov[i].dest_verts[j]].loc.x - ov[i].loc.x) +
                             (ov[ov[i].dest_verts[j]].loc.y - ov[i].loc.y) *
                             (ov[ov[i].dest_verts[j]].loc.y - ov[i].loc.y));
      } else { //it's a self loop
        if (optimize_positions_also == 1) {
          distance_to_target = (screen_width+screen_height)/3.0;
        } else {
          distance_to_target = (screen_width+screen_height)/3.0;
        }
        //printf("A bezier of self loop at %d is %f,%f\n", i, fg->verts[i].bezier[j].x, 
        //                                                    fg->verts[i].bezier[j].y);
        //printf("it's: %f, %f\n", fg->verts[i].bezier[j].x, fg->verts[i].bezier[j].y);
        //printf("which I'll multiply by %f\n", distance_to_target/2.0);
        
      }
      fg->verts[i].bezier[j].x *= distance_to_target/2.0;
      fg->verts[i].bezier[j].y *= distance_to_target/2.0;
    }
    fg->verts[i].loc.x = ov[i].loc.x;
    fg->verts[i].loc.y = ov[i].loc.y;
    //printf("vertex %d at: %f, %f\n", i, fg->verts[i].loc.x , fg->verts[i].loc.y);
  }
  
  
  vector2d scaledb1;
  vector2d scaledb2;
  vector2d vecSum;
  //vector2d vecProj;
  double sumDot;
  double mult;
  
  //if a vertex has only two edges and it's not a self loop, then straighten them
  for (i=0; i<num_verts; i++) {
    printf("Vertex %d has %d edges\n", i, fg->verts[i].num_edges);
    if (fg->verts[i].num_edges == 2 && fg->verts[i].edges[0] != fg->verts[i].edges[1]) {
      printf("(vertex %d: I'm straightening (%f,%f) and (%f,%f) to: ", i,fg->verts[i].bezier[0].x,
                                                           fg->verts[i].bezier[0].y,
                                                           fg->verts[i].bezier[1].x,
                                                           fg->verts[i].bezier[1].y);
      scaledb1 = fg->verts[i].bezier[0];
      mult = sqrt(vector2d_dot(scaledb1, scaledb1));
      scaledb1.x /= mult;
      scaledb1.y /= mult;
      scaledb2 = fg->verts[i].bezier[1];
      mult = sqrt(vector2d_dot(scaledb2, scaledb2));
      scaledb2.x /= mult;
      scaledb2.y /= mult;
      vecSum = vector2d_add(scaledb1, scaledb2);
      printf("vector sum direction: (%f, %f)\n", vecSum.x, vecSum.y);
      sumDot = vector2d_dot(vecSum, vecSum);
      if (fabs(sumDot) > 0.01) {
        mult = vector2d_dot(fg->verts[i].bezier[0], vecSum)/sumDot;
        printf("mult1: %f\n", mult);
        fg->verts[i].bezier[0].x -= mult*vecSum.x;
        fg->verts[i].bezier[0].y -= mult*vecSum.y;
        mult = vector2d_dot(fg->verts[i].bezier[1], vecSum)/sumDot;
        printf("mult2: %f\n", mult);
        fg->verts[i].bezier[1].x -= mult*vecSum.x;
        fg->verts[i].bezier[1].y -= mult*vecSum.y;
      }
      printf(" (%f,%f) and (%f,%f)\n", fg->verts[i].bezier[0].x,
                                                           fg->verts[i].bezier[0].y,
                                                           fg->verts[i].bezier[1].x,
                                                           fg->verts[i].bezier[1].y); 
                                                           fflush(stdout);
    }
  }
        
  

  
  //free memory
  for (i=0; i<num_verts; i++) {
    free(ov[i].dest_verts);
    free(ov[i].theta);
    free(grad[i].dest_verts);
    free(grad[i].theta);
    free(oldGrad[i].dest_verts);
    free(oldGrad[i].theta);
    free(temp[i].dest_verts);
    free(temp[i].theta);
  }
  free(ov);
  free(dir);
  free(grad);
  free(oldGrad);
  free(temp);
}

