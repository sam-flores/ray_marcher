

#define DEPTH 2

int quad(double A, double B, double C, double t[2]){


  int sols = -1;
  double det;

  t[0] = -1;
  t[1] = -1;
  det = B*B - 4*A*C;
  if(det > 0){
    sols = 2;
    t[0] = ((-B + sqrt(det))/(2*A));
    t[1] = ((-B - sqrt(det))/(2*A));
  }else if(det < 0){
    sols = 0;
  }else{
    sols = 1;
    t[0] = ((-B)/(2*A));
  }

  return sols;

}

int calc_plane(double A, double B, double t[2]){

  t[0] = -1;
  if(B != 0){
    t[0] = -A/B;
  }

}

double calc_mag(double N[3]){
  return sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
}

void compute_norm(double rxyz[3], int j_min, double N[3]){

    double px = partialX[j_min](rxyz);
    double py = partialY[j_min](rxyz);
    double pz = partialZ[j_min](rxyz);

    N[0] = px*obinv[j_min][0][0] + py*obinv[j_min][1][0] + pz*obinv[j_min][2][0];
    N[1] = px*obinv[j_min][0][1] + py*obinv[j_min][1][1] + pz*obinv[j_min][2][1];
    N[2] = px*obinv[j_min][0][2] + py*obinv[j_min][1][2] + pz*obinv[j_min][2][2];

    double mag = calc_mag(N);
    N[0] /= mag;
    N[1] /= mag;
    N[2] /= mag;

}

double calc_dist(double P[3], double Q[3]){

  double a = P[0] - Q[0];
  double b = P[1] - Q[1];
  double c = P[2] - Q[2]; // c is zero in 2D

  return sqrt(a*a + b*b + c*c);

}

void compute_reflection(double res[3],double src[3],
                          double tip[3], double N[3]){


      double V[3];
      double mag, vDotN;

      V[0] = src[0] - tip[0];
      V[1] = src[1] - tip[1];
      V[2] = src[2] - tip[2];

      // mag = calc_mag(V);
      // V[0] /= mag;
      // V[1] /= mag;
      // V[2] /= mag;

      vDotN = (V[0]*N[0] + V[1]*N[1] + V[2]*N[2]);

      res[0] = 2*vDotN*N[0] - V[0];
      res[1] = 2*vDotN*N[1] - V[1];
      res[2] = 2*vDotN*N[2] - V[2];


}


int shadow_ray(double rSrc[3], double rTip[3]){

  int i, j, j_min, sols = 0;
  double r_xyz[3], T[2], N[3], test[3];
  double xs, ys, zs, xe, ye, ze, dx, dy,dz, A, B, C;
  double temp_rSrc[3], temp_rTip[3];
  double min_t = 1e50, tval;

  // find t for general intersection equation
  // calculate what's below and then pass it
  // back through tranformation matrix


  j_min = -1 ;
  for(j = 0; j < num_objects; j++){ //loop through all objects

    temp_rSrc[0] = rSrc[0];
    temp_rSrc[1] = rSrc[1];
    temp_rSrc[2] = rSrc[2];

    temp_rTip[0] = rTip[0];
    temp_rTip[1] = rTip[1];
    temp_rTip[2] = rTip[2];

    M3d_mat_mult_pt(temp_rSrc, obinv[j], temp_rSrc) ;
    M3d_mat_mult_pt(temp_rTip, obinv[j], temp_rTip) ;

    xs = temp_rSrc[0]; ys = temp_rSrc[1]; zs = temp_rSrc[2];
    xe = temp_rTip[0]; ye = temp_rTip[1]; ze = temp_rTip[2];
    dx = temp_rTip[0] - temp_rSrc[0];
    dy = temp_rTip[1] - temp_rSrc[1];
    dz = temp_rTip[2] - temp_rSrc[2];

    if(type[j] == ELLIPSOID){

        A = dx*dx + dy*dy + dz*dz;
        B = 2*xs*dx + 2*ys*dy + 2*zs*dz;
        C = xs*xs + ys*ys + zs*zs - 1;
        sols = quad(A, B, C, T);

        tval = 1e50 ;
        for (i = 0 ; i < sols ; i++) {
          if (T[i] > 0 && T[i] < tval) { tval = T[i] ; }
        }
      }else if(type[j] == PLANE){
        A = ys;
        B = dy;
        sols = calc_plane(A, B, T);
        tval = 1e50 ;
        test[0] = xs + T[0]*dx;
        test[1] = zs + T[0]*dz;

        if( test[0] >=  domain[j][0]&&
            test[0] <= domain[j][1] &&
            test[1] >= range [j][0] &&
            test[1] <= range [j][1] &&
            T[0] < tval && T[0] > 0) { tval = T[0] ; }

      }else if(type[j] == HYPER){
        A = dx*dx - dy*dy + dz*dz;
        B = 2*xs*dx - 2*ys*dy + 2*zs*dz;
        C = xs*xs - ys*ys + zs*zs - 1;
        sols = quad(A, B, C, T);
        tval = 1e50 ;
        for (i = 0 ; i < sols ; i++) {
            test[0] = xs + T[i]*dx;
            test[1] = ys + T[i]*dy;
            test[2] = zs + T[i]*dz;
            if ( test[0] >= domain[j][0]  &&
                 test[1] >= range[j][0]   &&
                 test[2] >= rangeoid[j][0]&&
                 test[0] <= domain[j][1]  &&
                 test[1] <= range[j][1]   &&
                 test[2] <= rangeoid[j][1]&&
                 T[i] > 0 && T[i] < tval  ) { // not out of range
                 tval = T[i] ;
          }
        }
      }

     if (tval < min_t) {
       min_t = tval ;
       j_min = j ;
     }

  } // end for j

  //------------------------------------------------------------


  if (j_min == -1) {
    // no intersections
    return 0;
  } else {
      return 1;
    }
  }


void ray(double rSrc[3], double rTip[3], double argb[3], int d){

  int i, j, j_min, sols = 0;
  double r_xyz[3], T[2], N[3], test[3];
  double xs, ys, zs, xe, ye, ze, dx, dy,dz, A, B, C;
  double temp_rSrc[3], temp_rTip[3];
  double min_t = 1e50, tval;

  // find t for general intersection equation
  // calculate what's below and then pass it
  // back through tranformation matrix


  j_min = -1 ;
  for(j = 0; j < num_objects; j++){ //loop through all objects

    temp_rSrc[0] = rSrc[0];
    temp_rSrc[1] = rSrc[1];
    temp_rSrc[2] = rSrc[2];

    temp_rTip[0] = rTip[0];
    temp_rTip[1] = rTip[1];
    temp_rTip[2] = rTip[2];

    M3d_mat_mult_pt(temp_rSrc, obinv[j], temp_rSrc) ;
    M3d_mat_mult_pt(temp_rTip, obinv[j], temp_rTip) ;

    xs = temp_rSrc[0]; ys = temp_rSrc[1]; zs = temp_rSrc[2];
    xe = temp_rTip[0]; ye = temp_rTip[1]; ze = temp_rTip[2];
    dx = temp_rTip[0] - temp_rSrc[0];
    dy = temp_rTip[1] - temp_rSrc[1];
    dz = temp_rTip[2] - temp_rSrc[2];

    if(type[j] == ELLIPSOID){

        A = dx*dx + dy*dy + dz*dz;
        B = 2*xs*dx + 2*ys*dy + 2*zs*dz;
        C = xs*xs + ys*ys + zs*zs - 1;
        sols = quad(A, B, C, T);

        tval = 1e50 ;
        for (i = 0 ; i < sols ; i++) {
          if (T[i] > 0 && T[i] < tval) { tval = T[i] ; }
        }
      }else if(type[j] == PLANE){
        A = ys;
        B = dy;
        sols = calc_plane(A, B, T);
        tval = 1e50 ;
        test[0] = xs + T[0]*dx;
        test[1] = zs + T[0]*dz;

        if( test[0] >=  domain[j][0]&&
            test[0] <= domain[j][1] &&
            test[1] >= range [j][0] &&
            test[1] <= range [j][1] &&
            T[0] < tval && T[0] > 0) { tval = T[0] ; }

      }else if(type[j] == HYPER){
        A = dx*dx - dy*dy + dz*dz;
        B = 2*xs*dx - 2*ys*dy + 2*zs*dz;
        C = xs*xs - ys*ys + zs*zs - 1;
        sols = quad(A, B, C, T);
        tval = 1e50 ;
        for (i = 0 ; i < sols ; i++) {
            test[0] = xs + T[i]*dx;
            test[1] = ys + T[i]*dy;
            test[2] = zs + T[i]*dz;
            if ( test[0] >= domain[j][0]  &&
                 test[1] >= range[j][0]   &&
                 test[2] >= rangeoid[j][0]&&
                 test[0] <= domain[j][1]  &&
                 test[1] <= range[j][1]   &&
                 test[2] <= rangeoid[j][1]&&
                 T[i] > 0 && T[i] < tval  ) { // not out of range
                 tval = T[i] ;
          }
        }
      }

     if (tval < min_t) {
       min_t = tval ;
       j_min = j ;

       // point in object space
       r_xyz[0] = temp_rSrc[0] + min_t*(temp_rTip[0] - temp_rSrc[0]) ;
       r_xyz[1] = temp_rSrc[1] + min_t*(temp_rTip[1] - temp_rSrc[1]) ;
       r_xyz[2] = temp_rSrc[2] + min_t*(temp_rTip[2] - temp_rSrc[2]) ;
     }

  } // end for j

  //------------------------------------------------------------


  if (j_min == -1) {
    // no intersections
    // printf("no intersections\n") ;
    argb[0] = 0; argb[1] = 0; argb[2] = 0;

  } else {
    double res[3];
    // printf("intersections\n") ;
    // r_xyz needs object space at this point in code for compute_norm
    compute_norm(r_xyz, j_min, N); // N is the eyespace unit norm
    get_rgb(r_xyz, domain[j_min][0], domain[j_min][1],  // r_xyz in object space
                    range[j_min][0], range[j_min][1],
                    rangeoid[j_min][0], rangeoid[j_min][1],
                    type[j_min], color[j_min]); //color pixel based on u,v position

    // need point in eyespace :
    M3d_mat_mult_pt(r_xyz, obmat[j_min], r_xyz) ;
    compute_reflection(res, rSrc, rTip, N); //in world space
    double nSrc[3], nTip[3];

    nSrc[0] = r_xyz[0] + 0.001*res[0];
    nSrc[1] = r_xyz[1] + 0.001*res[1];
    nSrc[2] = r_xyz[2] + 0.001*res[2];
    nTip[0] = r_xyz[0] + res[0];
    nTip[1] = r_xyz[1] + res[1];
    nTip[2] = r_xyz[2] + res[2];
    double rgb[3];
    Light_Model (color[j_min], rSrc, r_xyz, N, rgb) ;

    if(shadow_ray(nSrc, light_in_eye_space) == 1){
      double f = AMBIENT / (AMBIENT + MAX_DIFFUSE) ;
      rgb[0] = f * color[j_min][0] ;
      rgb[1] = f * color[j_min][1] ;
      rgb[2] = f * color[j_min][2] ;
    } // do shadow ray trace


    if(d == DEPTH){ // base case
      argb[0] = rgb[0]*(1 - ref[j_min]);
      argb[1] = rgb[1]*(1 - ref[j_min]);
      argb[2] = rgb[2]*(1 - ref[j_min]);

    }else if(d < DEPTH){ // recurse case
          double aargb[3];
          d++;
          ray(nSrc, nTip, aargb, d);
          argb[0] = rgb[0]*(1 - ref[j_min]) + aargb[0]*ref[j_min]; // CONFUSING
          argb[1] = rgb[1]*(1 - ref[j_min]) + aargb[1]*ref[j_min];
          argb[2] = rgb[2]*(1 - ref[j_min]) + aargb[2]*ref[j_min];

      }
    }
  }
