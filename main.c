
double local_time = 0;

#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"
#include "light_model.c"

double eye[3], coi[3], up[3] ;
double vm[4][4], vi[4][4];
int width = 600;
int height = 600;
int window_width, window_height, window_square_size ;
double Half_window_size ;
double Half_angle_degrees;
double Tan_half_angle ;
int numobjects = 0;


#define MAXOBJ 10
#define EPSILON .01
#define PARTIAL .01
#define FAR 500
#define SPHERE 0
#define TORUS 1


struct Obj{
  double obmat[4][4];
  double obinv[4][4];
  double irgb[4];
  int type;
  double (*sdf)(double p[3], int id);
};

struct Scene{
  struct Obj obj[MAXOBJ];
};

struct Scene scene;

double torus(double point[3], int id){
  double r = .5; // unit torus
  double R = 1;
  double temp[3];
  M3d_mat_mult_pt(temp, scene.obj[id].obinv, point);
  double len = sqrt(temp[0]*temp[0] + temp[2]*temp[2]) - R;
  return sqrt(len*len + temp[1]*temp[1]) - r;
}

double sphere(double point[3], int id){
  double r = 1; // unit sphere
  double temp[3];
  M3d_mat_mult_pt(temp, scene.obj[id].obinv, point);
  return sqrt(  temp[0]*temp[0]  +
                temp[1]*temp[1]  +
                temp[2]*temp[2]) - r;
}

void set_objmat(double pos[3], double scale[3], double rot[3], int id){

  double A[4][4], Ai[4][4];
  int Tn, Ttypelist[100];
  double Tvlist[100];

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] = scale[0] ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] = scale[1] ;  Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] = scale[2] ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] = rot[0] ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] = rot[1] ;  Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] = rot[2] ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] = pos[0] ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] = pos[1] ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] = pos[2] ; Tn++ ;


  M3d_make_movement_sequence_matrix(A, Ai, Tn,
    Ttypelist, Tvlist); //move object

    //creat objmat and obinv matricies
    M3d_mat_mult(scene.obj[id].obmat, vm, A) ;
    M3d_mat_mult(scene.obj[id].obinv, Ai, vi) ;

}

void set_obj_sdf(int type, int id){

    if(type == TORUS){
      scene.obj[id].sdf = torus;
    }else if(type == SPHERE){
      scene.obj[id].sdf = sphere;
    }

}

void create_shape(int shape, int this_obj, double pos[3], double scale[3],
                    double rot[3], double rgba[4]){

  scene.obj[this_obj].type = shape;
  scene.obj[this_obj].irgb[0] = rgba[0];
  scene.obj[this_obj].irgb[1] = rgba[1];
  scene.obj[this_obj].irgb[2] = rgba[2];
  scene.obj[this_obj].irgb[3] = rgba[3]; // alpha
  set_obj_sdf(shape, this_obj);
  set_objmat(pos, scale, rot, this_obj);

}

void normalize(double v[3]){

    v[0] *= 1e9;
    v[1] *= 1e9;
    v[2] *= 1e9;
    double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      v[0] /= mag;
      v[1] /= mag;
      v[2] /= mag;
}

void compute_norm(double rxyz[3], int id, double N[3]){

    double dplus[3], dminus[3];

    dplus[0] = rxyz[0] + PARTIAL;
    dplus[1] = rxyz[1];
    dplus[2] = rxyz[2];

    dminus[0] = rxyz[0] - PARTIAL;
    dminus[1] = rxyz[1];
    dminus[2] = rxyz[2];
    N[0] = scene.obj[id].sdf(dplus, id) - scene.obj[id].sdf(dminus, id);

    dplus[0] = rxyz[0];
    dplus[1] = rxyz[1] + PARTIAL;
    dplus[2] = rxyz[2];

    dminus[0] = rxyz[0];
    dminus[1] = rxyz[1] - PARTIAL;
    dminus[2] = rxyz[2];
    N[1] = scene.obj[id].sdf(dplus, id) - scene.obj[id].sdf(dminus, id);

    dplus[0] = rxyz[0];
    dplus[1] = rxyz[1];
    dplus[2] = rxyz[2] + PARTIAL;

    dminus[0] = rxyz[0];
    dminus[1] = rxyz[1];
    dminus[2] = rxyz[2] - PARTIAL;
    N[2] = scene.obj[id].sdf(dplus, id) - scene.obj[id].sdf(dminus, id);

}

void ray(double ra[3], double rb[3], double rgb[3]){

    double ra_min[3], temp_ra[3];
    double d, d_min = 10;
    int id_min = -1;

    // find closest object to camera
    // set d_min to d if it's smaller than current d_min
    for(int id = 0; id < numobjects; id++){

      // reset temp_ra
      temp_ra[0] = ra[0];
      temp_ra[1] = ra[1];
      temp_ra[2] = ra[2];

      int step = 0;
      // get initial sdf val

      d = scene.obj[id].sdf(temp_ra, id);

      // if ray is very close to obj or very far from scene, kick out
      // else, advance the ray towards obj a magnitude of sdf(ra)
      while (fabs(d) > EPSILON && step < FAR) {
          //advance the ray
          temp_ra[0] = temp_ra[0] + rb[0] * d;
          temp_ra[1] = temp_ra[1] + rb[1] * d;
          temp_ra[2] = temp_ra[2] + rb[2] * d;

          // get sdf(ra);
          d = scene.obj[id].sdf(temp_ra, id);

          // kick out if our ray is off into oblivion
          if(fabs(d) > 10){
            step = FAR;
            break;
          }

          // else, go on to next step
          step++;
      }// end while loop

      if(step >= FAR){
          // do nothing
      }else if(d_min - d >= .001){
          ra_min[0] = temp_ra[0];
          ra_min[1] = temp_ra[1];
          ra_min[2] = temp_ra[2];
          id_min = id;
          d_min = d;
      }
    }
    // no min obj found, set color to black
    if(id_min == -1){
      rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;

    // else, set color to closest intersected object color
    }else{
      double norm[3];
      // both ra and norm are in eye space
      compute_norm(ra_min, id_min, norm);
      Light_Model(scene.obj[id_min].irgb, eye, ra_min, norm, rgb);
    }
}

void def_camera(){
  eye[0] = 0;
  eye[1] = 0;
  eye[2] = 0;

  coi[0] = eye[0] + 0;
  coi[1] = eye[1] + 0;
  coi[2] = eye[2] + 1;

  up[0] = eye[0] + 0;
  up[1] = eye[1] + 1;
  up[2] = eye[2] + 0;

  M3d_view (vm, vi,  eye, coi, up) ; //create view matrix


  // size of largest square INside window
  if (width < height) { window_square_size = width ; }
                               else { window_square_size = height ; }
  Half_window_size = 0.5*window_square_size ;
  Half_angle_degrees = 30 ;
  Tan_half_angle = tan(Half_angle_degrees*M_PI/180) ;

}

int main(){

  int q = '\0';

  //init graphics
  G_init_graphics(width, height);

  while(q != 'q'){

    //define the camera
    def_camera();

    //define the scene
    numobjects = 0;

    double pos[3], scale[3], rot[3], rgba[4];
    pos[0] = 3 - 3*sin(local_time); pos[1] = 0; pos[2] = 5;
    rot[0] = 360*sin(local_time); rot[1] = 0; rot[2] = 0;
    scale[0] = 1; scale[1] = 1; scale[2] = 1;
    rgba[0] = .4; rgba[1] = .2; rgba[2] = 1; rgba[3] = 1;
    create_shape(TORUS, numobjects, pos, scale, rot, rgba);
    numobjects++;

    pos[0] = -3 + 3*sin(local_time); pos[1] = 0; pos[2] = 5;
    rot[0] = -360*sin(local_time); rot[1] = 0; rot[2] = 0;
    scale[0] = 1; scale[1] = 1; scale[2] = 1;
    rgba[0] = 0; rgba[1] = .1; rgba[2] = 1; rgba[3] = 1;
    create_shape(TORUS, numobjects, pos, scale, rot, rgba);
    numobjects++;

    pos[0] = 0; pos[1] = 0; pos[2] = 5;
    rot[0] = 0; rot[1] = 0; rot[2] = 0;
    scale[0] = 1; scale[1] = cos(local_time); scale[2] = 1;
    rgba[0] = 1; rgba[1] = 0; rgba[2] = 1; rgba[3] = 1;
    create_shape(SPHERE, numobjects, pos, scale, rot, rgba);
    numobjects++;

    // screen space transform variables
    double a,b,c ;
    double H = Tan_half_angle;
    a = Half_window_size / H ;
    b = 0.5*width ;
    c = 0.5*height;
    double rgb[3];

    // define ray vector
    double ra[3], rb[3];

    //shoot rays
    for(double x = 0; x <= width; x++){
      for(double y = 0; y <= height; y++){

        // screen space transformation
        rb[0] = (x - b) / a - eye[0];
        rb[1] = (y - c) / a - eye[1];
        rb[2] = .5 - eye[2];

        // place ray_o @ eye
        ra[0] = eye[0];
        ra[1] = eye[1];
        ra[2] = eye[2];
        // ray march
        if((int)x % 1 == 0 && (int)y % 1 == 0) ray(ra, rb, rgb);
        else{ rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;}
        G_rgb(rgb[0], rgb[1], rgb[2]);
        G_point(x, y);
      }
    }
    // printf("done rendering!\n");

    //display
    G_display_image();

    //check for exit key
    q = G_no_wait_key();

    //increment time
    local_time += M_PI/50;

}

}
