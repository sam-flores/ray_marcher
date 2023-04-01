

#define MAXOBJ 100


double obmat[MAXOBJ][4][4] ;
double obinv[MAXOBJ][4][4] ;
double color[MAXOBJ][3] ;
double ref[MAXOBJ];
int    num_objects ;
int 	 type[MAXOBJ];
double domain[MAXOBJ][2], range[MAXOBJ][2], rangeoid[MAXOBJ][2];
double (*partialX[MAXOBJ])(double xyz[3]);
double (*partialY[MAXOBJ])(double xyz[3]);
double (*partialZ[MAXOBJ])(double xyz[3]);
double vm[4][4], vi[4][4];

#define ELLIPSOID 0
#define HYPER 1
#define PLANE 2


void sphere(double x, double y, double z, double scl){

	double A[4][4], Ai[4][4];
	int Tn, Ttypelist[MAXOBJ];
	double Tvlist[MAXOBJ];

	domain[num_objects][0] = -M_PI; domain[num_objects][1] = M_PI; //defined in world space
	range[num_objects][0] = 0; range[num_objects][1] = M_PI;

	Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] = scl ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] = x ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] = y ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] = z ; Tn++ ;

	M3d_make_movement_sequence_matrix(A, Ai, Tn,
		Ttypelist, Tvlist); //move object

    M3d_mat_mult(obmat[num_objects], vm, A) ;
    M3d_mat_mult(obinv[num_objects], Ai, vi) ;

    type[num_objects] = ELLIPSOID;

    partialX[num_objects] = uC_partialX;
    partialY[num_objects] = uC_partialY;
    partialZ[num_objects] = uC_partialZ;
}

void plane(double x, double y, double z, double scl, double r){

	double A[4][4], Ai[4][4];
	int Tn, Ttypelist[100];
	double Tvlist[100];

	domain[num_objects][0] = -2; domain[num_objects][1] = 2; //defined in world space
	range[num_objects][0] = -2 ; range[num_objects][1] = 2 ;

	Tn = 0 ;
	Ttypelist[Tn] = RX ; Tvlist[Tn] = r ; Tn++ ;

	Ttypelist[Tn] = TX ; Tvlist[Tn] = x ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] = y ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] = z ; Tn++ ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] = scl ; Tn++ ;




	M3d_make_movement_sequence_matrix(A, Ai, Tn,
		Ttypelist, Tvlist); //move object

    M3d_mat_mult(obmat[num_objects], vm, A) ;
    M3d_mat_mult(obinv[num_objects], Ai, vi) ;
    type[num_objects] = PLANE;
    partialX[num_objects] = plane_partialX;
    partialY[num_objects] = plane_partialY;
    partialZ[num_objects] = plane_partialZ;
}

void hloid(double x, double y, double z, double scl){

	double A[4][4], Ai[4][4];
	int Tn, Ttypelist[100];
	double Tvlist[100];

	domain[num_objects][0] = -M_PI; domain[num_objects][1] = M_PI; // theta
	range[num_objects][0] = -1 ; range[num_objects][1] = 1 ; // height
	rangeoid[num_objects][0] = -M_PI/2 ; rangeoid[num_objects][1] = M_PI/2 ; // phi

	Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] = scl ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] = scl ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] = x   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] = y   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] = z   ; Tn++ ;

	M3d_make_movement_sequence_matrix(A, Ai, Tn,
		Ttypelist, Tvlist); //move object

    M3d_mat_mult(obmat[num_objects], vm, A) ;
    M3d_mat_mult(obinv[num_objects], Ai, vi) ;
    type[num_objects] = HYPER;
    partialX[num_objects] = hyper_partialX;
    partialY[num_objects] = hyper_partialY;
    partialZ[num_objects] = hyper_partialZ;
}
