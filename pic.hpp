#ifndef PICLEAP
#define PICLEAP

namespace pic {



  #define __d {prueba (__LINE__);}
  void prueba(int d);

  void   initialize_Particles(double *pos_e, double *vel_e, double *pos_i,
      double *vel_i, int li, int le);

  double create_Velocities_X(double fmax, double vphi);

  double create_Velocities_Y(double fmax, double vphi);

  void   Concentration(double *pos_x, double *pos_y, double *n, int NSP, double hx);

  void   poisson2D_dirichletX_periodicY(double *phi, std::complex <double> *rho,
      double hx);

  void   electric_field(double *phi, double *E_X, double *E_Y, double hx);

  void Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int &NSP, int especie,
     double *E_X, double *E_Y, double hx, int &total_perdidos, double &mv2perdidas );

  void   Funcion_Distribucion(double *pos, double *vel, int NSP, char *archivo_X,
      char *archivo_Y);
}
#endif
