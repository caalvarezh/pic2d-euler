#ifndef PICUDA
#define PICUDA



namespace pic_cuda {
  using namespace std;
#define __d {prueba (__LINE__);}

  void prueba(int d);

  __device__ double atomicAdd(double* address, double val);

  void initialize_Particles (double *pos_x, double *pos_y, double *vel_x, double *vel_y,
      int NSP, int fmax_x, int fmax_y, int vphi_x, int vphi_y);

  double create_Velocities_X(double fmax, double vphi);

  double create_Velocities_Y(double fmax, double vphi);

  void   H_Concentration(double *h_pos_x, double *h_pos_y, double *h_n, int NSP, double hx);

  __global__ void   D_Concentration(double *d_pos_x, double *d_pos_y, double *d_n, int NSP, double hx);

  void H_rhoKernel (double *ne, double *ni, complex<double> *rho_h);

   __global__ void D_rhoKernel(double *ne, double *ni, complex<double> *rho_d, double cte_rho);

  void   poisson2D_dirichletX_periodicY(double *phi, std::complex<double> *rho, double hx);

  void   H_electric_field (double *h_phi, double *h_E_X, double *h_E_Y, double hx);

  __global__ void D_electric_field (double *d_phi, double *d_E_X, double *d_E_Y, double hx);

  __global__ void D_Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y,
      int &NSP, int fact, double *E_X, double *E_Y, double hx, double L_MAX_X, double L_MAX_Y);

  void H_Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int &NSP, int especie,
      double *E_X, double *E_Y, double hx, int &total_perdidos, double &mv2perdidas);

  void   Funcion_Distribucion(double *pos, double *vel, int NSP, char *archivo_X,
      char *archivo_Y);

}
#endif
