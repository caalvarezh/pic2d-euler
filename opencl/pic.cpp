#define BLOCK_SIZE 1024
#define BLOCK_SIZE2 32
#include <CL/cl.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>
#include <fftw3.h>

using namespace std;
namespace pic_cl {
  int db = 0;
  const int MAX_SPE     = 100000000;           // Limite (computacional) de Superpartículas electrónicas
  //  const int MAX_SPI     = 1000;           // Limite (computacional) de Superpartículas iónicas
  const int J_X         = 2048;           // Número de puntos de malla X. Recomendado: Del orden 2^n+1
  const int J_Y         = 1024;         // Número de puntos de malla Y. Recomendado: Del orden 2^n
  const int ELECTRONS   = 0;
  const int IONS        = 1;
  //  const int X           = 0;
  //  const int Y           = 1;
  const int RAZON_MASAS = 1.98e5;    // M_I/E_MASS (Plata)

  const int FACTOR_CARGA_E = 10;
  const int FACTOR_CARGA_I = 10;       //Número de partículas por superpartícula.
  const int FACT_EL = (-1);

  const double E_MASS      = 9.10938291e-31;  // Masa del Electrón
  const double E_CHARGE    = 1.6021e-19;      // Carga del Electrón
  const double K_BOLTZMANN = 1.3806504e-23;   // Constante de Boltzmann
  const double EPSILON_0   = 8.854187e-12;    // Permitividad eléctrica del vacío
  const double FACT_I = (1. / RAZON_MASAS);
  const double T0 = 1e-13;                   //Escala de tiempo: Tiempo de vaporización
  const double VFLUX_I = 1e3;  // Componentes de Velocidad de flujo iónico
  const double VFLUX_E_X = (std::sqrt(RAZON_MASAS) * VFLUX_I);
  const double VFLUX_E_Y = (std::sqrt(RAZON_MASAS) * VFLUX_I);
  const double VPHI_I_X = (VFLUX_I / VFLUX_I);    // Velocidad térmica Iónica (X)
  const double VPHI_I_Y = (VFLUX_I / VFLUX_I);    // Velocidad térmica Iónica (Y)
  const double VPHI_E_X = (VFLUX_E_X / VFLUX_I);    // Velocidad térmica Electrónica (X)
  const double VPHI_E_Y = (VFLUX_E_Y / VFLUX_I);    // Velocidad térmica Electrónica (Y)
  const double FLUJO_INICIAL = 4.5e33;        // Flujo inicial (# part/m^2*s)
  const double FI_MAXWELL_X = (2. / (M_PI * VPHI_I_X));     // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica (X)
  const double FI_MAXWELL_Y = (1. / (M_PI * VPHI_I_Y));     // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica
  const double FE_MAXWELL_X =  (2. / (M_PI * VPHI_E_X));    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
  const double FE_MAXWELL_Y = (1. / (M_PI * VPHI_E_Y));     // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
  const double M_I = (RAZON_MASAS * E_MASS);                // masa Ión
  const double X0 = (VFLUX_I * T0);                         //Escala de longitud: Distancia recorrida en x por un ión en el tiempo t_0
  const double CTE_E = (RAZON_MASAS * X0 / (VFLUX_I * T0));
  const double CTE_ER = RAZON_MASAS;
  const double DT = 1.e-5;                                  // Paso temporal
  const double VFLUX_I_magnitud =  std::sqrt(VFLUX_I * VFLUX_I + VFLUX_I * VFLUX_I); // Velocidad de flujo iónico (m/s)  =  sqrt(2*K_BOLTZMANN*Te/(M_PI*M_I))
  const double vflux_e_magnitud =  std::sqrt(VFLUX_E_X * VFLUX_E_X + VFLUX_E_Y * VFLUX_E_Y);
  const double Te = (M_PI * 0.5 * E_MASS * (std::pow(VFLUX_E_X, 2) / K_BOLTZMANN));    // Temperatura electrónica inicial (°K)

  const double NI03D = (FLUJO_INICIAL / VFLUX_I);
  const double NE03D = (FLUJO_INICIAL / VFLUX_E_X);
  const double LAMBDA_D = sqrt(EPSILON_0 * K_BOLTZMANN * Te / (NE03D * pow(E_CHARGE, 2)));  //Longitud de Debye
  const double DELTA_X = (LAMBDA_D);   //Paso espacial
  const double L_MAX_X = (((J_X-1) * DELTA_X) / X0);                      // Longitud región de simulación
  const double L_MAX_Y = (((J_Y-1) * DELTA_X) / X0);                      // Longitud región de simulación

  const double cte_rho = pow(E_CHARGE * T0, 2) / (M_I * EPSILON_0 * pow(X0, 3)); //Normalización de EPSILON_0
  const int    NTe = 1e5;
  const int    NTI = 1e5;                                  //Número de partículas "reales"
  const double n_0 = double(NTe);                   // Densidad de partículas

  inline void checkErr(cl_int err, const char * name) {
    if (err != CL_SUCCESS) {
      std::cerr << "ERROR: " << name  << " (" << err << ")" << std::endl;
      exit(EXIT_FAILURE);
    }
  }



  //*********************
  //Velocidades Iniciales
  //*********************

  double create_Velocities_X(double fmax,double vphi) {// función para generar distribución semi-maxwelliana de velocidades de las particulas
    // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
    double sigma = vphi;                           // sigma = vflujo = vth    ( "dispersión" de la distribución Maxweliana)
    double vmin =  0. ;                            // Rapidez mínima
    double vmax =  4. * sigma;                       // Rapidez máxima
    double v, f, f_random;

    static int flag  = 1;
    if (flag  ==  0) {
      int seed  =  time (NULL);
      srand (seed);
      flag  =  1;
    }

    v = vmin+(vmax-vmin)*double(rand())/double(RAND_MAX); // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
    f  = fmax*exp(-(1.0/M_PI)*pow(v/vphi,2));

    f_random  =  fmax*double(rand())/double(RAND_MAX);    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

    if (f_random > f)
      return create_Velocities_X(fmax,vphi);
    else
      return  v;
  }

  double create_Velocities_Y(double fmax,double vphi) {// función para generar distribución semi-maxwelliana de velocidades de las particulas
    // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
    double sigma = vphi;                           // sigma = vflujo = vth    ( "dispersión" de la distribución Maxweliana)
    double vmin =  -3.*sigma;                            // Rapidez mínima
    double vmax =  3.*sigma;                       // Rapidez máxima
    double v,f,f_random;

    static int flag  =  1;
    if (flag  ==  0) {
      int seed  =  time (NULL);
      srand (seed);
      flag  =  1;
    }

    v = vmin + (vmax - vmin) * double(rand())/ double(RAND_MAX); // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
    f = fmax * exp(-(1.0 / M_PI) * pow(v / vphi, 2));

    f_random  =  fmax * double(rand()) / double(RAND_MAX);    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

    if (f_random > f) return create_Velocities_Y(fmax, vphi);
    else return  v;
  }

  void initialize_Particles (double *pos_x, double *pos_y, double *vel_x, double *vel_y,
      int NSP, int fmax_x, int fmax_y, int vphi_x, int vphi_y) {
    for (int i = 0; i < MAX_SPE; i++) {
      pos_x[i + NSP] = 0;
      vel_x[i + NSP] = create_Velocities_X (fmax_x, vphi_x);
      pos_y[i + NSP] = L_MAX_Y / 2.0;
      vel_y[i + NSP] = create_Velocities_Y(fmax_y, vphi_y);
    }
  }

  //**************************************************************************************
  //Determinación del aporte de carga de cada superpartícula sobre las 4 celdas adyacentes
  //**************************************************************************************

  void Concentration (double *pos_x, double *pos_y, double *n, int NSP, double hx) {
    int j_x,j_y;
    double temp_x,temp_y;
    double jr_x,jr_y;
    for(int i = 0; i < J_X * J_Y; i++) {
      n[i] = 0.;
    } // Inicializar densidad de carga

    for (int i = 0; i < NSP;i++) {
      jr_x = pos_x[i] / hx; // indice (real) de la posición de la superpartícula
      j_x  = (int) jr_x;    // indice  inferior (entero) de la celda que contiene a la superpartícula
      temp_x  =  jr_x - j_x;
      jr_y = pos_y[i] / hx; // indice (real) de la posición de la superpartícula
      j_y  = (int) jr_y;    // indice  inferior (entero) de la celda que contiene a la superpartícula
      temp_y  =  jr_y - j_y;

      n[j_y + j_x * J_Y] += (1. - temp_x) * (1. - temp_y) / (hx * hx * hx);
      n[j_y + (j_x + 1) * J_Y] += temp_x * (1. - temp_y) / (hx * hx * hx);
      n[(j_y + 1) + j_x * J_Y] += (1. - temp_x) * temp_y / (hx * hx * hx);
      n[(j_y + 1) + (j_x + 1) * J_Y] += temp_x * temp_y / (hx * hx * hx);
    }
  }

  //***********************************************************************
  //Cálculo del Potencial Electrostático en cada uno de los puntos de malla
  //***********************************************************************
  //

  void poisson2D_dirichletX_periodicY(double *phi, complex<double> *rho, double hx) {
    int M = J_X - 2, N = J_Y;
    double h = hx;
    double hy = hx;
    double *f;
    fftw_complex  *f2;
    fftw_plan p,p_y,p_i,p_yi;
    f= (double*) fftw_malloc(sizeof(double)* M);
    f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_y = fftw_plan_dft_1d(N, f2, f2, FFTW_FORWARD, FFTW_ESTIMATE);
    p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_yi = fftw_plan_dft_1d(N, f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Columnas FFT
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < M; j++)
        f[j] = rho[(j + 1) * N + k].real();
      fftw_execute(p);
      for (int j = 0; j < M; j++) {
        rho[(j + 1) * N + k].real(f[j]);
      }
    }

    // Filas FFT
    for (int j = 0; j < M; j++) {
      for (int k = 0; k < N; k++)
        memcpy( &f2[k], &rho[(j + 1) * N + k], sizeof( fftw_complex ) );
      fftw_execute(p_y);
      for (int k = 0; k < N; k++)
        memcpy( &rho[(j + 1) * N + k], &f2[k], sizeof( fftw_complex ) );
    }

    // Resolver en el espacio de Fourier
    complex<double> i(0.0, 1.0);
    double pi = M_PI;
    complex<double> Wy = exp(2.0 * pi * i / double(N));
    complex<double> Wn = 1.0;
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        complex<double> denom = h * h * 2.0 + hy * hy * 2.0;
        denom -= hy * hy * (2 * cos((m + 1) * pi / (M + 1))) + h * h *(Wn + 1.0 / Wn);
        if (denom != 0.0)
          rho[(m + 1) * N + n] *= h * h * hy * hy / denom;
        Wn *= Wy;
      }
    }

    // Inversa de las filas
    for (int j = 0; j < M; j++) {
      for (int k = 0; k < N; k++)
        memcpy( &f2[k], &rho[(j + 1) * N + k], sizeof( fftw_complex ) );
      fftw_execute(p_yi);
      for (int k = 0; k < N; k++) {
        memcpy( &rho[(j + 1) * N + k], &f2[k], sizeof( fftw_complex ) );
        rho[(j + 1) * N + k] /= double(N); //La transformada debe ser normalizada.
      }
    }


    //Inversa C lumnas FFT
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < M; j++)
        f[j]=rho[(j + 1) * N + k].real();
      fftw_execute(p_i);
      for (int j = 0; j < M; j++)
        phi[(j + 1) * N + k] = f[j] / double(2 * (M + 1));
    }

    for (int k = 0; k < N; k++) {
      phi[0 * N + k]=0;
      phi[(J_X - 1) * N + k]=0;
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(p_i);
    fftw_destroy_plan(p_y);
    fftw_destroy_plan(p_yi);
    fftw_free(f); fftw_free(f2);

  }

  void electric_field(double *phi, double *E_X, double *E_Y, double hx) {

    for (int j = 1; j < J_X - 1; j++) {
      for (int k = 0; k < J_Y; k++) {
        E_X[j * J_Y + k] = (phi[(j - 1) * J_Y + k]
            - phi[(j + 1) * J_Y + k]) / (2. * hx);
        E_Y[j * J_Y + k] = (phi[j * J_Y + ((J_Y + k - 1) % J_Y)]
            - phi[j * J_Y + ((k + 1) % J_Y)]) / (2. * hx);

        E_X[k] = 0.0;  //Cero en las fronteras X
        E_Y[k] = 0.0;
        E_X[(J_X - 1) * J_Y + k] = 0.0;
        E_Y[(J_X - 1) * J_Y + k] = 0.0;
      }
    }
  }

  void H_electric_field (double *h_phi, double *h_E_X, double *h_E_Y, double hx,
      cl::Context context, cl::Program program, cl::CommandQueue queue) {
    int size = J_X * J_Y * sizeof(double);
    //prepare memory
    // Create the memory buffers
    cl_int error;
    cl::Buffer d_phi = cl::Buffer(context, CL_MEM_READ_ONLY, size);
    cl::Buffer d_E_X = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    cl::Buffer d_E_Y = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    // Copy the input data to the input buffers using the
    // command queue for the first device
    error = queue.enqueueWriteBuffer(d_phi, CL_TRUE, 0, size, h_phi);
    checkErr(error, "queue 1");
    error = queue.enqueueWriteBuffer(d_E_X, CL_TRUE, 0, size, h_E_X);
    checkErr(error, "queue 2");
    error = queue.enqueueWriteBuffer(d_E_Y, CL_TRUE, 0, size, h_E_Y);
    checkErr(error, "queue 3");
    //launch kernel
    // Make kernel
    cl::Kernel D_electric_field(program, "D_electric_field", &error);
    checkErr(error, "kernel");
    // Set the kernel arguments
    error = D_electric_field.setArg(0, d_phi);
    checkErr(error, "set arg 1");
    error = D_electric_field.setArg(1, d_E_X);
    checkErr(error, "set arg 2");
    error = D_electric_field.setArg(2, d_E_Y);
    checkErr(error, "set arg 3");
    error = D_electric_field.setArg(3, hx);
    checkErr(error, "set arg 4");
    error = D_electric_field.setArg(4, J_X);
    checkErr(error, "set arg 5");
    error = D_electric_field.setArg(5, J_Y);
    checkErr(error, "set arg 6");
    // Execute the kernel
    cl::NDRange global(J_X, J_Y);
    //checkErr(error, "global");
    error = queue.enqueueNDRangeKernel(D_electric_field, cl::NullRange, global, cl::NullRange);
    checkErr(error, "call kernel");

/*
    cl::Kernel D_electric_field_border(program, "D_electric_field_border");
    // Set the kernel arguments
    error = D_electric_field_border.setArg(0, d_E_X);
    checkErr(error, "set arg 1");
    error = D_electric_field_border.setArg(1, d_E_Y);
    checkErr(error, "set arg 2");
    error = D_electric_field_border.setArg(2, hx);
    checkErr(error, "set arg 3");
    error = D_electric_field_border.setArg(3, J_X);
    checkErr(error, "set arg 4");
    error = D_electric_field_border.setArg(4, J_Y);
    checkErr(error, "set arg 4");

    cl::NDRange global1(BLOCK_SIZE2, 1, 1);
    cl::NDRange local1(ceil(double(J_Y) / BLOCK_SIZE2), 1, 1);
    error = queue.enqueueNDRangeKernel(D_electric_field_border, cl::NullRange, global1, local1);
    checkErr(error, "call kernel");
*/
    // get the answer and free memory
    // Copy the output data back to the host
    error = queue.enqueueReadBuffer(d_E_X, CL_TRUE, 0, size, h_E_X);
    checkErr(error, "read buffer");
    error = queue.enqueueReadBuffer(d_E_Y, CL_TRUE, 0, size, h_E_Y);
    checkErr(error, "read buffer2");
  }

void Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int &NSP, int especie,
     double *E_X, double *E_Y, double hx, int &total_perdidos, double &mv2perdidas) {
    int j_x,j_y;
    double temp_x,temp_y,Ep_X, Ep_Y,fact;
    double jr_x,jr_y;
    int kk1 = 0;

    if (especie ==  ELECTRONS)
      fact = FACT_EL * CTE_E * FACTOR_CARGA_E;
    else
      fact = FACT_I * CTE_E * FACTOR_CARGA_E;

    for (int i = 0; i < NSP; i++) {
      jr_x = pos_x[i] / hx;     // Índice (real) de la posición de la superpartícula (X)
      j_x  = int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
      temp_x = jr_x - double(j_x);
      jr_y = pos_y[i] / hx;     // Índice (real) de la posición de la superpartícula (Y)
      j_y  = int(jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
      temp_y  =  jr_y-double(j_y);

      Ep_X = (1 - temp_x) * (1 - temp_y) * E_X[j_x * J_Y + j_y] +
        temp_x * (1 - temp_y) * E_X[(j_x + 1) * J_Y + j_y] +
        (1 - temp_x) * temp_y * E_X[j_x * J_Y + (j_y + 1)] +
        temp_x * temp_y * E_X[(j_x + 1) * J_Y + (j_y + 1)];

      Ep_Y = (1 - temp_x) * (1 - temp_y) * E_Y[j_x * J_Y + j_y] +
        temp_x * (1 - temp_y) * E_Y[(j_x + 1) * J_Y + j_y] +
        (1 - temp_x) * temp_y * E_Y[j_x * J_Y + (j_y + 1)] +
        temp_x * temp_y * E_Y[(j_x + 1) * J_Y + (j_y + 1)];


      vel_x[i] += (DT * fact) * Ep_X;
      vel_y[i] += (DT * fact) * Ep_Y;

      pos_x[i] += vel_x[i] * DT;
      pos_y[i] += vel_y[i] * DT;

      if(pos_x[i] < 0) {//Rebote en la pared del material.
        pos_x[i] = -pos_x[i];
        vel_x[i] = -vel_x[i];
      }

      while(pos_y[i] > L_MAX_Y) //Ciclo en el eje Y.
        pos_y[i] = pos_y[i] - L_MAX_Y;

      while(pos_y[i]<0.0) //Ciclo en el eje Y.
        pos_y[i] = L_MAX_Y + pos_y[i];

      if(pos_y[i] < 0.0) //Ciclo en el eje Y.
        pos_y[i] += L_MAX_Y;

      if(pos_x[i] >= 0 && pos_x[i] <= L_MAX_X) {
        pos_x[kk1] = pos_x[i];
        pos_y[kk1] = pos_y[i];
        vel_x[kk1] = vel_x[i];
        vel_y[kk1] = vel_y[i];
        kk1++;
      }

      //Salida espacio de Fase
    }

  }


  void H_Motion(double *h_pos_x, double *h_pos_y, double *h_vel_x, double *h_vel_y, int &NSP,
      int especie, double *h_E_X, double *h_E_Y, double hx, int &total_perdidos, double &mv2perdidas,
      cl::Context &context, cl::Program &program, cl::CommandQueue &queue) {
    double fact;
    int k = 0;
    if (especie ==  ELECTRONS)
      fact = FACT_EL * CTE_E * FACTOR_CARGA_E;
    else
      fact = FACT_I * CTE_E * FACTOR_CARGA_E;
    //LANZAR KERNEL

    int size  = NSP * sizeof(double);
    int size1 = J_X * J_Y * sizeof(double);
    //prepare memory
    // Create the memory buffers
    cl::Buffer d_pos_x = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    cl::Buffer d_pos_y = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    cl::Buffer d_vel_x = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    cl::Buffer d_vel_y = cl::Buffer(context, CL_MEM_READ_WRITE, size);
    cl::Buffer d_E_X = cl::Buffer(context, CL_MEM_READ_ONLY, size1);
    cl::Buffer d_E_Y = cl::Buffer(context, CL_MEM_READ_ONLY, size1);
    // Copy the input data to the input buffers using the
    // command queue for the first device
    queue.enqueueWriteBuffer(d_pos_x, CL_TRUE, 0, size, h_pos_x);
    queue.enqueueWriteBuffer(d_pos_y, CL_TRUE, 0, size, h_pos_y);
    queue.enqueueWriteBuffer(d_vel_x, CL_TRUE, 0, size, h_vel_x);
    queue.enqueueWriteBuffer(d_vel_y, CL_TRUE, 0, size, h_vel_y);
    queue.enqueueWriteBuffer(d_E_X, CL_TRUE, 0, size1, h_E_X);
    queue.enqueueWriteBuffer(d_E_Y, CL_TRUE, 0, size1, h_E_Y);

    //launch kernel
    // Make kernel
    cl_int err;
    cl::Kernel D_Motion(program, "D_Motion", &err);
    checkErr(err, "kernel execution");
    // Set the kernel arguments
    D_Motion.setArg(0, d_pos_x);
    D_Motion.setArg(1, d_pos_y);
    D_Motion.setArg(2, d_vel_x);
    D_Motion.setArg(3, d_vel_y);
    D_Motion.setArg(4, NSP);
    D_Motion.setArg(5, fact);
    D_Motion.setArg(6, d_E_X);
    D_Motion.setArg(7, d_E_Y);
    D_Motion.setArg(8, hx);
    D_Motion.setArg(9, L_MAX_X);
    D_Motion.setArg(10, L_MAX_Y);
    D_Motion.setArg(11, DT);
    D_Motion.setArg(12, J_X);
    D_Motion.setArg(13, J_Y);
    // Execute the kernel
    cl::NDRange global(NSP);
    //cl::NDRange local(ceil(float(NSP) / BLOCK_SIZE));
    queue.enqueueNDRangeKernel(D_Motion, cl::NullRange, global, cl::NullRange);

    // get the answer and free memory
    queue.enqueueReadBuffer(d_pos_x, CL_TRUE, 0, size, h_pos_x);
    queue.enqueueReadBuffer(d_pos_y, CL_TRUE, 0, size, h_pos_y);
    queue.enqueueReadBuffer(d_vel_x, CL_TRUE, 0, size, h_vel_x);
    queue.enqueueReadBuffer(d_vel_y, CL_TRUE, 0, size, h_vel_y);

    // FIN
    for (int i = 0; i < NSP; i++) {
      if (h_pos_x[i] >= 0 && h_pos_x[i] <= L_MAX_X) {
        h_pos_x[k] = h_pos_x[i];
        h_pos_y[k] = h_pos_y[i];
        h_vel_x[k] = h_vel_x[i];
        h_vel_y[k] = h_vel_y[i];
        k++;
      }
    }

  }

}
