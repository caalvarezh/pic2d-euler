#include "pic_cuda.hpp"
#define BLOCK_SIZE 1024


using namespace std;
namespace pic_cuda {
  void prueba(int d){
    cout<<"LINE: "<< d <<endl;
  }

  __device__
    double atomicAdd(double* address, double val) {
      unsigned long long int* address_as_ull =
        (unsigned long long int*)address;
      unsigned long long int old = *address_as_ull, assumed;
      do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val +
              __longlong_as_double(assumed)));
      } while (assumed != old);
      return __longlong_as_double(old);
    }

  void initialize_Particles (double *pos_x, double *pos_y, double *vel_x, double *vel_y,
      int NSP, int fmax_x, int fmax_y, int vphi_x, int vphi_y) {
    for (int i = 0; i < MAX_SPE;i++) {
      pos_x[i + NSP] = 0;
      vel_x[i + NSP] = create_Velocities_X (fmax_x, vphi_x);
      pos_y[i + NSP] = L_MAX_Y / 2.0;
      vel_y[i + NSP] = create_Velocities_Y(fmax_y, vphi_y);
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

    static int flag  =  0;
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

    static int flag  =  0;
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

  //**************************************************************************************
  //Determinación del aporte de carga de cada superpartícula sobre las 4 celdas adyacentes
  //**************************************************************************************
   void H_Concentration (double *h_pos_x, double *h_pos_y, double *h_n, int NSP, double hx) {
    int size  = NSP * sizeof(double);
    int size1 = J_X * J_Y * sizeof(double);
    double *d_pos_x , *d_pos_y, *d_n;
    cudaMalloc(&d_pos_x, size);
    cudaMalloc(&d_pos_y, size);
    cudaMalloc(&d_n, size1);
    cudaMemcpy(d_pos_x, h_pos_x, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_y, h_pos_y, size, cudaMemcpyHostToDevice);
    cudaMemset(d_n, 0, size1);

    dim3 dimBlock(BLOCK_SIZE,1,1);
    dim3 dimGrid(ceil(NSP/float(BLOCK_SIZE)),1, 1);
    D_Concentration<<< dimGrid, dimBlock >>>(d_pos_x, d_pos_y, d_n, NSP, hx);
    cudaDeviceSynchronize();
    cudaMemcpy(h_n, d_n, size1, cudaMemcpyDeviceToHost);
    cudaFree(d_pos_x);
    cudaFree(d_pos_y);
    cudaFree(d_n);
  }

  __global__
    void D_Concentration (double *d_pos_x, double *d_pos_y, double *d_n, int NSP,double hx) {
      int i = blockIdx.x*blockDim.x+threadIdx.x;
      int j_x, j_y;
      double temp_x, temp_y, tmp;
      double jr_x, jr_y;
      if(i < NSP) {
        jr_x = d_pos_x[i] / hx; // indice (real) de la posición de la superpartícula
        j_x  = (int) jr_x;    // indice  inferior (entero) de la celda que contiene a la superpartícula
        temp_x  =  jr_x - j_x;
        jr_y = d_pos_y[i] / hx; // indice (real) de la posición de la superpartícula
        j_y  = (int) jr_y;    // indice  inferior (entero) de la celda que contiene a la superpartícula
        temp_y  =  jr_y - j_y;
        __threadfence();
        tmp = (1. - temp_x) * (1. - temp_y) / (hx * hx * hx);
        atomicAdd(&d_n[j_y + j_x * J_Y], tmp);
        tmp = temp_x * (1. - temp_y) / (hx * hx * hx);
        atomicAdd(&d_n[j_y + (j_x + 1) * J_Y], tmp);
        tmp = (1. - temp_x) * temp_y / (hx * hx * hx);
        atomicAdd(&d_n[(j_y + 1) + j_x * J_Y], tmp);
        tmp = temp_x * temp_y / (hx * hx * hx);
        atomicAdd(&d_n[(j_y + 1) + (j_x + 1) * J_Y], tmp);
      }
    }

  __global__
  void D_rhoKernel(double *ne, double *ni, double *rho_d, double cte_rho) {
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    int i = blockIdx.y*blockDim.y+threadIdx.y;
    int index;
    index = (i*J_Y+j);
    if((i < J_X) && (j < J_Y)) {
      rho_d[index] = cte_rho * FACTOR_CARGA_E * (ni[index] - ne[index]) / n_0;
    }
  }

  void H_rhoKernel (double *ne, double *ni, double *rho_h) {
    int size1 = J_X * J_Y * sizeof(double);
    double *rho_d;
    cudaMalloc(&rho_d, size1);
    cudaMemcpy(rho_d, rho_h, size1, cudaMemcpyHostToDevice);
    dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE,1);
    dim3 dimGrid(ceil(J_Y / (double)BLOCK_SIZE) ,ceil(J_Y / (double) BLOCK_SIZE), 1);
    D_rhoKernel<<< dimGrid, dimBlock >>>(ne, ni, rho_d, cte_rho);
    cudaDeviceSynchronize();
    cudaMemcpy(rho_h,rho_d, size1, cudaMemcpyDeviceToHost);
    cudaFree(rho_d);
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
    //fftw_complex  *f2;

    //f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_plan p, p_i;

    f = (double*) fftw_malloc(sizeof(double)* N);
    complex<double> *buff_rho = (complex<double> *) malloc( sizeof(complex<double>) * M);

    p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);


    // Columnas FFT
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < M; j++)
        f[j] = rho[(j + 1) * N + k].real();
      fftw_execute(p);
      for (int j = 0; j < M; j++)
        rho[(j + 1) * N + k].real() = f[j];
    }


    //   printVectorToFile(f,M);

    // Filas FFT
    cufftDoubleComplex *f2_d;
    cufftHandle p_y, p_yi;
    cudaMalloc((void**)&f2_d, sizeof(cufftDoubleComplex) * N);
    int size2 = sizeof(cufftDoubleComplex) * N;
    cufftPlan1d(&p_y, N,CUFFT_Z2Z,1);
    cufftPlan1d(&p_yi,N,CUFFT_Z2Z,1);

    for (int j = 0; j < M; j++) {
      for (int k = 0; k < N; k++) {
        buff_rho[k] = rho[(j + 1) * N + k];
      }

      cudaMemcpy(&f2_d, &buff_rho, size2, cudaMemcpyHostToDevice);
      cufftExecZ2Z(p_y,f2_d,f2_d,CUFFT_FORWARD);
      cudaMemcpy(&buff_rho, &f2_d, size2, cudaMemcpyDeviceToHost);

      for (int k = 0; k < N; k++) {
        rho[(j + 1) * N + k] = buff_rho[k];
      }
    }

    //   printComplexVectorToFile(f2, N);

    // Resolver en el espacio de Fourier
    complex<double> i(0.0, 1.0);
    double pi = M_PI;
    complex<double> Wy = exp(2.0 * pi * i / double(N));
    complex<double> Wn = 1.0;
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        complex<double> denom = h * h * 2.0 + hy * hy * 2.0;
        denom -= hy * hy * (2 * cos((m + 1) * pi / (M + 1))) + h * h * (Wn + 1.0 / Wn);
        if (denom != 0.0)
          rho[(m + 1) * N + n] *= h * h * hy * hy / denom;
        Wn *= Wy;
      }
    }

    // Inversa de las filas
    for (int j = 0; j < M; j++) {
      for (int k = 0; k < N; k++)
        buff_rho[k] = rho[(j + 1) * N + k];
      cudaMemcpy(&f2_d, &buff_rho, size2, cudaMemcpyHostToDevice);
      //fftw_execute(p_yi);
      cufftExecZ2Z(p_yi,f2_d,f2_d,CUFFT_INVERSE);
      cudaMemcpy(&buff_rho, &f2_d, size2, cudaMemcpyDeviceToHost);
      for (int k = 0; k < N; k++) {
        rho[(j + 1) * N + k] = buff_rho[k] / double(N); //La transformada debe ser normalizada.
      }
    }

    //Inversa Columnas FFT
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < M; j++)
        f[j]=rho[(j + 1) * N + k].real();
      fftw_execute(p_i);
      for (int j = 0; j < M; j++)
        phi[(j + 1) * N + k] = f[j] / double(2 * (M + 1));
    } // pasar a cuda..

    //Aqui es importante pasar el phi de la función para el potencial

    for (int k = 0; k < N; k++) {
      phi[k]=0;
      phi[(J_X - 1) * N + k]=0;
    }

    fftw_destroy_plan(p);
    cufftDestroy(p_y);
    cufftDestroy(p_yi);
    fftw_destroy_plan(p_i);

    fftw_free(f);
    cudaFree(f2_d);
  }

  //*********************************************************
  __global__
    void D_electric_field (double *d_phi, double *d_E_X, double *d_E_Y, double hx) {
      int j = (blockIdx.x * blockDim.x + threadIdx.x)+1;
      if(j < J_X - 1) {
        for (int k = 1; k < J_Y - 1; k++) {
          d_E_X[j * J_Y + k] = (d_phi[(j - 1) * J_Y + k] - d_phi[(j + 1) * J_Y + k]) / (2. * hx);
          d_E_Y[j * J_Y + k] = (d_phi[j * J_Y + (k - 1)] - d_phi[j * J_Y + (k + 1)]) / (2. * hx);
          d_E_X[k] = 0.0; //Cero en las fronteras X
          d_E_Y[k] = 0.0;
          d_E_X[(J_X - 1) * J_Y + k] = 0.0;
          d_E_Y[(J_X - 1) * J_Y + k] = 0.0;
        }
        d_E_X[j * J_Y + 0] = (d_phi[(j - 1) * J_Y] - d_phi[((j + 1) * J_Y)]) / (2. * hx);
        d_E_Y[j * J_Y + 0] = (d_phi[j * J_Y + (J_Y - 1)] - d_phi[j * J_Y + 1]) / (2. * hx);
        d_E_X[j * J_Y + (J_Y - 1)] = (d_phi[(j - 1) * J_Y + (J_Y - 1)] - d_phi[(j + 1) * J_Y + (J_Y - 1)]) / (2. * hx);
        d_E_Y[j * J_Y + (J_Y - 1)] = (d_phi[j * J_Y + (J_Y - 2)] - d_phi[j * J_Y]) / (2. * hx);
      }
    }

  void H_electric_field (double *h_phi, double *h_E_X, double *h_E_Y, double hx) {
    int size1 = J_X * J_Y * sizeof(double);
    double *d_E_X, *d_E_Y, *d_phi;

    cudaMalloc(&d_E_X, size1);
    cudaMalloc(&d_E_Y, size1);
    cudaMalloc(&d_phi, size1);
    cudaMemcpy(d_phi, h_phi, size1, cudaMemcpyHostToDevice);
    cudaMemset(d_E_X, 0, size1);
    cudaMemset(d_E_Y, 0, size1);

    dim3 dimBlock(BLOCK_SIZE, 1, 1);
    dim3 dimGrid3(ceil(float(J_X) / BLOCK_SIZE), 1, 1);

    D_electric_field<<< dimGrid3, dimBlock >>> (d_phi, d_E_X, d_E_Y, hx);
    cudaDeviceSynchronize();

    cudaMemcpy(h_E_X,d_E_X, size1, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_E_Y,d_E_Y, size1, cudaMemcpyDeviceToHost);
    cudaFree(d_E_X);
    cudaFree(d_E_Y);
    cudaFree(d_phi);
  }

  __global__
    void D_Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, bool *visit,
        int NSP, int fact, double *E_X, double *E_Y, double hx, double L_MAX_X, double L_MAX_Y) {
      int j_x,j_y;
      double temp_x,temp_y,Ep_X, Ep_Y;
      double jr_x,jr_y;
      int i = blockIdx.x * blockDim.x + threadIdx.x;
      if ( i < NSP) {
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


        vel_x[i] += CTE_E * FACTOR_CARGA_E * fact * Ep_X * DT;
        vel_y[i] += CTE_E * FACTOR_CARGA_E * fact * Ep_Y * DT;

        pos_x[i] += vel_x[i] * DT;
        pos_y[i] += vel_y[i] * DT;

        if(pos_x[i] < 0) {//Rebote en la pared del material.
          pos_x[i] *= -1;
          vel_x[i] *= -1;
        }

        if (pos_x[i] >= L_MAX_X) {//Partícula fuera del espacio de Simulación
          visit[i] = false;
        }

        pos_y[i] = fmod(pos_y[i], L_MAX_Y);

        if(pos_y[i] < 0.0) //Ciclo en el eje Y.
          pos_y[i] += L_MAX_Y;

      }
    }

  void H_Motion(double *h_pos_x, double *h_pos_y, double *h_vel_x, double *h_vel_y, int &NSP,
      int especie, double *h_E_X, double *h_E_Y, double hx, int &total_perdidos, double &mv2perdidas) {
    double fact;
    int k = 0;
    int conteo_perdidas = 0;

    if (especie ==  ELECTRONS)
      fact = FACT_EL;
    else
      fact = FACT_I;

    //LANZAR KERNEL

    int size  = NSP * sizeof(double);
    int size1 = J_X * J_Y * sizeof(double);
    double *d_pos_x, *d_pos_y, *d_vel_x, *d_vel_y, *d_E_X, *d_E_Y;
    bool *d_visit, *h_visit;
    h_visit = (bool *) malloc((size / 2));
    cudaMalloc(&d_pos_x, size);
    cudaMalloc(&d_pos_y, size);
    cudaMalloc(&d_vel_x, size);
    cudaMalloc(&d_vel_y, size);
    cudaMalloc(&d_visit, size / 2);
    cudaMalloc(&d_E_X, size1);
    cudaMalloc(&d_E_Y, size1);
    cudaMemcpy(d_pos_x, h_pos_x, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_y, h_pos_y, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel_x, h_vel_x, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel_y, h_vel_y, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_E_X, h_E_X, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(d_E_Y, h_E_Y, size1, cudaMemcpyHostToDevice);
    cudaMemset(d_visit, true, size / 2);

    dim3 dimBlock(BLOCK_SIZE,1,1);
    dim3 dimGrid(ceil(NSP/float(BLOCK_SIZE)),1, 1);
    D_Motion<<< dimGrid, dimBlock >>>(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_visit, NSP, fact, d_E_X, d_E_Y, hx, L_MAX_X, L_MAX_Y);
    cudaDeviceSynchronize();
    cudaMemcpy(h_pos_x, d_pos_x, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_pos_y, d_pos_y, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel_x, d_vel_x, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel_y, d_vel_y, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_visit, d_visit, size / 2, cudaMemcpyDeviceToHost);
    cudaFree(d_pos_x);
    cudaFree(d_pos_y);
    cudaFree(d_vel_x);
    cudaFree(d_vel_y);
    cudaFree(d_E_X);
    cudaFree(d_E_Y);
    cudaFree(d_visit);


    // FIN
    for (int i = 0; i < NSP; i++) {
      if (!h_visit[i]) {//Partícula fuera del espacio de Simulación
        conteo_perdidas++;
        double tmp = pow( sqrt(h_vel_x[i] * h_vel_x[i] + h_vel_y[i] * h_vel_y[i]) , 2);
        if(especie  ==  ELECTRONS) {
          mv2perdidas+= tmp;
        } else {
          mv2perdidas+= tmp / (RAZON_MASAS);
        }
      }

      if (h_pos_x[i] >= 0 && h_pos_x[i] <= L_MAX_X) {
        h_pos_x[k] = h_pos_x[i];
        h_pos_y[k] = h_pos_y[i];
        h_vel_x[k] = h_vel_x[i];
        h_vel_y[k] = h_vel_y[i];
        k++;
      }
    }
    free(h_visit);
    //Salida espacio de Fase
    total_perdidos += conteo_perdidas;
    NSP -= conteo_perdidas;
  }

  void Funcion_Distribucion(double *pos, double *vel, int NSP, char *archivo_X, char *archivo_Y) {
    double Nc = 100;
    FILE *pFile[2];
    pFile[0]  =  fopen(archivo_X,"w");
    pFile[1]  =  fopen(archivo_Y,"w");
    int suma = 0;
    int ind = 0;
    double a;

    for(int i = 0; i < 2 ;i++) {//Max. & min. velocity values.
      double max = 0;
      double min = 0;
      double dv;
      for (int h = 0; h < NSP; h++) {
        if(vel[h + (i * NSP)] < min)
          min = vel[h + (i * NSP)];//Min. Velocity.

        if(vel[h + (i * NSP)] > max)
          max = vel[h + (i * NSP)];//Max. Velocity.
      }

      dv  =  (max - min) / Nc;
      a = min;//Start velocity counter.

      //printf("min = %e max = %e dv =  %e kt = %d #Particulas  =  %d ", min,max, dv,kt, NSP);
      for(ind = 0; ind < Nc; ind++) {
        suma  = 0;
        for (int j = 0; j < NSP; j++) {
          if(a <=  vel[j + (i * NSP)] && vel[j + (i * NSP)] < a + dv)
            suma++;
        }
        fprintf(pFile[i]," %e  %d  \n", a, suma);
        a  =  a + dv;
      }
      fclose(pFile[i]);
    }
  }

  /*
     void Funcion_Distribucion(double pos[MAX_SPE][2], double vel[MAX_SPE][2] , int NSP, char *archivo_X, char *archivo_Y) {
     double Nc = 100;
     FILE *pFile[2];
//pFile[0]  =  fopen(archivo_X,"w"); pFile[1]  =  fopen(archivo_Y,"w");
int suma = 0;
int ind = 0;
double a;

for(int i = 0;i<2;i++) {//Max. & min. velocity values.
double max = 0;
double min = 0;
double dv;
for (int h = 0;h<NSP;h++) {
if(vel[h][i]<min)
min = vel[h][i];//Min. Velocity.

if(vel[h][i]>max)
max = vel[h][i];//Max. Velocity.

}

dv  =  (max-min)/Nc;
a = min;//Start velocity counter.

//printf("min = %e max = %e dv =  %e kt = %d #Particulas  =  %d ", min,max, dv,kt, NSP);
for(ind = 0; ind < Nc;ind++) {
suma  = 0;
for (int j = 0;j<NSP;j++) {
if(a <=  vel[j][i] && vel[j][i] < a + dv)
suma++;
}
//fprintf(pFile[i]," %e  %d  \n", a, suma);
a  =  a + dv;
}
//fclose(pFile[i]);
}
}*/

}
