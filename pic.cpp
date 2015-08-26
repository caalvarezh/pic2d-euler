#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>
#include <fftw3.h>

#define __d {prueba (__LINE__);}

//Magnitudes físicas

using namespace std;

void prueba(int d){
  cout<<"LINE: "<< d <<endl;
}

void initialize_Particles (double *pos_e, double *vel_e, double *pos_i, double *vel_i
int main() {
  //************************
  // Parámetros del sistema
  //************************

  int le = 0, li = 0;
  double  t_0, x_0;
  int  total_e_perdidos = 0;
  int  total_i_perdidos = 0;
  double  mv2perdidas = 0;

  double  ND = NE03D * pow(LAMBDA_D,3);                          //Parámetro del plasma
  int     NTe = 1e5, NTI = 1e5;                                  //Número de partículas "reales"
  int     k_MAX_inj;   //Tiempo máximo de inyección
  int     K_total;     //Tiempo total
  int     Ntv = 8;
  int     NTSPe, NTSPI, MAX_SPE_dt, MAX_SPI_dt;
  double  cte_rho = pow(E_CHARGE * T0, 2) / (M_I * EPSILON_0 * pow(X0, 3)); //Normalización de EPSILON_0
  double  phi0 = 2. * K_BOLTZMANN * Te / (M_PI * E_CHARGE ), E0 = phi0 / X0;
  FILE    *outEnergia;


  //***************************
  //Constantes de normalización
  //***************************

  //double  X0 = LAMBDA_D;                //Escala de longitud: Longitud de Debye
  double  n0  =  double(NTe) / (X0 * X0 * X0);
  double  ni0_3D  =  NI03D * pow(X0, 3);
  double  ne0_3D  =  NE03D * pow(X0, 3);
  double  om_p  =  VFLUX_E_X / LAMBDA_D;                    //Frecuencia del plasma
  double hx;
  int seed  =  time (NULL);
  srand (seed);  // Semilla para generar números aleatorios dependiendo del reloj interno.

  //******************
  //ARCHIVOS DE SALIDA
  //******************
  outEnergia = fopen("Energia","w");
  char buffer[40], bufferb[40], bufferc[40], bufferd[40];
  /*
     for(int i  =  0; i<= 80; i++) {
     sprintf(buffer, "fase_ele%d", i);
     outFase_ele[i] = fopen(buffer, "w");
     }

     for(int i  =  0; i<= 80; i++) {
     sprintf(buffer, "fase_ion%d", i);
     outFase_ion[i] = fopen(buffer, "w");
     }
     */

  //printf("NI03D = %e \nNE03D = %e \nTemperatura electronica = %e eV \n", NI03D, NE03D, Te * K_BOLTZMANN / (1.602e-19));
  //printf("Longitud de Debye = %e  \nFrecuencia del plasma = %e \n", LAMBDA_D,om_p);
  //printf("Tp = %e \nND = %e \nLX = %e \nLY = %e \n", 2 * M_PI / om_p, ND, Lmax[0], Lmax[1]);

  //printf("CTE_E = %e  \ncte_rho = %e  \nTe  =  %e  \nhx*Ld  =  %e  \n",CTE_E,cte_rho, Te, hx*LAMBDA_D );

  //printf("dt/T0 = %e    \ndt/T = %e   \nhx/LAMBDA_D = %e \nTiempo vapor. = %d dt \n",dt/t_0,dt/T,hx/(LAMBDA_D/X0), k_MAX_inj);

  //****************************************
  // Inicialización de variables del sistema
  //****************************************
  int size = MAX_SPE * 2 * sizeof(double);
  int size1 = MAX_SPI * 2 * sizeof(double);
  int size2 = J_X * J_Y * sizeof(double);
  int size3 = J_X * J_Y * sizeof(complex<double>);

  double *pos_e, *pos_i, *vel_e, *vel_i, *ne, *ni;
  double *phi, *E_X, *E_Y;
  complex<double> *rho;
  double  E_i,E_e,E_field,E_total,E_perdida;

  pos_e = (double *) malloc(size);
  pos_i = (double *) malloc(size1);
  vel_e = (double *) malloc(size);
  vel_i = (double *) malloc(size1);
  ne    = (double *) malloc(size2);
  ni    = (double *) malloc(size2);
  phi   = (double *) malloc(size2);
  E_X   = (double *) malloc(size2);
  E_Y   = (double *) malloc(size2);
  rho   = (complex<double> *) malloc(size3);

  //***************************
  // Normalización de variables
  //***************************

  t_0 = 1;
  x_0 = 1;
  hx = DELTA_X / X0;                            // Paso espacial
  double n_0 = n0 * X0 * X0 * X0;                   // Densidad de partículas
  NTSPe = NTe / FACTOR_CARGA_E;
  NTSPI = NTI / FACTOR_CARGA_I; // Número total de superpartículas
  // Inyectadas en un tiempo T0.
  // ( =  número de superpartículas
  // Inyectadas por unidad de tiempo,
  // puesto que T0*(normalizado) = 1.


  //printf("X0^3 = %e \nn0i = %e \nlambda/hx = %e \nTemp  =  %e\n", X0*X0*X0, NI03D, LAMBDA_D/X0,K_BOLTZMANN*Te);
  //printf("dt = %e \nMAX_SPE_dt = %d  \n",dt_emision/t_0,MAX_SPI_dt);
  //printf("Energia = %e \n",ET0);

  int Kemision = 20;  //Pasos para liberar partículas
  double dt_emision = Kemision * DT; //Tiempo para liberar partículas

  MAX_SPE_dt = NTSPe * dt_emision;   //Número de Superpartículas el. liberadas cada vez.
  MAX_SPI_dt = MAX_SPE_dt;


  // Ciclo de tiempo

  k_MAX_inj = t_0 / DT;
  K_total = Ntv * k_MAX_inj;

  initialize_Particles (pos_e, vel_e, pos_i, vel_i, li, le);//Velocidades y posiciones iniciales de las partículas (no las libera).
  clock_t tiempo0  =  clock();

  for(int kk  =  0, kt  =  0; kt <= K_total; kt++) {
    if(kt % 50000 == 0) {
      printf("kt = %d\n", kt);
      printf("le = %d   li = %d \n",le, li );
    }

    if(kt <= k_MAX_inj && kt == kk) {// Inyectar superpartículas (i-e)
      le+= MAX_SPE_dt;
      li+= MAX_SPI_dt;
      kk = kk + Kemision;
    }
    //-----------------------------------------------
    // Calculo de "densidad de carga 2D del plasma"

    Concentration (pos_e, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    Concentration (pos_i, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    poisson2D_dirichletX_periodicY(phi, rho, hx);
    // Calcular campo eléctrico en puntos de malla

    electric_field(phi, E_X, E_Y, hx);

    // imprimir el potencial electroestatico.
    if(kt % 50000  ==  0) {
      sprintf(buffer,"Poisson%d.data", kt);
      ofstream dataFile(buffer);
      for (int j  =  0; j < J_X; j++) {
        double thisx  =  j * hx;
        for (int k  =  0; k < J_Y; k++) {
          double thisy  =  k * hx;
          dataFile << thisx << '\t' << thisy << '\t' << phi[(j * J_Y) + k] << '\n';
        }
        dataFile << '\n';
      }
      dataFile.close();
    }

    /* //imprimit la densidad
       if(kt % 50000 == 0) {
    // Escribir a archivo
    sprintf(buffer,"n%d.data", kt);
    ofstream dataFile(buffer);
    for (int j  =  0; j < J_X; j++) {
    double thisx  =  j * hx;
    for (int k  =  0; k < J_Y; k++) {
    double thisy  =  k * hx;
    dataFile << thisx << '\t' << thisy << '\t' << ni[j][k] << '\t'<< ne[j][k] << '\t' << E_X[j][k]<< '\t' << E_Y[j][k] <<'\n';
    }
    dataFile << '\n';
    }
    dataFile.close();
    }
    */
    // Avanzar posiciones de superpartículas electrónicas e Iónicas

    Motion(pos_e, vel_e, le, ELECTRONS, E_X, E_Y, kt, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_i, vel_i, li, IONS, E_X, E_Y, kt, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    //Motion_e(pos_e,vel_e,le, E_X, E_Y, total_e_perdidos, mv2perdidas);
    //Motion_i(pos_i,vel_i,li, E_X, E_Y, total_i_perdidos, mv2perdidas);

    //Cálculo de energías.
    if(kt % 2000 == 0 && kt > 0) {
      E_i = 0; //Se inicializan las variables de la energía
      E_e = 0;
      E_field = 0;
      E_total = 0;
      E_perdida  =  0;

      for (int y = 0; y < li; y++)
        E_i  =  E_i + pow(sqrt(vel_i[y] * vel_i[y] + vel_i[y + MAX_SPI] * vel_i[y + MAX_SPI]), 2) / (M_PI);

      for (int y = 0; y < le; y++)
        E_e  =  E_e + pow(sqrt(vel_e[y] * vel_e[y] + vel_e[y + MAX_SPE] * vel_e[y + MAX_SPE]), 2) / (M_PI * RAZON_MASAS);

      //Energía del campo
      for (int r = 0; r < J_X; r++)
        for(int ry = 0; ry < J_Y; ry++)
          E_field  =  E_field +  (ni[(r * J_Y) + ry] - ne[(r * J_Y) + ry]) * phi[(r * J_Y) + ry] * hx * hx * hx;

      E_field  =  2 * E_field / (M_PI);

      E_total  =  E_field + E_i + E_e;

      //Energia perdida por partícula pérdida
      E_perdida  =   mv2perdidas / M_PI;

      fprintf(outEnergia,"%e %e %e %e %e %d  %e \n", kt * DT, E_total, E_i, E_e, E_field, total_e_perdidos + total_i_perdidos, E_perdida );
    }//Cierre de calculo de energia

    clock_t tiempo1  =  clock();
    if(kt % 50000 == 0) {
      cout << " CPU time " << kt / 50000 << "  =  " << double(tiempo1 - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;
      tiempo0  =  clock();
    }

    //Salida de función de distribución

/*    if(kt % 50000  ==  0) {
      sprintf(buffer,"fdist_ele%dx.data", kt);
      sprintf(bufferb,"fdist_ele%dy.data", kt);
      sprintf(bufferc,"fdist_ion%dx.data", kt);
      sprintf(bufferd,"fdist_ion%dy.data", kt);
      Funcion_Distribucion(pos_e,vel_e,le, buffer,  bufferb);
      Funcion_Distribucion(pos_i,vel_i,li, bufferc, bufferd);
    }
*/
  } //Cierre del ciclo principal

  fclose(outEnergia);
  /*for(int i = 0;i<= 80;i++) {
    fclose(outFase_ele[i]);
    fclose(outFase_ion[i]);
  }*/

  free(pos_e);
  free(pos_i);
  free(vel_e);
  free(vel_i);
  free(ne);
  free(ni);
  free(phi);
  free(E_X);
  free(E_Y);
  free(rho);

  return (0);
}// FINAL MAIN

,
    int li, int le) {
  for (int i = 0;i<MAX_SPE;i++) {
    pos_e[i + le] = 0;
    vel_e[i + le] =  create_Velocities_X (FE_MAXWELL_X, VPHI_E_X);
    pos_e[i + le + MAX_SPE] = L_MAX_Y / 2.0;
    vel_e[i + le + MAX_SPE] = create_Velocities_Y(FE_MAXWELL_Y, VPHI_E_Y);

    pos_i[i + li] = 0;
    vel_i[i + li] = create_Velocities_X (FI_MAXWELL_X, VPHI_I_X);
    pos_i[i + li + MAX_SPI] = L_MAX_Y / 2.0;
    vel_i[i + li + MAX_SPI] = create_Velocities_Y (FI_MAXWELL_Y, VPHI_I_Y);
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
void Concentration (double *pos, double *n, int NSP, double hx) {
  int j_x,j_y;
  double temp_x,temp_y;
  double jr_x,jr_y;
  for(int i = 0; i < J_X * J_Y; i++) {
    n[i] = 0.;
  } // Inicializar densidad de carga

  for (int i = 0; i < NSP;i++) {
    jr_x = pos[i] / hx; // indice (real) de la posición de la superpartícula
    j_x  = (int) jr_x;    // indice  inferior (entero) de la celda que contiene a la superpartícula
    temp_x  =  jr_x - j_x;
    jr_y = pos[i + MAX_SPE] / hx; // indice (real) de la posición de la superpartícula
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
  fftw_plan p, p_y, p_i, p_yi;
  f = (double*) fftw_malloc(sizeof(double)* M);
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
    for (int j = 0; j < M; j++)
      rho[(j + 1) * N + k].real() = f[j];
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
      denom -= hy * hy * (2 * cos((m + 1) * pi / (M + 1))) + h * h * (Wn + 1.0 / Wn);
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

  //Inversa Columnas FFT
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < M; j++)
      f[j]=rho[(j + 1) * N + k].real();
    fftw_execute(p_i);
    for (int j = 0; j < M; j++)
      phi[(j + 1) * N + k] = f[j] / double(2 * (M + 1));
  }

  for (int k = 0; k < N; k++) {
    phi[k]=0;
    phi[(J_X - 1) * N + k]=0;
  }

  fftw_destroy_plan(p);
  fftw_destroy_plan(p_i);
  fftw_destroy_plan(p_y);
  fftw_destroy_plan(p_yi);
  fftw_free(f); fftw_free(f2);
}

//*********************************************************
void electric_field(double *phi, double *E_X, double *E_Y, double hx) {

  for (int j = 1; j < J_X - 1; j++) {
    for (int k = 1; k < J_Y - 1; k++) {
      E_X[j * J_Y + k] = (phi[(j - 1) * J_Y + k] - phi[(j + 1) * J_Y + k]) / (2. * hx);
      E_Y[j * J_Y + k] = (phi[j * J_Y + (k - 1)] - phi[j * J_Y + (k + 1)]) / (2. * hx);

      E_X[k] = 0.0;  //Cero en las fronteras X
      E_Y[k] = 0.0;
      E_X[(J_X - 1) * J_Y + k] = 0.0;
      E_Y[(J_X - 1) * J_Y + k] = 0.0;
    }

    E_X[j * J_Y] = (phi[(j - 1) * J_Y] - phi[((j + 1) * J_Y + 0)]) / (2. * hx);
    E_Y[j * J_Y] = (phi[j * J_Y + (J_Y - 1)] - phi[j * J_Y + 1]) / (2. * hx);

    E_X[j * J_Y + (J_Y - 1)] = (phi[(j - 1) * J_Y + (J_Y - 1)] - phi[(j + 1) * J_Y + (J_Y - 1)]) / (2. * hx);
    E_Y[j * J_Y + (J_Y-1)] = (phi[j * J_Y + (J_Y - 2)] - phi[j * J_Y]) / (2. * hx);
  }

}

//*******************************************************

void Motion(double *pos, double *vel, int &NSP, int especie, double *E_X,
    double *E_Y, int kt, double hx, int &total_perdidos, double &mv2perdidas) {
  int j_x,j_y;
  double temp_x,temp_y,Ep_X, Ep_Y,fact;
  double jr_x,jr_y;
  int kk1 = 0;
  int conteo_perdidas = 0;

  if(especie ==  ELECTRONS)
    fact = FACT_EL;
  else
    fact = FACT_I;

  for (int i = 0; i < NSP; i++) {
    jr_x = pos[i] / hx;     // Índice (real) de la posición de la superpartícula (X)
    j_x  = int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
    temp_x = jr_x - double(j_x);
    jr_y = pos[i + MAX_SPE] / hx;     // Índice (real) de la posición de la superpartícula (Y)
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


    vel[i] = vel[i] + CTE_E * FACTOR_CARGA_E * fact * Ep_X * DT;
    vel[i + MAX_SPE] = vel[i + MAX_SPE] + CTE_E * FACTOR_CARGA_E * fact * Ep_Y * DT;

    pos[i] += vel[i] * DT;
    pos[i + MAX_SPE] += vel[i + MAX_SPE] * DT;

    if(pos[i]<0) {//Rebote en la pared del material.
      pos[i] = -pos[i];
      vel[i] = -vel[i];
    }

    if (pos[i] >= L_MAX_X) {//Partícula fuera del espacio de Simulación
      conteo_perdidas++;
      total_perdidos++;
      if(especie  ==  ELECTRONS) {
        //printf("Electron perdido No. %d,  i = %d, kt = %d \n",total_perdidos, i ,kt);
        mv2perdidas+= pow( sqrt(vel[i] * vel[i] + vel[i + MAX_SPE] * vel[i + MAX_SPE]) , 2);
      }
      else {
        //printf("Ion perdido No. %d,  i = %d, kt = %d \n",total_perdidos, i ,kt);
        mv2perdidas+= pow( sqrt(vel[i] * vel[i] + vel[i + MAX_SPE] * vel[i + MAX_SPE]), 2) / (RAZON_MASAS);
      }
    }

    pos[i + MAX_SPE] = fmod(pos[i + MAX_SPE], L_MAX_Y);

    if(pos[i + MAX_SPE] < 0.0) //Ciclo en el eje Y.
      pos[i + MAX_SPE] += L_MAX_Y;

    if(pos[i] >= 0 && pos[i] <= L_MAX_X) {
      pos[kk1] = pos[i];
      pos[kk1 + MAX_SPE] = pos[i + MAX_SPE];
      vel[kk1] = vel[i];
      vel[kk1 + MAX_SPE] = vel[i + MAX_SPE];
      kk1++;
    }

    //Salida espacio de Fase
 }

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
