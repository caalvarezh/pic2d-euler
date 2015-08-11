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

#define E_MASS       9.10938291e-31  // Masa del Electrón
#define E_CHARGE     1.6021e-19      // Carga del Electrón
#define K_BOLTZMANN  1.3806504e-23   // Constante de Boltzmann
#define EPSILON_0    8.854187e-12    // Permitividad eléctrica del vacío

#define MAX_SPE      10000           // Limite (computacional) de Superpartículas electrónicas
#define MAX_SPI      10000           // Limite (computacional) de Superpartículas iónicas
//#define J_X         65          cambiar?
#define J_X          64           // Número de puntos de malla X. Recomendado: Del orden 2^n+1
//#define J_Y         16          cambiar?
#define J_Y          16           // Número de puntos de malla Y. Recomendado: Del orden 2^n
#define ELECTRONS    0
#define IONS         1
//#define X           0           cambiar?
#define X            0
//#define Y           1
#define Y            1
#define RAZON_MASAS  1.98e5    // M_I/E_MASS (Plata)

#define FACTOR_CARGA_E 10
#define FACTOR_CARGA_I 10       //Número de partículas por superpartícula.
#define FACT_EL (-1)
#define FACT_I (1. / RAZON_MASAS)
#define T0  1e-13                   //Escala de tiempo: Tiempo de vaporización
#define VFLUX_I 1e3  // Componentes de Velocidad de flujo iónico
#define VFLUX_E_X (sqrt(RAZON_MASAS) * VFLUX_I)
#define VFLUX_E_Y (sqrt(RAZON_MASAS) * VFLUX_I)
#define VPHI_I_X (VFLUX_I / VFLUX_I)    // Velocidad térmica Iónica (X)
#define VPHI_I_Y (VFLUX_I / VFLUX_I)    // Velocidad térmica Iónica (Y)
#define VPHI_E_X (VFLUX_E_X / VFLUX_I)    // Velocidad térmica Electrónica (X)
#define VPHI_E_Y (VFLUX_E_Y / VFLUX_I)    // Velocidad térmica Electrónica (Y)
#define FLUJO_INICIAL 4.5e33        // Flujo inicial (# part/m^2*s)
#define FI_MAXWELL_X (2. / (M_PI * VPHI_I_X))    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica (X)
#define FI_MAXWELL_Y (1. / (M_PI * VPHI_I_Y))    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica
#define FE_MAXWELL_X (2. / (M_PI * VPHI_E_X))    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
#define FE_MAXWELL_Y (1. / (M_PI * VPHI_E_Y))    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
#define M_I (RAZON_MASAS * E_MASS)    // masa Ión
#define X0 (VFLUX_I * T0)   //Escala de longitud: Distancia recorrida en x por un ión en el tiempo t_0
#define CTE_E (RAZON_MASAS * X0 / (VFLUX_I * T0))
#define DT 1.e-5                                 // Paso temporal
#define VFLUX_I_magnitud sqrt(VFLUX_I * VFLUX_I + VFLUX_I * VFLUX_I) // Velocidad de flujo iónico (m/s)  =  sqrt(2*K_BOLTZMANN*Te/(M_PI*M_I))
#define vflux_e_magnitud sqrt(VFLUX_E_X * VFLUX_E_X + VFLUX_E_Y * VFLUX_E_Y)
#define Te (M_PI * 0.5 * E_MASS * pow(VFLUX_E_X, 2) / K_BOLTZMANN)    // Temperatura electrónica inicial (°K)

#define NI03D  (FLUJO_INICIAL / VFLUX_I)
#define NE03D  (FLUJO_INICIAL / VFLUX_E_X)
#define LAMBDA_D sqrt(EPSILON_0 * K_BOLTZMANN * Te / (NE03D * pow(E_CHARGE, 2)))  //Longitud de Debye
#define DELTA_X (LAMBDA_D)   //Paso espacial
#define L_MAX_X (((J_X-1) * DELTA_X) / X0)                      // Longitud región de simulación
#define L_MAX_Y (((J_Y-1) * DELTA_X) / X0)                      // Longitud región de simulación


// La mayoria son constantes
//double  RAZON_MASAS = 1.98e5;     // M_I/E_MASS (Plata)
// (vflujo = sqrt(2K_BOLTZMANNTe/(pi*me))

using namespace std;
void prueba(int d){
          cout<<"LINE: "<< d <<endl;
}
/*
void   initialize_Particles(double *pos_e, double *vel_e, double *pos_i,
                            double *vel_i, int li, int le);

void   Concentration(double *pos, double *n, int NSP, double hx);

void   poisson2D_dirichletX_periodicY(double *phi, complex <double> *rho,
                                      double hx);

void   electric_field(double *phi, double *E_X, double *E_Y, double hx);

void   Motion(double *pos, double *vel, int &NSP, int especie, double *E_X,
              double *E_Y, int kt, double hx, int &total_perdidos,
              double mv2perdidas);

void   Funcion_Distribucion(double *pos, double *vel, int NSP, char *archivo_X,
                            char *archivo_Y);

*/
void   initialize_Particles(double pos_e[MAX_SPE][2], double vel_e[MAX_SPE][2],
                            double pos_i[MAX_SPI][2], double vel_i[MAX_SPI][2],
                            int li, int le);

double create_Velocities_X(double fmax, double vphi);

double create_Velocities_Y(double fmax, double vphi);

void   Concentration(double pos[MAX_SPE][2], double n[J_X][J_Y], int NSP, double hx);

void   poisson2D_dirichletX_periodicY(double phi[J_X][J_Y],
                                      complex <double> rho[J_X][J_Y], double hx);

void   electric_field(double phi[J_X][J_Y], double E_X[J_X][J_Y],
                      double E_Y[J_X][J_Y], double hx);

void  Motion(double pos[MAX_SPE][2],  double vel[MAX_SPE][2],  int &NSP,
             int especie, double E_X[J_X][J_Y], double E_Y[J_X][J_Y], int kt, double hx,
             int &total_perdidos, double &mv2perdidas);

  //void   Motion(double pos[MAX_SPE][2],  double vel[MAX_SPE][2],  int NSP,
//              int especie, double E_X[J_X][J_Y], double E_Y[J_X][J_Y], int kt, double hx);//,
//int total_e_perdidos, int total_i_perdidos, double mv2perdidas);

/*void   Motion_e(double pos[MAX_SPE][2],  double vel[MAX_SPE][2],
                int NSP, double E_X[J_X][J_Y], double E_Y[J_X][J_Y], int kt);//,
//int total_e_perdidos, double mv2perdidas);

void   Motion_i(double pos[MAX_SPE][2],  double vel[MAX_SPE][2],
                int NSP, double E_X[J_X][J_Y], double E_Y[J_X][J_Y],
                int total_i_perdidos, double mv2perdidas, int kt);
*/
void   Funcion_Distribucion(double pos[MAX_SPE][2], double vel[MAX_SPE][2],
                            int NSP, char *archivo_X, char *archivo_Y);



//************************
//Parámetros de simulación
//************************



FILE    *outFase_ele[81];
FILE    *outFase_ion[81];


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

  for(int i  =  0; i<= 80; i++) {
    sprintf(buffer, "fase_ele%d", i);
    outFase_ele[i] = fopen(buffer, "w");
  }

  for(int i  =  0; i<= 80; i++) {
    sprintf(buffer, "fase_ion%d", i);
    outFase_ion[i] = fopen(buffer, "w");
  }


  //printf("NI03D = %e \nNE03D = %e \nTemperatura electronica = %e eV \n", NI03D, NE03D, Te * K_BOLTZMANN / (1.602e-19));
  //printf("Longitud de Debye = %e  \nFrecuencia del plasma = %e \n", LAMBDA_D,om_p);
  //printf("Tp = %e \nND = %e \nLX = %e \nLY = %e \n", 2 * M_PI / om_p, ND, Lmax[0], Lmax[1]);

  printf("CTE_E = %e  \ncte_rho = %e  \nTe  =  %e  \nhx*Ld  =  %e  \n",CTE_E,cte_rho, Te, hx*LAMBDA_D );

  //printf("dt/T0 = %e    \ndt/T = %e   \nhx/LAMBDA_D = %e \nTiempo vapor. = %d dt \n",dt/t_0,dt/T,hx/(LAMBDA_D/X0), k_MAX_inj);

  //****************************************
  // Inicialización de variables del sistema
  //****************************************

  double  pos_e[MAX_SPE][2],pos_i[MAX_SPI][2];  //Vectores de Posición, 0:X 1:Y.
  double  vel_e[MAX_SPE][2], vel_i[MAX_SPI][2]; //Vectores de Velocidad, 0:X 1:Y.
  double  ne[J_X][J_Y],ni[J_X][J_Y];            //Densidad de partículas en puntos de malla.
  complex <double> rho[J_X][J_Y];               //Densidad de carga en los puntos de malla.
  double  phi[J_X][J_Y], E_X[J_X][J_Y], E_Y[J_X][J_Y];
  double  E_i,E_e,E_field,E_total,E_perdida;


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


  printf("X0^3 = %e \nn0i = %e \nlambda/hx = %e \nTemp  =  %e\n", X0*X0*X0, NI03D, LAMBDA_D/X0,K_BOLTZMANN*Te);
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
    if(kt % 10000 == 0) {
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
        rho[j][k] =  cte_rho * FACTOR_CARGA_E * (ni[j][k]- ne[j][k]) / n_0;

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
          dataFile << thisx << '\t' << thisy << '\t' << phi[j][k] << '\n';
        }
        dataFile << '\n';
      }
      dataFile.close();
    }

    //imprimit la densidad
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
        E_i  =  E_i + pow(sqrt(vel_i[y][0] * vel_i[y][0] + vel_i[y][1] * vel_i[y][1]), 2) / (M_PI);

      for (int y = 0; y < le; y++)
        E_e  =  E_e + pow(sqrt(vel_e[y][0] * vel_e[y][0] + vel_e[y][1] * vel_e[y][1]), 2) / (M_PI * RAZON_MASAS);

      //Energía del campo
      for (int r = 0; r < J_X; r++)
        for(int ry = 0; ry < J_Y; ry++)
          E_field  =  E_field +  (ni[r][ry] - ne[r][ry]) * phi[r][ry] * hx * hx * hx;

      E_field  =  2 * E_field / (M_PI);

      E_total  =  E_field + E_i + E_e;

      //Energia perdida por partícula pérdida
      E_perdida  =   mv2perdidas / M_PI;

      fprintf(outEnergia,"%e %e %e %e %e %d  %e \n", kt * DT, E_total, E_i, E_e, E_field, total_e_perdidos+total_i_perdidos, E_perdida );
    }//Cierre de calculo de energia

    clock_t tiempo1  =  clock();
    if(kt % 1000 == 0) {
      cout << " CPU time 1000  =  " << double(tiempo1 - tiempo0) / CLOCKS_PER_SEC<< " sec" << endl;
      tiempo0  =  clock();
    }

    //Salida de función de distribución

    if(kt % 50000  ==  0) {
      sprintf(buffer,"fdist_ele%dx.data", kt);
      sprintf(bufferb,"fdist_ele%dy.data", kt);
      sprintf(bufferc,"fdist_ion%dx.data", kt);
      sprintf(bufferd,"fdist_ion%dy.data", kt);
      Funcion_Distribucion(pos_e,vel_e,le, buffer,  bufferb);
      Funcion_Distribucion(pos_i,vel_i,li, bufferc, bufferd);
    }

  } //Cierre del ciclo principal

  fclose(outEnergia);
  for(int i = 0;i<= 80;i++) {
    fclose(outFase_ele[i]);
    fclose(outFase_ion[i]);
  }

  return (0);
}// FINAL MAIN

//**********************************
//Función de Inyección de partículas
//**********************************

void initialize_Particles (double *pos_e, double *vel_e, double *pos_i, double *vel_i,
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


//**********************************
//Función de Inyección de partículas
//**********************************

void initialize_Particles (double pos_e[MAX_SPE][2], double vel_e[MAX_SPE][2], double pos_i[MAX_SPI][2], double vel_i[MAX_SPI][2],
    int li, int le) {
  for (int i = 0;i<MAX_SPE;i++) {
    pos_e[i+le][X] = 0;
    vel_e[i+le][X] =  create_Velocities_X (FE_MAXWELL_X, VPHI_E_X);
    pos_e[i+le][Y] = L_MAX_Y / 2.0;
    vel_e[i+le][Y] = create_Velocities_Y(FE_MAXWELL_Y, VPHI_E_Y);
    pos_i[i+li][X] = 0;
    vel_i[i+li][X] = create_Velocities_X (FI_MAXWELL_X, VPHI_I_X);
    pos_i[i+li][Y] = L_MAX_Y / 2.0;
    vel_i[i+li][Y] = create_Velocities_Y (FI_MAXWELL_Y, VPHI_I_Y);
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
/*
  if (flag  ==  0) {
    int seed  =  time (NULL);
    srand (seed);
    flag  =  1;
  }
*/
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
/*
  if (flag  ==  0) {
    int seed  =  time (NULL);
    srand (seed);
    flag  =  1;
  }
*/
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

  for (int i = 0;i<NSP;i++) {
    jr_x = pos[i]/hx; // indice (real) de la posición de la superpartícula
    j_x  = int(jr_x);    // indice  inferior (entero) de la celda que contiene a la superpartícula
    temp_x  =  jr_x-j_x;
    jr_y = pos[i * MAX_SPE] / hx; // indice (real) de la posición de la superpartícula
    j_y  = int(jr_y);    // indice  inferior (entero) de la celda que contiene a la superpartícula
    temp_y  =  jr_y-j_y;

    n[j_x + (j_y * J_X)] +=  (1.-temp_x)*(1.-temp_y)/(hx*hx*hx);
    n[(j_x+1) + (j_y * J_X)] +=  temp_x*(1.-temp_y)/(hx*hx*hx);
    n[j_x + ((j_y + 1) * J_X)] +=  (1.-temp_x)*temp_y/(hx*hx*hx);
    n[(j_x + 1) + ((j_y + 1) * J_X)] +=  temp_x*temp_y/(hx*hx*hx);
  }
}

void Concentration (double pos[MAX_SPE][2], double n[J_X][J_Y],int NSP, double hx) {
  int j_x,j_y;
  double temp_x,temp_y;
  double jr_x,jr_y;
  for (j_x = 0; j_x<J_X; j_x++) {
    for (j_y = 0; j_y<J_Y; j_y++) {
      n[j_x][j_y]  =  0.;
    }
  } // Inicializar densidad de carga

  for (int i = 0;i<NSP;i++) {
    jr_x = pos[i][0]/hx; // indice (real) de la posición de la superpartícula
    j_x  = int(jr_x);    // indice  inferior (entero) de la celda que contiene a la superpartícula
    temp_x  =  jr_x-j_x;
    jr_y = pos[i][1]/hx; // indice (real) de la posición de la superpartícula
    j_y  = int(jr_y);    // indice  inferior (entero) de la celda que contiene a la superpartícula
    temp_y  =  jr_y-j_y;

    n[j_x][j_y] +=  (1.-temp_x)*(1.-temp_y)/(hx*hx*hx);
    n[j_x+1][j_y] +=  temp_x*(1.-temp_y)/(hx*hx*hx);
    n[j_x][j_y+1] +=  (1.-temp_x)*temp_y/(hx*hx*hx);
    n[j_x+1][j_y+1] +=  temp_x*temp_y/(hx*hx*hx);

  }
}



//***********************************************************************
//Cálculo del Potencial Electrostático en cada uno de los puntos de malla
//***********************************************************************
void poisson2D_dirichletX_periodicY(double phi[J_X][J_Y],complex <double> rho[J_X][J_Y], double hx) {
  int M = J_X-2,N = J_Y;
  double h  =  hx;
  double hy  =  hx;
  double *f;
  fftw_complex  *f2;
  fftw_plan p,p_y,p_i,p_yi;
  f =  (double*) fftw_malloc(sizeof(double)* M);
  f2  =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  p  =  fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
  p_y  =  fftw_plan_dft_1d(N, f2, f2, FFTW_FORWARD, FFTW_ESTIMATE);
  p_i  =  fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
  p_yi  =  fftw_plan_dft_1d(N, f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);


  // Columnas FFT
  for (int k  =  0; k < N; k++) {
    for (int j  =  0; j < M; j++)
      f[j] = rho[j+1][k].real();
    fftw_execute(p);
    for (int j  =  0; j < M; j++)
      rho[j+1][k].real()  =  f[j];
  }

  // Filas FFT
  for (int j  =  0; j < M; j++) {
    for (int k  =  0; k < N; k++)
      memcpy( &f2[k], &rho[j+1][k], sizeof( fftw_complex ) );
    fftw_execute(p_y);
    for (int k  =  0; k < N; k++)
      memcpy( &rho[j+1][k], &f2[k], sizeof( fftw_complex ) );
  }

  // Resolver en el espacio de Fourier
  complex<double> i(0.0, 1.0);
  double pi  =  M_PI;
  complex<double> Wy  =  exp(2.0 * pi * i / double(N));
  complex<double> Wn  =  1.0;
  for (int m  =  0; m < M; m++) {
    for (int n  =  0; n < N; n++) {
      complex<double> denom  =  h*h*2.0+hy*hy*2.0;
      denom -=  hy*hy*(2*cos((m+1)*pi/(M+1))) + h*h*(Wn + 1.0 / Wn);
      if (denom !=  0.0)
        rho[m+1][n] *=  h*h*hy*hy / denom;
      Wn *=  Wy;
    }
  }

  // Inversa de las filas
  for (int j  =  0; j < M; j++) {
    for (int k  =  0; k < N; k++)
      memcpy( &f2[k], &rho[j+1][k], sizeof( fftw_complex ) );
    fftw_execute(p_yi);
    for (int k  =  0; k < N; k++) {
      memcpy( &rho[j+1][k], &f2[k], sizeof( fftw_complex ) );
      rho[j+1][k] /=  double(N); //La transformada debe ser normalizada.
    }
  }

  //Inversa Columnas FFT
  for (int k  =  0; k < N; k++) {
    for (int j  =  0; j < M; j++)
      f[j] = rho[j+1][k].real();
    fftw_execute(p_i);
    for (int j  =  0; j < M; j++)
      phi[j+1][k] = f[j]/double(2*(M+1));
  }

  for (int k  =  0; k < N; k++) {
    phi[0][k] = 0;
    phi[J_X-1][k] = 0;
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
      E_X[j + (k * J_X)] = (phi[(j-1) + (k * J_X)] - phi[(j + 1) + (k * J_X)]) / (2. * hx);
      E_Y[j + (k * J_X)] = (phi[j + ((k - 1) * J_X)] - phi[j + ((k + 1) * J_X)]) / (2. * hx);

      E_X[k * J_X] = 0;  //Cero en las fronteras X
      E_Y[k * J_X] = 0;
      E_X[(J_X-1) + (k * J_X)] = 0;
      E_Y[(J_X-1) + (k * J_X)] = 0;
    }

    E_X[j] = (phi[j - 1] - phi[j + 1]) / (2. * hx);
    E_Y[j] = (phi[j + ((J_Y - 1) * J_X)] - phi[j + J_X]) / (2. * hx);

    E_X[j + ((J_Y - 1) * J_X)] = (phi[(j - 1) + ((J_Y - 1) * J_X)] - phi[(j + 1) + ((J_Y - 1) * J_X)]) / (2. * hx);
    E_Y[j + ((J_Y - 1) * J_X)] = (phi[j + ((J_Y - 2) * J_X)] - phi[j]) / (2. * hx);

  }
}

void electric_field(double phi[J_X][J_Y], double E_X[J_X][J_Y], double E_Y[J_X][J_Y], double hx) {

  for (int j = 1;j<J_X-1;j++) {

    for (int k = 1;k<J_Y-1;k++) {
      E_X[j][k] = (phi[j-1][k]-phi[j+1][k])/(2.*hx);
      E_Y[j][k] = (phi[j][k-1]-phi[j][k+1])/(2.*hx);

      E_X[0][k] = 0;  //Cero en las fronteras X
      E_Y[0][k] = 0;
      E_X[J_X-1][k] = 0;
      E_Y[J_X-1][k] = 0;
    }

    E_X[j][0] = (phi[j-1][0]-phi[j+1][0])/(2.*hx);
    E_Y[j][0] = (phi[j][J_Y-1]-phi[j][1])/(2.*hx);

    E_X[j][J_Y-1] = (phi[j-1][J_Y-1]-phi[j+1][J_Y-1])/(2.*hx);
    E_Y[j][J_Y-1] = (phi[j][J_Y-2]-phi[j][0])/(2.*hx);

  }
}

//*******************************************************
// total_e_perdidos, total_i_perdidos, mv2perdidas,
void  Motion(double pos[MAX_SPE][2],  double vel[MAX_SPE][2],  int &NSP,
    int especie, double E_X[J_X][J_Y], double E_Y[J_X][J_Y], int kt, double hx,
    int &total_perdidos, double &mv2perdidas) {
  int j_x,j_y;
  double temp_x,temp_y,Ep_X, Ep_Y,fact;
  double jr_x,jr_y;
  int kk1 = 0;
  int conteo_perdidas = 0;

  if(especie ==  ELECTRONS)
    fact = FACT_EL;
  else
    fact = FACT_I;

  for (int i = 0;i<NSP;i++) {
    jr_x = pos[i][X]/hx;     // Índice (real) de la posición de la superpartícula (X)
    j_x  = int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
    temp_x  =  jr_x-double(j_x);
    jr_y = pos[i][Y]/hx;     // Índice (real) de la posición de la superpartícula (Y)
    j_y  = int(jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
    temp_y  =  jr_y-double(j_y);

    Ep_X = (1-temp_x)*(1-temp_y)*E_X[j_x][j_y]+
      temp_x*(1-temp_y)*E_X[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_X[j_x][j_y+1]+
      temp_x*temp_y*E_X[j_x+1][j_y+1];

    Ep_Y = (1-temp_x)*(1-temp_y)*E_Y[j_x][j_y]+
      temp_x*(1-temp_y)*E_Y[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_Y[j_x][j_y+1]+
      temp_x*temp_y*E_Y[j_x+1][j_y+1];

    vel[i][X] = vel[i][X] + CTE_E * FACTOR_CARGA_E * fact * Ep_X * DT;
    vel[i][Y] = vel[i][Y] + CTE_E * FACTOR_CARGA_E * fact * Ep_Y * DT;

    pos[i][X] = pos[i][X] + vel[i][X] * DT;
    pos[i][Y] = pos[i][Y] + vel[i][Y] * DT;

    if(pos[i][X]<0) {//Rebote en la pared del material.
      pos[i][X] = -pos[i][X];
      vel[i][X] = -vel[i][X];
    }

    if (pos[i][X]>= L_MAX_X) {//Partícula fuera del espacio de Simulación
      conteo_perdidas++;
      total_perdidos++;
      if(especie  ==  ELECTRONS) {
        printf("Electron perdido No. %d,  i = %d, kt = %d \n",total_perdidos, i ,kt);
        mv2perdidas+= pow( sqrt(vel[i][X]*vel[i][X]+vel[i][Y]*vel[i][Y]) , 2);
      }
      else {
        printf("Ion perdido No. %d,  i = %d, kt = %d \n",total_perdidos, i ,kt);
        mv2perdidas+= pow( sqrt(vel[i][X]*vel[i][X]+vel[i][Y]*vel[i][Y]) , 2)/(RAZON_MASAS);
      }
    }

    while(pos[i][Y] > L_MAX_Y) //Ciclo en el eje Y.
      pos[i][Y] = pos[i][Y]-L_MAX_Y;

    while(pos[i][Y] < 0.0) //Ciclo en el eje Y.
      pos[i][Y] = L_MAX_Y + pos[i][Y];

    if(pos[i][X]>= 0 && pos[i][X]<= L_MAX_X) {
      pos[kk1][X] = pos[i][X];
      pos[kk1][Y] = pos[i][Y];
      vel[kk1][X] = vel[i][X];
      vel[kk1][Y] = vel[i][Y];
      kk1++;
    }

    //Salida espacio de Fase
    if(kt % 10000  ==  0 && especie  ==  ELECTRONS);
    fprintf(outFase_ele[kt/10000]," %e   %e  %e  %e  %e \n",
        kt * DT, pos[i][X],vel[i][X],pos[i][Y],vel[i][Y]);

    if(kt % 10000  ==  0 && especie  ==  IONS)
      fprintf(outFase_ion[kt/10000]," %e   %e  %e  %e  %e \n",
          kt * DT, pos[i][X],vel[i][X],pos[i][Y],vel[i][Y]);

  }

  NSP -= conteo_perdidas;
  /*if(especie  ==  ELECTRONS)
    le = le-conteo_perdidas;
  else
    li = li-conteo_perdidas;
*/
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


void Funcion_Distribucion(double pos[MAX_SPE][2], double vel[MAX_SPE][2] , int NSP, char *archivo_X, char *archivo_Y) {
  double Nc = 100;
  FILE *pFile[2];
  pFile[0]  =  fopen(archivo_X,"w"); pFile[1]  =  fopen(archivo_Y,"w");
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
      fprintf(pFile[i]," %e  %d  \n", a, suma);
      a  =  a + dv;
    }
    fclose(pFile[i]);
  }
}
/*
//hx FACTOR_CARGA_I dt total_e_perdidos mv2perdidas
void  Motion_e(double pos[MAX_SPE][2], double vel[MAX_SPE][2], int NSP,
    double E_X[J_X][J_Y], double E_Y[J_X][J_Y], int kt){//,
  //int total_e_perdidos, double mv2perdidas) {
  int j_x,j_y;
  double temp_x,temp_y,Ep_X, Ep_Y,fact;
  double jr_x,jr_y;
  int kk1 = 0;
  int conteo_perdidas = 0;

  for (int i = 0;i<NSP;i++) {
    jr_x = pos[i][X]/hx;     // Índice (real) de la posición de la superpartícula (X)
    j_x  = int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
    temp_x  =  jr_x-double(j_x);
    jr_y = pos[i][Y]/hx;     // Índice (real) de la posición de la superpartícula (Y)
    j_y  = int(jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
    temp_y  =  jr_y-double(j_y);

    Ep_X = (1-temp_x)*(1-temp_y)*E_X[j_x][j_y]+
      temp_x*(1-temp_y)*E_X[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_X[j_x][j_y+1]+
      temp_x*temp_y*E_X[j_x+1][j_y+1];

    Ep_Y = (1-temp_x)*(1-temp_y)*E_Y[j_x][j_y]+
      temp_x*(1-temp_y)*E_Y[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_Y[j_x][j_y+1]+
      temp_x*temp_y*E_Y[j_x+1][j_y+1];

    vel[i][X] = vel[i][X]+CTE_E*FACTOR_CARGA_E*FACT_EL*Ep_X*dt;
    vel[i][Y] = vel[i][Y]+CTE_E*FACTOR_CARGA_E*FACT_EL*Ep_Y*dt;

    pos[i][X] = pos[i][X]+vel[i][X]*dt;
    pos[i][Y] = pos[i][Y]+vel[i][Y]*dt;

    if(pos[i][X]<0) {//Rebote en la pared del material.
      pos[i][X] = -pos[i][X];
      vel[i][X] = -vel[i][X];
    }

    if (pos[i][X]>= L_MAX_X) {//Partícula fuera del espacio de Simulación
      conteo_perdidas++;
      total_e_perdidos++;
      printf("Electron perdido No. %d,  i = %d, kt = %d \n",total_e_perdidos, i ,kt);
      mv2perdidas+= pow( sqrt(vel[i][X]*vel[i][X]+vel[i][Y]*vel[i][Y]) , 2);
    }

    while(pos[i][Y]>L_MAX_Y) //Ciclo en el eje Y.
      pos[i][Y] = pos[i][Y]-L_MAX_Y;

    while(pos[i][Y]<0.0) //Ciclo en el eje Y.
      pos[i][Y] = L_MAX_Y+pos[i][Y];

    if(pos[i][X]>= 0.0 && pos[i][X]<= L_MAX_X) {
      pos[kk1][X] = pos[i][X];
      pos[kk1][Y] = pos[i][Y];
      vel[kk1][X] = vel[i][X];
      vel[kk1][Y] = vel[i][Y];
      kk1++;
    }

    //Salida espacio de Fase
    if(kt%10000 == 0)
      fprintf(outFase_ele[kt/10000]," %e   %e  %e  %e  %e \n",kt*dt,pos[i][X],vel[i][X],pos[i][Y],vel[i][Y]);
  }

  le = le-conteo_perdidas;
}

// FACTOR_CARGA_I FACT_I dt  mv2perdidas
void  Motion_i(double pos[MAX_SPE][2], double vel[MAX_SPE][2], int NSP,
    double E_X[J_X][J_Y], double E_Y[J_X][J_Y],
    int total_i_perdidos, double mv2perdidas, int kt) {

  int j_x,j_y;
  double temp_x,temp_y,Ep_X, Ep_Y,fact;
  double jr_x,jr_y;
  int kk1 = 0;
  int conteo_perdidas = 0;


  for (int i = 0;i<NSP;i++) {
    jr_x = pos[i][X]/hx;     // Índice (real) de la posición de la superpartícula (X)
    j_x  = int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
    temp_x  =  jr_x-double(j_x);
    jr_y = pos[i][Y]/hx;     // Índice (real) de la posición de la superpartícula (Y)
    j_y  = int(jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
    temp_y  =  jr_y-double(j_y);

    Ep_X = (1-temp_x)*(1-temp_y)*E_X[j_x][j_y]+
      temp_x*(1-temp_y)*E_X[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_X[j_x][j_y+1]+
      temp_x*temp_y*E_X[j_x+1][j_y+1];

    Ep_Y = (1-temp_x)*(1-temp_y)*E_Y[j_x][j_y]+
      temp_x*(1-temp_y)*E_Y[j_x+1][j_y]+
      (1-temp_x)*temp_y*E_Y[j_x][j_y+1]+
      temp_x*temp_y*E_Y[j_x+1][j_y+1];

    vel[i][X] = vel[i][X]+CTE_E*FACTOR_CARGA_E*FACT_I*Ep_X*dt;
    vel[i][Y] = vel[i][Y]+CTE_E*FACTOR_CARGA_E*FACT_I*Ep_Y*dt;

    pos[i][X] = pos[i][X]+vel[i][X]*dt;
    pos[i][Y] = pos[i][Y]+vel[i][Y]*dt;

    if(pos[i][X]<0) {//Rebote en la pared del material.
      pos[i][X] = -pos[i][X];
      vel[i][X] = -vel[i][X];
    }

    if (pos[i][X]>= L_MAX_X) {//Partícula fuera del espacio de Simulación
      conteo_perdidas++;
      total_i_perdidos++;
      printf("Ion perdido No. %d,  i = %d, kt = %d \n",total_i_perdidos, i ,kt);
      mv2perdidas+= pow( sqrt(vel[i][X]*vel[i][X]+vel[i][Y]*vel[i][Y]) , 2)/(RAZON_MASAS);
    }

    while(pos[i][Y]>L_MAX_Y) //Ciclo en el eje Y.
      pos[i][Y] = pos[i][Y]-L_MAX_Y;

    while(pos[i][Y]<0.0) //Ciclo en el eje Y.
      pos[i][Y] = L_MAX_Y+pos[i][Y];

    if(pos[i][X]>= 0.0 && pos[i][X]<= L_MAX_X) {
      pos[kk1][X] = pos[i][X];
      pos[kk1][Y] = pos[i][Y];
      vel[kk1][X] = vel[i][X];
      vel[kk1][Y] = vel[i][Y];
      kk1++;
    }

    //Salida espacio de Fase
    if(kt%10000 == 0)
      fprintf(outFase_ion[kt/10000]," %e   %e  %e  %e  %e \n",kt*dt,pos[i][X],vel[i][X],pos[i][Y],vel[i][Y]);
  }

  li = li-conteo_perdidas;
}
*/

