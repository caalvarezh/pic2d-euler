#ifndef PICLEAP
#define PICLEAP

#include <cmath>

namespace pic {
  #define __d {prueba (__LINE__);}

  const int MAX_SPE     = 10000;           // Limite (computacional) de Superpartículas electrónicas
  const int MAX_SPI     = 10000;           // Limite (computacional) de Superpartículas iónicas
  const int J_X         = 129;           // Número de puntos de malla X. Recomendado: Del orden 2^n+1
  const int J_Y         = 64;           // Número de puntos de malla Y. Recomendado: Del orden 2^n
  const int ELECTRONS   = 0;
  const int IONS        = 1;
  const int X           = 0;
  const int Y           = 1;
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
  const double DT = 1.e-5;                                  // Paso temporal
  const double VFLUX_I_magnitud =  std::sqrt(VFLUX_I * VFLUX_I + VFLUX_I * VFLUX_I); // Velocidad de flujo iónico (m/s)  =  sqrt(2*K_BOLTZMANN*Te/(M_PI*M_I))
  const double vflux_e_magnitud =  std::sqrt(VFLUX_E_X * VFLUX_E_X + VFLUX_E_Y * VFLUX_E_Y);
  const double Te = (M_PI * 0.5 * E_MASS * (std::pow(VFLUX_E_X, 2) / K_BOLTZMANN));    // Temperatura electrónica inicial (°K)

  const double NI03D = (FLUJO_INICIAL / VFLUX_I);
  const double NE03D = (FLUJO_INICIAL / VFLUX_E_X);
  const double LAMBDA_D = std::sqrt(EPSILON_0 * K_BOLTZMANN * Te / (NE03D * std::pow(E_CHARGE, 2)));  //Longitud de Debye
  const double DELTA_X = (LAMBDA_D);   //Paso espacial
  const double L_MAX_X = (((J_X-1) * DELTA_X) / X0);                      // Longitud región de simulación
  const double L_MAX_Y = (((J_Y-1) * DELTA_X) / X0);                      // Longitud región de simulación

  void prueba(int d);

  void   initialize_Particles(double *pos_e, double *vel_e, double *pos_i,
      double *vel_i, int li, int le);

  double create_Velocities_X(double fmax, double vphi);

  double create_Velocities_Y(double fmax, double vphi);

  void   Concentration(double *pos, double *n, int NSP, double hx);

  void   poisson2D_dirichletX_periodicY(double *phi, complex <double> *rho,
      double hx);

  void   electric_field(double *phi, double *E_X, double *E_Y, double hx);

  void Motion(double *pos, double *vel, int &NSP, int especie, double *E_X,
      double *E_Y, int kt, double hx, int &total_perdidos, double &mv2perdidas);

  void   Funcion_Distribucion(double *pos, double *vel, int NSP, char *archivo_X,
      char *archivo_Y);

  void   Funcion_Distribucion(double pos[MAX_SPE][2], double vel[MAX_SPE][2],
      int NSP, char *archivo_X, char *archivo_Y);
}
#endif
