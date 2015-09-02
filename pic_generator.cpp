#include "pic.cpp"

using namespace std;
using namespace pic;

void write_test_concentration(double *pos, double *n, int NSP, double hx);
void write_test_poisson(double *phi, complex<double> *rho, double hx);
void write_test_out_poisson(double *phi, complex<double> *rho);
void write_test_electric(double *phi, double *E_X, double *E_Y, double hx);
void write_test_Motion(double *pos, double *vel, int &NSP, int especie, double *E_X,
    double *E_Y, int kt, double hx, int &total_perdidos, double &mv2perdidas);
void write_test_out_Motion(double *pos, double *vel, int &NSP, double *E_X, double *E_Y,
    int &total_perdidos, double &mv2perdidas);

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
  // FILE    *outEnergia;


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
  //outEnergia = fopen("Energia","w");
  char buffer[40], bufferb[40], bufferc[40], bufferd[40];
  /*
     for(int i  =  0; i<= 80; i++) {
     sprintf(buffer, "fase_ele%d", i);
     }

     for(int i  =  0; i<= 80; i++) {
     sprintf(buffer, "fase_ion%d", i);
     outFase_ion[i] = fopen(buffer, "w");
     }
     */

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
    if(kt <= k_MAX_inj && kt == kk) {// Inyectar superpartículas (i-e)
      le+= MAX_SPE_dt;
      li+= MAX_SPI_dt;
      kk = kk + Kemision;
    }
    //-----------------------------------------------
    // Calculo de "densidad de carga 2D del plasma"

    Concentration (pos_e, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    Concentration (pos_i, ni, li, hx);// Calcular concentración de superpartículas Iónicas

    if(kt % 25000 == 0) {
      write_test_concentration(pos_e, ne, le, hx);
      write_test_concentration(pos_i, ni, li, hx);
    }

    for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    if(kt % 25000 == 0) {
      write_test_poisson(phi, rho, hx);
    }
    poisson2D_dirichletX_periodicY(phi, rho, hx);
    if(kt % 25000 == 0) {
      write_test_out_poisson(phi, rho);
    }
    // Calcular campo eléctrico en puntos de malla
    electric_field(phi, E_X, E_Y, hx);
    if(kt % 25000 == 0) {
      write_test_electric(phi, E_X, E_Y, hx);
    }

    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    if(kt % 25000 == 0) {
      write_test_Motion(pos_e, vel_e, le, ELECTRONS, E_X, E_Y, kt, hx, total_e_perdidos, mv2perdidas);
      write_test_Motion(pos_i, vel_i, li, IONS, E_X, E_Y, kt, hx, total_i_perdidos, mv2perdidas);
    }
    Motion(pos_e, vel_e, le, ELECTRONS, E_X, E_Y, kt, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_i, vel_i, li, IONS, E_X, E_Y, kt, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    if(kt % 25000 == 0) {
      write_test_out_Motion(pos_e, vel_e, le, E_X, E_Y, total_e_perdidos, mv2perdidas);
      write_test_out_Motion(pos_i, vel_i, li, E_X, E_Y, total_i_perdidos, mv2perdidas);
    }

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

    }//Cierre de calculo de energia

  } //Cierre del ciclo principal
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

void write_test_concentration(double *pos, double *n, int NSP, double hx) {
  static int num = 0;
  char name[40] = "test_concentration";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << MAX_SPE * 2 << endl;
  for(int i = 0; i < MAX_SPE * 2; i++)
    cout << pos[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << n[i] << " ";
  cout << endl << NSP << " " << hx << endl;
  fclose(stdout);
}

void write_test_poisson(double *phi, complex<double> *rho, double hx) {
  static int num = 0;
  char name[40] = "test_poisson";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << phi[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rho[i] << " ";
  cout << endl << hx << endl;
  fclose(stdout);
}

void write_test_out_poisson(double *phi, complex<double> *rho) {
  static int num = 0;
  char name[40] = "test_out_poisson";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << phi[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rho[i] << " ";
  fclose(stdout);
}

void write_test_electric(double *phi, double *E_X, double *E_Y, double hx) {
  static int num = 0;
  char name[40] = "test_electric";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << phi[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_X[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_Y[i] << " ";
  cout << endl << hx << endl;
  fclose(stdout);
}

void write_test_Motion(double *pos, double *vel, int &NSP, int especie, double *E_X,
    double *E_Y, int kt, double hx, int &total_perdidos, double &mv2perdidas) {
  static int num = 0;
  char name[40] = "test_Motion";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << pos[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << vel[i] << " ";
  cout << endl << NSP << endl << especie << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_X[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_Y[i] << " ";
  cout << endl << kt << endl << hx << endl << total_perdidos << endl << mv2perdidas << endl;
  fclose(stdout);
}

void write_test_out_Motion(double *pos, double *vel, int &NSP, double *E_X, double *E_Y,
    int &total_perdidos, double &mv2perdidas) {
  static int num = 0;
  char name[40] = "test_out_Motion";
  sprintf(name, "%d.data\0", num++);
  freopen (name, "w", stdout);
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << pos[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << vel[i] << " ";
  cout << endl << NSP << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_X[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_Y[i] << " ";
  cout << endl << total_perdidos << endl << mv2perdidas << endl;
  fclose(stdout);
}
