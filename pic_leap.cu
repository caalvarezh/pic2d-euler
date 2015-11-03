#include "pic.cu"
#include "pic.cpp"

using namespace std;
using namespace pic_cuda;

void write_test_concentration(double *pos_x, double *pos_y, double *n, int NSP, double hx) {
  static int num = 0;
  char name[40];
  sprintf(name, "test_concentration%d.data", num++);
  freopen (name, "w", stdout);
  cout << NSP << endl;
  for(int i = 0; i < NSP; i++)
    cout << pos_x[i] << " ";
  for(int i = 0; i < NSP; i++)
    cout << pos_y[i] << " ";

  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << n[i] << " ";
  cout << endl << NSP << " " << hx << endl;
  fclose(stdout);
}

void write_test_electric(double *phi, double *E_X, double *E_Y, double hx) {
  static int num = 0;
  char name[40];
  sprintf(name, "test_electric%d.data", num++);
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

void write_test_Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int &NSP, int especie,
     double *E_X, double *E_Y, double hx, int &total_perdidos, double &mv2perdidas) {
  static int num = 0;
  char name[40];
  sprintf(name, "test_Motion%d.data", num++);
  freopen (name, "w", stdout);

  cout << NSP << endl;
  for(int i = 0; i < NSP; i++)
    cout << pos_x[i] << " ";
  for(int i = 0; i < NSP; i++)
    cout << pos_y[i] << " ";

  for(int i = 0; i < NSP; i++)
    cout << vel_x[i] << " ";
  for(int i = 0; i < NSP; i++)
    cout << vel_y[i] << " ";

  cout << endl << especie << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_X[i] << " ";
  cout << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << E_Y[i] << " ";
  cout << endl << hx << endl << total_perdidos << endl << mv2perdidas << endl;
  fclose(stdout);
}

void write_test_out_Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int NSP,
    int NSP_out) {
  static int num = 0;
  char name[40] = "";
  sprintf(name, "test_out_Motion%d.data", num++);
  freopen (name, "w", stdout);
  cout << NSP << endl;
  for(int i = 0; i < NSP; i++)
    cout << pos_x[i] << " ";
  for(int i = 0; i < NSP; i++)
    cout << pos_y[i] << " ";

  for(int i = 0; i < NSP; i++)
    cout << vel_x[i] << " ";
  for(int i = 0; i < NSP; i++)
    cout << vel_y[i] << " ";

  cout << endl << NSP_out << endl;
  fclose(stdout);
}

int main() {
  //************************
  // Parámetros del sistema
  //************************
  clock_t tiempo0  =  clock();

  int le = 0, li = 0;
  double  t_0, x_0;
  int  total_e_perdidos = 0;
  int  total_i_perdidos = 0;
  double  mv2perdidas = 0;

  double  ND = NE03D * pow(LAMBDA_D,3);                          //Parámetro del plasma
  int     k_MAX_inj;   //Tiempo máximo de inyección
  int     K_total;     //Tiempo total
  int     Ntv = 8;
  int     NTSPe, NTSPI, MAX_SPE_dt, MAX_SPI_dt;
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
  //int seed  =  time (NULL);
  //srand (seed);  // Semilla para generar números aleatorios dependiendo del reloj interno.
  //******************
  //ARCHIVOS DE SALIDA
  //******************
  //outEnergia = fopen("Energia","w");
  char buffer[40];//, bufferb[40], bufferc[40], bufferd[40];

  //****************************************
  // Inicialización de variables del sistema
  //****************************************
  int size = MAX_SPE * sizeof(double);
  int size1 = J_X * J_Y * sizeof(double);
  int size2 = J_X * J_Y * sizeof(complex<double>);

  double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y, *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y, *ne, *ni;
  double *phi, *E_X, *E_Y;
  complex<double> *rho;
  //double  E_i,E_e,E_field,E_total,E_perdida;

  pos_e_x = (double *) malloc(size);
  pos_e_y = (double *) malloc(size);
  pos_i_x = (double *) malloc(size);
  pos_i_y = (double *) malloc(size);

  vel_e_x = (double *) malloc(size);
  vel_e_y = (double *) malloc(size);
  vel_i_x = (double *) malloc(size);
  vel_i_y = (double *) malloc(size);
  ne    = (double *) malloc(size1);
  ni    = (double *) malloc(size1);
  phi   = (double *) malloc(size1);
  E_X   = (double *) malloc(size1);
  E_Y   = (double *) malloc(size1);
  phi   = (double *) malloc(size1);
  rho   = (complex<double> *) malloc(size2);

  //***************************
  // Normalización de variables
  //***************************

  t_0 = 1;
  x_0 = 1;
  hx = DELTA_X / X0;                            // Paso espacial
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

  initialize_Particles (pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, FE_MAXWELL_X, FE_MAXWELL_Y, VPHI_E_X, VPHI_E_Y);//Velocidades y posiciones iniciales de las partículas (no las libera).
  initialize_Particles (pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, FI_MAXWELL_X, FI_MAXWELL_Y, VPHI_I_X, VPHI_I_Y);//Velocidades y posiciones iniciales de las partículas (no las libera).


  double tacum = 0;
  for(int kk  =  0, kt  =  0; kt <= K_total; kt++) {
    cout << hx << endl;
    break;
    //cerr << kt << endl;
    /*if(kt % 50000 == 0) {
      printf("kt = %d\n", kt);
      printf("le = %d   li = %d \n",le, li );
    }*/
    if(kt <= k_MAX_inj && kt == kk) {// Inyectar superpartículas (i-e)
      le+= MAX_SPE_dt;
      li+= MAX_SPI_dt;
      kk = kk + Kemision;
    }
    //-----------------------------------------------
    // Calculo de "densidad de carga 2D del plasma"
    //H_Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    //H_Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    if(kt % 2000 == 0) {
      write_test_concentration(pos_e_x, pos_e_y, ne, le, hx);
      write_test_concentration(pos_i_x, pos_i_y, ni, le, hx);
    }
    pic::Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    pic::Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas

    //H_rhoKernel(ne, ni, rho);
     for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    //poisson2D_dirichletX_periodicY(phi, rho, hx);
    pic::poisson2D_dirichletX_periodicY(phi, rho, hx);
    // Calcular campo eléctrico en puntos de malla
    if(kt % 2000 == 0)
      write_test_electric(phi, E_X, E_Y, hx);
    H_electric_field(phi, E_X, E_Y, hx);
    //pic::electric_field(phi, E_X, E_Y, hx);

    // imprimir el potencial electroestatico.
    if(kt % 1  !=  0) {
      //cout << "le: " << le << " li: " << li << endl;
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
    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    //H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    //H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    if(kt % 2000 == 0) {
      write_test_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);
      write_test_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);
    }
    pic::Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    pic::Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);

    if(kt % 1 != 0) {
      cout << " CPU time " << kt / 1000 << "  =  " << tacum / CLOCKS_PER_SEC << " sec" << endl;
    }

    //Salida de función de distribución

  } //Cierre del ciclo principal
  free(pos_e_x);
  free(pos_e_y);
  free(pos_i_x);
  free(pos_i_y);
  free(vel_e_x);
  free(vel_e_y);
  free(vel_i_x);
  free(vel_i_y);
  free(ne);
  free(ni);
  free(phi);
  free(E_X);
  free(E_Y);
  free(rho);
  //cout << " CPU time  = " << double(clock() - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;

  return (0);
}// FINAL MAIN

