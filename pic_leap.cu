#include "pic.cu"
#include "pic.cpp"

using namespace std;
using namespace pic_cuda;

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
  int seed  =  time (NULL);
  srand (seed);  // Semilla para generar números aleatorios dependiendo del reloj interno.

  //******************
  //ARCHIVOS DE SALIDA
  //******************
  //outEnergia = fopen("Energia","w");
  char buffer[40], bufferb[40], bufferc[40], bufferd[40];

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

  for(int kk  =  0, kt  =  0; kt <= K_total; kt++) {
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

    H_Concentration (pos_e, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    H_Concentration (pos_i, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    //pic::Concentration (pos_e, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    //pic::Concentration (pos_i, ni, li, hx);// Calcular concentración de superpartículas Iónicas

    //H_rhoKernel(ne, ni, rho);
     for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    //poisson2D_dirichletX_periodicY(phi, rho, hx);
    pic::poisson2D_dirichletX_periodicY(phi, rho, hx);
    // Calcular campo eléctrico en puntos de malla

    H_electric_field(phi, E_X, E_Y, hx);
    //pic::electric_field(phi, E_X, E_Y, hx);

    // imprimir el potencial electroestatico.
    if(kt % 50000  ==  0) {
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

    H_Motion(pos_e, vel_e, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    H_Motion(pos_i, vel_i, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    //pic::Motion(pos_e, vel_e, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    //pic::Motion(pos_i, vel_i, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);

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

 //     fprintf(outEnergia,"%e %e %e %e %e %d  %e \n", kt * DT, E_total, E_i, E_e, E_field, total_e_perdidos + total_i_perdidos, E_perdida );
    }//Cierre de calculo de energia

    /*clock_t tiempo1  =  clock();
    if(kt % 50000 == 0) {
      cout << " CPU time " << kt / 50000 << "  =  " << double(tiempo1 - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;
      tiempo0  =  clock();
    }*/

    //Salida de función de distribución

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
  cout << " CPU time  = " << double(clock() - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;

  return (0);
}// FINAL MAIN

