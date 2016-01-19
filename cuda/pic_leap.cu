#include "pic_cuda.cu"

using namespace std;
using namespace pic_cuda;

int main() {
  clock_t tiempo0  =  clock();

  int le = 0, li = 0;
  double  t_0;
  int  total_e_perdidos = 0;
  int  total_i_perdidos = 0;
  double  mv2perdidas = 0;

  double  ND = NE03D * pow(LAMBDA_D,3);                          //Parámetro del plasma
  int     NTe = 1e5;                                  //Número de partículas "reales"
  int     k_MAX_inj;   //Tiempo máximo de inyección
  int     K_total;     //Tiempo total
  int     Ntv = 8;
  int     NTSPe, MAX_SPE_dt, MAX_SPI_dt;
  double  cte_rho = pow(E_CHARGE * T0, 2) / (M_I * EPSILON_0 * pow(X0, 3)); //Normalización de EPSILON_0
 // FILE    *outEnergia;


  //***************************
  //Constantes de normalización
  //***************************

  //double  X0 = LAMBDA_D;                //Escala de longitud: Longitud de Debye
  double  n0  =  double(NTe) / (X0 * X0 * X0);
  double  ni0_3D  =  NI03D * pow(X0, 3);
  double  ne0_3D  =  NE03D * pow(X0, 3);
  double hx;
  int seed  =  time (NULL);
  srand (seed);  // Semilla para generar números aleatorios dependiendo del reloj interno.

  //******************
  //ARCHIVOS DE SALIDA
  //******************
  //outEnergia = fopen("Energia","w");
  char buffer[40];//, bufferb[40], bufferc[40], bufferd[40];
  /*
     for(int i  =  0; i<= 80; i++) {
     sprintf(buffer, "fase_ele%d", i);
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
 //    outFase_ele[i] = fopen(buffer, "w");

  //****************************************
  // Inicialización de variables del sistema
  //****************************************
  int size = MAX_SPE * 2 * sizeof(double);
  int size1 = J_X * J_Y * sizeof(double);
  int size2 = J_X * J_Y * sizeof(complex<double>);

  double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y, *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y,
         *ne, *ni, *phi, *E_X, *E_Y;
  complex<double> *rho;
  //double  E_i, E_e, E_field;


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
  rho   = (complex<double> *) malloc(size2);



  //***************************
  // Normalización de variables
  //***************************

  t_0 = 1;
  hx = DELTA_X / X0;                            // Paso espacial
  double n_0 = n0 * X0 * X0 * X0;                   // Densidad de partículas
  NTSPe = NTe / FACTOR_CARGA_E;
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

  //Velocidades y posiciones iniciales de las partículas (no las libera).

  initialize_Particles (pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, FE_MAXWELL_X, FE_MAXWELL_Y, VPHI_E_X, VPHI_E_Y);
  initialize_Particles (pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, FI_MAXWELL_X, FI_MAXWELL_Y, VPHI_I_X, VPHI_I_Y);

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


    Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
/*
    H_Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    gpu_error(cudaGetLastError());
    H_Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    gpu_error(cudaGetLastError());
*/
    for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    poisson2D_dirichletX_periodicY(phi, rho, hx);
    // Calcular campo eléctrico en puntos de malla


    electric_field(phi, E_X, E_Y, hx);
  /*
    H_electric_field(phi, E_X, E_Y, hx);
    gpu_error(cudaGetLastError());
*/
    // imprimir el potencial electroestatico.
    if(kt % 100)
      cout << "100" << endl;
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
    Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);

 /*   H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
*/
    //Cálculo de energías.
 /*   if(kt % 2000 == 0 && kt > 0) {
      E_i = 0; //Se inicializan las variables de la energía
      E_e = 0;
      E_field = 0;

      for (int y = 0; y < li; y++)
        E_i  =  E_i + pow(sqrt(vel_i[y] * vel_i[y] + vel_i[y + MAX_SPI] * vel_i[y + MAX_SPI]), 2) / (M_PI);

      for (int y = 0; y < le; y++)
        E_e  =  E_e + pow(sqrt(vel_e[y] * vel_e[y] + vel_e[y + MAX_SPE] * vel_e[y + MAX_SPE]), 2) / (M_PI * RAZON_MASAS);

      //Energía del campo
      for (int r = 0; r < J_X; r++)
        for(int ry = 0; ry < J_Y; ry++)
          E_field  =  E_field +  (ni[(r * J_Y) + ry] - ne[(r * J_Y) + ry]) * phi[(r * J_Y) + ry] * hx * hx * hx;

      E_field  =  2 * E_field / (M_PI);


      //Energia perdida por partícula pérdida

 //     fprintf(outEnergia,"%e %e %e %e %e %d  %e \n", kt * DT, E_total, E_i, E_e, E_field, total_e_perdidos + total_i_perdidos, E_perdida );
    }//Cierre de calculo de energia*/

    clock_t tiempo1  =  clock();
    /*if(kt % 50000 == 0) {
      cout << " CPU time " << kt / 50000 << "  =  " << double(tiempo1 - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;
      tiempo0  =  clock();
    }*/

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

  //fclose(outEnergia);
  /*for(int i = 0;i<= 80;i++) {
    fclose(outFase_ele[i]);
    fclose(outFase_ion[i]);
  }*/

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
  cout << " CPU time = " << double(clock() - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;

  return (0);
}// FINAL MAIN

