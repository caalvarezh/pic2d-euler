#include "pic_cuda.cu"

using namespace std;
using namespace pic_cuda;

#define gpu_error(ans) { gpu_assert((ans), __LINE__); }

inline void gpu_assert(cudaError_t code, int line){
  if (code != cudaSuccess)
    cerr<<"GPUerror: "<<cudaGetErrorString(code)<<" in "<< line<<endl;
}

void initialize_vectors (double *pos_x, double *pos_y, double *vel_x, double *vel_y, double *E_X) {
  int NSP = MAX_SPE;
  for(int i = 0; i < NSP; i++) {
    pos_x[i] = rand() % int(L_MAX_X);
    pos_y[i] = rand() % int(L_MAX_Y);
    vel_x[i] = (rand() % 100) / 2.0;
    vel_y[i] = (rand() % 100) / 2.0;
  }
  for(int i = 0; i < J_X * J_Y; i++) {
    E_X[i] = rand() % 2342;
  }
}

int main() {
  //************************
  // Parámetros del sistema
  //************************

  int le = MAX_SPE, li = MAX_SPE;
  int  total_e_perdidos = 0;
  int  total_i_perdidos = 0;
  double  mv2perdidas = 0;

  double  ND = NE03D * pow(LAMBDA_D,3);                          //Parámetro del plasma
 // FILE    *outEnergia;


  //***************************
  //Constantes de normalización
  //***************************

  //double  X0 = LAMBDA_D;                //Escala de longitud: Longitud de Debye
  double  ni0_3D  =  NI03D * pow(X0, 3);
  double  ne0_3D  =  NE03D * pow(X0, 3);

  int size = MAX_SPE * sizeof(double);
  int size1 = J_X * J_Y * sizeof(double);

  double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y, *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y, *ne, *ni;
  double *phi, *E_X, *E_Y;
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

  //***************************
  // Normalización de variables
  //***************************

  double hx = DELTA_X / X0;                            // Paso espacial
  int max_it = 20;

  double tcon, tscon, telec, tselec, tmot, tsmot;
  tcon = tscon = telec = tselec = tmot = tsmot = 0.0 ;
  clock_t tiempo;
  cout << "start " << endl;
  for(int it  =  0; it <= max_it; it++) {
    cout << it << endl;
    initialize_vectors (pos_e_x, pos_e_y, vel_e_x, vel_e_y, E_X);
    initialize_vectors (pos_i_x, pos_i_y, vel_i_x, vel_i_y, E_Y);
    for(int i = 0; i < J_X * J_Y; i++)
      phi[i] =  rand() % 8234;

    // Calculo de "densidad de carga 2D del plasma"
    tiempo = clock();
    H_Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    gpu_error(cudaGetLastError());
    H_Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    gpu_error(cudaGetLastError());
    tcon += clock() - tiempo;

    tiempo = clock();
    Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    tscon += clock() - tiempo;

    // Calcular campo eléctrico en puntos de malla
    tiempo = clock();
    H_electric_field(phi, E_X, E_Y, hx);
    gpu_error(cudaGetLastError());
    telec += clock() - tiempo;

    tiempo = clock();
    electric_field(phi, E_X, E_Y, hx);
    tselec += clock() - tiempo;

    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    tiempo = clock();
    H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
    tmot += clock() - tiempo;

    tiempo = clock();
    Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tsmot += clock() - tiempo;

  } //Cierre del ciclo principal
  cout << "Concentration\nGPU = " << tcon / CLOCKS_PER_SEC << " sec  CPU = " << tscon / CLOCKS_PER_SEC << endl;
  cout << "Electric field\nGPU = " << telec / CLOCKS_PER_SEC << " sec  CPU = " << tsmot / CLOCKS_PER_SEC << endl;
  cout << "Motion\nGPU = " << tmot / CLOCKS_PER_SEC << " sec CPU = " << tsmot / CLOCKS_PER_SEC << endl;
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

  return (0);
}// FINAL MAIN

