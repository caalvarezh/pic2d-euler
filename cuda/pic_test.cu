#include "pic_cuda.cu"
#define ELE 1
#define CON 0
#define MOT 0
using namespace std;
using namespace pic_cuda;

#define gpu_error(ans) { gpu_assert((ans), __LINE__); }

inline void gpu_assert(cudaError_t code, int line){
  if (code != cudaSuccess)
    cerr<<"GPUerror: "<<cudaGetErrorString(code)<<" in "<< line<<endl;
}

bool compare(double *c1, double *c2, int size) {
  int cont = 0;
  for(int i = 0; i < size; i++) {
    if(fabs(c1[i] - c2[i]) > 1e-2) {
      cout << "Fail in " << i << endl;
      cout << c1[i] << " " << c2[i] << endl;
      cont++;
    }
  }
  if(cont > 0) {
    cout << "total " << cont << " fails" << endl;
    return true;
  }
  return false;
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

  double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y, *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y,
         *ne, *ne1, *ni, *ni1, *phi, *E_X, *E_Y, *E_X1, *E_Y1, *pos_ex, *pos_ey, *vel_ex, *vel_ey,
         *pos_ix, *pos_iy, *vel_ix, *vel_iy;
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

  pos_ex = (double *) malloc(size);
  pos_ey = (double *) malloc(size);
  pos_ix = (double *) malloc(size);
  pos_iy = (double *) malloc(size);

  vel_ex = (double *) malloc(size);
  vel_ey = (double *) malloc(size);
  vel_ix = (double *) malloc(size);
  vel_iy = (double *) malloc(size);
  ne1   = (double *) malloc(size1);
  ni1   = (double *) malloc(size1);
  E_X1  = (double *) malloc(size1);
  E_Y1  = (double *) malloc(size1);

  //***************************
  // Normalización de variables
  //***************************

  double hx = DELTA_X / X0;                            // Paso espacial
  int max_it = 1;

  double tcon, tscon, telec, tselec, tmot, tsmot;
  tcon = tscon = telec = tselec = tmot = tsmot = 0.0 ;
  clock_t tiempo;
  cout << "start " << endl;
  for(int it  =  0; it < max_it; it++) {
    initialize_vectors (pos_e_x, pos_e_y, vel_e_x, vel_e_y, E_X);
    initialize_vectors (pos_i_x, pos_i_y, vel_i_x, vel_i_y, E_Y);
    memcpy(E_X1, E_X, size1);
    memcpy(E_Y1, E_Y, size1);
    memcpy(pos_ex, pos_e_x, size);
    memcpy(pos_ix, pos_i_x, size);
    memcpy(vel_ex, vel_e_x, size);
    memcpy(vel_ix, vel_i_x, size);
    memcpy(pos_ey, pos_e_y, size);
    memcpy(pos_iy, pos_i_y, size);
    memcpy(vel_ey, vel_e_y, size);
    memcpy(vel_iy, vel_i_y, size);
    for(int i = 0; i < J_X * J_Y; i++)
      phi[i] =  rand() % 8234;

    // Calculo de "densidad de carga 2D del plasma"
#if CON
    tiempo = clock();
    H_Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    gpu_error(cudaGetLastError());
    H_Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    gpu_error(cudaGetLastError());
    tcon += clock() - tiempo;

    tiempo = clock();
    Concentration (pos_e_x, pos_e_y, ne1, le, hx);// Calcular concentración de superpartículas electrónicas
    Concentration (pos_i_x, pos_i_y, ni1, li, hx);// Calcular concentración de superpartículas Iónicas
    tscon += clock() - tiempo;

    if(compare(ne, ne1, J_X * J_Y))
      cout << "fail ne" << endl;
    if(compare(ni, ni1, J_X * J_Y))
      cout << "fail ni" << endl;

    // Calcular campo eléctrico en puntos de malla
#endif

#if ELE
    tiempo = clock();
    H_electric_field(phi, E_X, E_Y, hx);
    gpu_error(cudaGetLastError());
    telec += clock() - tiempo;

    tiempo = clock();
    electric_field(phi, E_X1, E_Y1, hx);
    tselec += clock() - tiempo;

    if(compare(E_X, E_X1, J_X * J_Y))
      cout << "----- fail ex -----" << endl << endl;

    if(compare(E_Y, E_Y1, J_X * J_Y))
      cout << "----- fail ey -----" << endl << endl;
#endif

#if MOT
    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    tiempo = clock();
    H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    gpu_error(cudaGetLastError());
    tmot += clock() - tiempo;

    tiempo = clock();
    Motion(pos_ex, pos_ey, vel_ex, vel_ey, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_ix, pos_iy, vel_ix, vel_iy, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tsmot += clock() - tiempo;

    if(compare(pos_e_x, pos_ex, MAX_SPE))
      cout << "fail posex" << endl;
    if(compare(pos_e_y, pos_ey, MAX_SPE))
      cout << "fail posey" << endl;
    if(compare(pos_i_x, pos_ix, MAX_SPE))
      cout << "fail posix" << endl;
    if(compare(pos_i_y, pos_iy, MAX_SPE))
      cout << "fail posiy" << endl;

    if(compare(vel_e_x, vel_ex, MAX_SPE))
      cout << "fail velex" << endl;
    if(compare(vel_e_y, vel_ey, MAX_SPE))
      cout << "fail veley" << endl;
    if(compare(vel_i_x, vel_ix, MAX_SPE))
      cout << "fail velix" << endl;
    if(compare(vel_i_y, vel_iy, MAX_SPE))
      cout << "fail veliy" << endl;
#endif
  } //Cierre del ciclo principal
  int div = max_it * CLOCKS_PER_SEC;
  cout << std::fixed;
  cout << "Concentration\nGPU = \n\t" << tcon / div << "\nsec  CPU = \n\t" << tscon / div << endl;
  cout << "Electric field\nGPU = \n\t" << telec / div << "\nsec  CPU = \n\t" << tselec / div << endl;
  cout << "Motion\nGPU =\n\t" << tmot / div << "\nsec CPU =\n\t" << tsmot / div << endl;
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

