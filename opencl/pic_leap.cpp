#include "pic.cpp"

using namespace std;
using namespace pic_cl;

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

inline void checkErr(cl_int err, const char * name) {
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name  << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
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

  /*
  ** OpenCL precompute
  */
  // Query for platforms
  std::vector < cl::Platform > platforms;
  cl::Platform::get(&platforms);
  // Get a list of devices on this platform
  std::vector < cl::Device > devices;
  platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
  // Create a context for the devices 
  cl::Context context(devices);
  // Create a command queue for the first device
  cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);
  // Read the program source
  std::ifstream sourceFile("pic_cl.cl");
  std::string sourceCode( std::istreambuf_iterator < char > (sourceFile), (std::istreambuf_iterator < char > ()));
  cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length() + 1));
  // Make program from the source code
  cl::Program program = cl::Program(context, source);
  // Build the program for the devices
  program.build(devices);
  /*
  ** End OpenCL
  */

  //***************************
  // Normalización de variables
  //***************************

  double hx = DELTA_X / X0;                            // Paso espacial
  int max_it = 100000;

  initialize_vectors (pos_e_x, pos_e_y, vel_e_x, vel_e_y, E_X);
  initialize_vectors (pos_i_x, pos_i_y, vel_i_x, vel_i_y, E_Y);
  for(int i = 0; i < J_X * J_Y; i++)
    phi[i] =  rand() % 8234;

  double telec = 0.0, tmot = 0.0 ;
  clock_t tiempo;
  for(int it  =  0; it <= max_it; it++) {

    // Calcular campo eléctrico en puntos de malla
    tiempo = clock();
    H_electric_field(phi, E_X, E_Y, hx, context, program, queue);
    telec += clock() - tiempo;

    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    tiempo = clock();
    H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx,
        total_e_perdidos, mv2perdidas, context, program, queue);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx,
        total_i_perdidos, mv2perdidas, context, program, queue);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tmot += clock() - tiempo;

  } //Cierre del ciclo principal
  cout << " CPU time Electric field =  " << telec / CLOCKS_PER_SEC << " sec" << endl;
  cout << " CPU time Motion         =  " << tmot / CLOCKS_PER_SEC << " sec" << endl;
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

