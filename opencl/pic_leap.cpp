#include "pic_cuda.cu"

using namespace std;
using namespace pic_cuda;

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
  ** OpenCL Device Query
  */

// Use this to check the output of each API call
    cl_int status;

    //-----------------------------------------------------
    // STEP 1: Discover and initialize the platforms
    //-----------------------------------------------------

    cl_uint numPlatforms = 0;
    cl_platform_id *platforms = NULL;

    // Use clGetPlatformIDs() to retrieve the number of
    // platforms
    status = clGetPlatformIDs(0, NULL, &numPlatforms);

    // Allocate enough space for each platform
    platforms =
        (cl_platform_id*)malloc(
                numPlatforms*sizeof(cl_platform_id));

    // Fill in platforms with clGetPlatformIDs()
    status = clGetPlatformIDs(numPlatforms, platforms,
            NULL);

    //-----------------------------------------------------
    // STEP 2: Discover and initialize the devices
    //-----------------------------------------------------

    cl_uint numDevices = 0;
    cl_device_id *devices = NULL;

    // Use clGetDeviceIDs() to retrieve the number of
    // devices present
    status = clGetDeviceIDs(
            platforms[0],
            CL_DEVICE_TYPE_ALL,
            0,
            NULL,
            &numDevices);

    // Allocate enough space for each device
    devices =
        (cl_device_id*)malloc(
                numDevices*sizeof(cl_device_id));

    // Fill in devices with clGetDeviceIDs()
    status = clGetDeviceIDs(
            platforms[0],
            CL_DEVICE_TYPE_ALL,
            numDevices,
            devices,
            NULL);

    //-----------------------------------------------------
    // STEP 3: Create a context
    //-----------------------------------------------------

    cl_context context = NULL;

    // Create a context using clCreateContext() and
    // associate it with the devices
    context = clCreateContext(
            NULL,
            numDevices,
            devices,
            NULL,
            NULL,
            &status);

    //-----------------------------------------------------
    // STEP 4: Create a command queue
    //-----------------------------------------------------

    cl_command_queue cmdQueue;

    // Create a command queue using clCreateCommandQueue(),
    // and associate it with the device you want to execute
    // on
    cmdQueue = clCreateCommandQueue(
            context,
            devices[0],
            0,
            &status);




  /*
  ** End Device Query
  */

  //***************************
  // Normalización de variables
  //***************************

  double hx = DELTA_X / X0;                            // Paso espacial
  int max_it = 1000;

  initialize_vectors (pos_e_x, pos_e_y, vel_e_x, vel_e_y, E_X);
  initialize_vectors (pos_i_x, pos_i_y, vel_i_x, vel_i_y, E_Y);
  for(int i = 0; i < J_X * J_Y; i++)
    phi[i] =  rand() % 8234;

  double tcon = 0.0, telec = 0.0, tmot = 0.0 ;
  clock_t tiempo;
  for(int it  =  0; it <= max_it; it++) {

    // Calculo de "densidad de carga 2D del plasma"
    tiempo = clock();
    H_Concentration (pos_e_x, pos_e_y, ne, le, hx);// Calcular concentración de superpartículas electrónicas
    H_Concentration (pos_i_x, pos_i_y, ni, li, hx);// Calcular concentración de superpartículas Iónicas
    tcon += clock() - tiempo;

    // Calcular campo eléctrico en puntos de malla
    tiempo = clock();
    H_electric_field(phi, E_X, E_Y, hx);
    telec += clock() - tiempo;

    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    tiempo = clock();
    H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tmot += clock() - tiempo;

  } //Cierre del ciclo principal
  cout << " CPU time Concentration  =  " << tcon / CLOCKS_PER_SEC << " sec" << endl;
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

