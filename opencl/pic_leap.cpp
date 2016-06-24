#include "pic.cpp"

using namespace std;
using namespace pic_cl;

bool compare(double *c1, double *c2, int size) {
  bool flag = false;
  int cont = 0;
  for(int i = 0; i < size; i++) {
    if(fabs(c1[i] - c2[i]) > 1) {
      cout << i << " : " << c1[i] << " -  " << c2[i] << " = " << fabs(c1[i] - c2[i])<< endl;
      //cont++;
      flag = true;
      return flag;
    }
  }
  if(cont > 0)
    cout << "fail in " << cont << " of " << size << " = " << size - cont << endl;
  return flag;
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
  
                       //Parámetro del plasma
 // FILE    *outEnergia;


  //***************************
  //Constantes de normalización
  //***************************

  int size = MAX_SPE * sizeof(double);
  int size1 = J_X * J_Y * sizeof(double);

  double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y, *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y,
         *phi, *E_X, *E_Y, *E_X1, *E_Y1, *pos_ex, *pos_ey, *vel_ex, *vel_ey,
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
  phi   = (double *) malloc(size1);
  E_X   = (double *) malloc(size1);
  E_Y   = (double *) malloc(size1);
  phi   = (double *) malloc(size1);

  pos_ex = (double *) malloc(size);
  pos_ey = (double *) malloc(size);
  pos_ix = (double *) malloc(size);
  pos_iy = (double *) malloc(size);

  vel_ex = (double *) malloc(size);
  vel_ey = (double *) malloc(size);
  vel_ix = (double *) malloc(size);
  vel_iy = (double *) malloc(size);
  E_X1  = (double *) malloc(size1);
  E_Y1  = (double *) malloc(size1);

  /*
  ** OpenCL precompute
  */
  // Query for platforms
  cl_int error;
  vector < cl::Platform > platforms;
  error = cl::Platform::get(&platforms);
  checkErr(error, "Platform");
  // Get a list of devices on this platform
  vector < cl::Device > devices;
  error = platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
  checkErr(error, "device");
  cl::Device default_device = devices[0];
  // Create a context for the devices
  cl::Context context(devices, NULL, NULL, NULL, &error);
  checkErr(error, "context");
  // Create a command queue for the first device
  cl::CommandQueue queue = cl::CommandQueue(context, default_device, 0, &error);
  checkErr(error, "command queue");
  // Read the program source
  cl::Program::Sources source;
  ifstream sourceFile;

  sourceFile.open("pic.cl");
  string s_elec( istreambuf_iterator < char > (sourceFile), (istreambuf_iterator < char > ()));
  sourceFile.close();

  source.push_back(make_pair(s_elec.c_str(),   s_elec.length() + 1) );
  //checkErr(error, "program source");

  // Make program from the source code
  cl::Program program = cl::Program(context, source, &error);
  checkErr(error, "program");

  // Build the program for the devices
  error = program.build(devices);
  //cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";

  checkErr(error, "program build");
  /*
  ** End OpenCL
  */

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

  initialize_Particles (pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, FE_MAXWELL_X, FE_MAXWELL_Y, VPHI_E_X, VPHI_E_Y);//Velocidades y posiciones iniciales de las partículas>
  initialize_Particles (pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, FI_MAXWELL_X, FI_MAXWELL_Y, VPHI_I_X, VPHI_I_Y);//Velocidades y posiciones iniciales de las partículas>

  double tacum = 0;
  for(int kk  =  0, kt  =  0; kt <= K_total; kt++) {
    cout << kt << endl;
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

    for (int j  =  0; j < J_X; j++)
      for (int k  =  0; k < J_Y; k++)
        rho[j * J_Y + k] = cte_rho * FACTOR_CARGA_E * (ni[j * J_Y + k] - ne[j * J_Y + k]) / n_0;

    // Calcular potencial eléctrico en puntos de malla
    poisson2D_dirichletX_periodicY(phi, rho, hx);
    // Calcular campo eléctrico en puntos de malla

    electric_field(phi, E_X, E_Y, hx);

    // imprimir el potencial electroestatico.
    if(kt % 5000  ==  0) {
      cout << "Poisson #" << kt << endl; 
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

    Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);

    clock_t tiempo1  =  clock();
    if(kt % 5000 == 0) {
      cout << " CPU time " << kt / 5000 << "  =  " << double(tiempo1 - tiempo0) / CLOCKS_PER_SEC << " sec" << endl;
      tiempo0  =  clock();
    }
} //Cierre del ciclo principal
/*
  double hx = DELTA_X / X0;                            // Paso espacial
  int max_it = 0;

  double telec, tselec, tmot, tsmot;
  telec = tselec = tmot = tsmot = 0.0 ;
  clock_t tiempo;
  for(int it  =  0; it <= max_it; it++) {
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

    // Calcular campo eléctrico en puntos de malla
    tiempo = clock();
    H_electric_field(phi, E_X, E_Y, hx, context, program, queue);
    telec += clock() - tiempo;

    tiempo = clock();
    electric_field(phi, E_X1, E_Y1, hx);
    tselec += clock() - tiempo;

    if(compare(E_X, E_X1, J_X * J_Y))
      cout << "fail ex" << endl;
    if(compare(E_Y, E_Y1, J_X * J_Y))
      cout << "fail ey" << endl;

    // Avanzar posiciones de superpartículas electrónicas e Iónicas
    tiempo = clock();
    cout << "S ELEC" << endl;
    H_Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y, le, ELECTRONS, E_X, E_Y, hx,
        total_e_perdidos, mv2perdidas, context, program, queue);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    cout << "S ION" << endl;
    H_Motion(pos_i_x, pos_i_y, vel_i_x, vel_i_y, li, IONS, E_X, E_Y, hx,
        total_i_perdidos, mv2perdidas, context, program, queue);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tmot += clock() - tiempo;

    tiempo = clock();
    cout << "H ELEC" << endl;
    Motion(pos_ex, pos_ey, vel_ex, vel_ey, le, ELECTRONS, E_X, E_Y, hx, total_e_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    cout << "H ION" << endl;
    Motion(pos_ix, pos_iy, vel_ix, vel_iy, li, IONS, E_X, E_Y, hx, total_i_perdidos, mv2perdidas);//, total_elec_perdidos, total_ion_perdidos, mv2_perdidas);
    tsmot += clock() - tiempo;

    if(compare(pos_e_x, pos_ex, MAX_SPE))
      cout << "fail posex" << endl;
    else
      cout << "success posex" << endl;
    if(compare(pos_e_y, pos_ey, MAX_SPE))
      cout << "fail posey" << endl;
    else
      cout << "success posey" << endl;
    if(compare(vel_e_x, vel_ex, MAX_SPE))
      cout << "fail velex" << endl;
    else
      cout << "success velex" << endl;
    if(compare(vel_e_y, vel_ey, MAX_SPE))
      cout << "fail veley" << endl;
    else
      cout << "success veley" << endl;

    if(compare(pos_i_x, pos_ix, MAX_SPE))
      cout << "fail posix" << endl;
    else
      cout << "success posix" << endl;
    if(compare(pos_i_y, pos_iy, MAX_SPE))
      cout << "fail posiy" << endl;
    else
      cout << "success posiy" << endl;

    if(compare(vel_i_x, vel_ix, MAX_SPE))
      cout << "fail velix" << endl;
    else
      cout << "success velix" << endl;
    if(compare(vel_i_y, vel_iy, MAX_SPE))
      cout << "fail veliy" << endl;
    else
      cout << "success veliy" << endl;
  } //Cierre del ciclo principal*/
  cout << fixed;
  /*cout << " GPU time Electric field =  " << telec / CLOCKS_PER_SEC << " sec" << endl;
  cout << " CPU time Electric field =  " << tselec / CLOCKS_PER_SEC << " sec" << endl;
  cout << " Electric field X = " << tselec / telec << endl;
  cout << " GPU time Motion         =  " << tmot / CLOCKS_PER_SEC << " sec" << endl;
  cout << " CPU time Motion         =  " << tsmot / CLOCKS_PER_SEC << " sec" << endl;
  cout << " Motion X = " << tsmot / tmot << endl;*/
  cout << "Termino" << endl;
  free(pos_e_x);
  free(pos_e_y);
  free(pos_i_x);
  free(pos_i_y);
  free(vel_e_x);
  free(vel_e_y);
  free(vel_i_x);
  free(vel_i_y);
  free(phi);
  free(E_X);
  free(E_Y);

  return (0);
}// FINAL MAIN

