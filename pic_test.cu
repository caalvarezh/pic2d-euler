#include "pic_cuda.cu"

using namespace std;

double test_electric (int test_num, bool parallel = true) {
  char name[40];
  sprintf(name, "test_electric%d.data", test_num);
  freopen (name, "r", stdin);
  //cout << "Reading from " << name << endl;
  int size;
  cin >> size;
  double *phi = (double *) malloc(size * sizeof(double));
  for(int i = 0; i < size; i++)
    cin >> phi[i];
  double *E_X = (double *) malloc(size * sizeof(double));
  for(int i = 0; i < size; i++)
    cin >> E_X[i];
  double *E_Y = (double *) malloc(size * sizeof(double));
  for(int i = 0; i < size; i++)
    cin >> E_Y[i];
  int hx;
  cin >> hx;
  clock_t t = clock();
  if(parallel){
    pic_cuda::H_electric_field(phi, E_X, E_Y, hx);
  }else{
    pic_cuda::electric_field(phi, E_X, E_Y, hx);
  }
  t = clock() - t;
  return (double)t;
}

double test_concentration(int test_num, bool parallel = true) {
  char name[40];
  sprintf(name, "test_concentration%d.data", test_num);
  freopen (name, "r", stdin);
  //cout << "Reading from " << name << endl;
  int pos_size;
  int J_X = pic_cuda::J_X;
  int J_Y = pic_cuda::J_Y;
  cin >> pos_size;
  double *pos_x, *pos_y;
  pos_x = (double *) malloc(pos_size * sizeof(double));
  pos_y = (double *) malloc(pos_size * sizeof(double));
  for(int i = 0; i < pos_size; i++)
    cin >> pos_x[i];
  for(int i = 0; i < pos_size; i++)
    cin >> pos_y[i];
  int NSP = pos_size;
  double hx;
  double *n_out = (double *) malloc( J_X * J_Y * sizeof(double));
  cin >>  hx;
  //cout << "..." << endl;
  clock_t t;
  if(parallel) {
    t = clock();
    pic_cuda::H_Concentration(pos_x, pos_y, n_out, NSP, hx);
    t = clock() - t;
    /*for(int i = 0; i < n_size; i++) {
      if(n_in[i] - n_out[i] > 10e-3) {
        cout << "Error at " << i << " : " << n_in[i] << " " << n_out[i] << endl;
        flag = false;
      }
    }*/
  } else {
    t = clock();
    pic_cuda::Concentration(pos_x, pos_y, n_out, NSP, hx);
    t = clock() - t;
  }
  fclose(stdin);
  return (double)t;
}

double test_motion(int test_num, bool parallel = true) {
  char name[40];
  sprintf(name, "test_Motion%d.data", test_num);
  freopen (name, "r", stdin);
  //cout << "Reading from " << name << endl;

  int pos_size;
  cin >> pos_size;
  double *pos_x, *pos_y;
  double *vel_x, *vel_y;
  pos_x = (double *) malloc(pos_size * sizeof(double));
  pos_y = (double *) malloc(pos_size * sizeof(double));
  vel_x = (double *) malloc(pos_size * sizeof(double));
  vel_y = (double *) malloc(pos_size * sizeof(double));
  for(int i = 0; i < pos_size; i++)
    cin >> pos_x[i];

  for(int i = 0; i < pos_size; i++)
    cin >> pos_y[i];

  for(int i = 0; i < pos_size; i++)
    cin >> vel_x[i];
  for(int i = 0; i < pos_size; i++)
    cin >> vel_y[i];

  int NSP = pos_size, especie;
  cin >> especie;

  int E_X_size;
  cin >> E_X_size;
  double *E_X;
  E_X = (double *) malloc(E_X_size * sizeof(double));
  for(int i = 0; i < E_X_size; i++)
    cin >> E_X[i];

  double *E_Y;
  E_Y = (double *) malloc(E_X_size * sizeof(double));
  for(int i = 0; i < E_X_size; i++)
    cin >> E_Y[i];

  double hx, mv2perdidas = 0;
  int total_perdidos = 0;
  cin >> hx;

  fclose(stdin);
  clock_t t = clock();
  if(parallel) {
    pic_cuda::H_Motion(pos_x, pos_y, vel_x, vel_y, NSP, especie, E_X, E_Y, hx, total_perdidos, mv2perdidas);
  } else {
    pic_cuda::Motion(pos_x, pos_y, vel_x, vel_y, NSP, especie, E_X, E_Y, hx, total_perdidos, mv2perdidas);
  }
  t = clock() - t;
  free(pos_x);
  free(pos_y);
  free(vel_x);
  free(vel_y);
  free(E_X);
  free(E_Y);
  return double(t);
}


int main() {
  int max_test = 25;
  clock_t t = clock();
  double time = 0;
  for(int i = 0; i < max_test; i++) {
    //cout << "Running test Sec " << i << ":" << endl;
    time += test_concentration(i, false);
    //time += test_electric(i, false);
    //time += test_motion(i, false);
    //cout << "Test " << i << " passed." << endl;
  }
  double b = time /CLOCKS_PER_SEC;
  cout << "Secuential: " << b << endl;
  cout << "--------------------------------------------------" << endl;
  time = 0;
  for(int i = 0; i < max_test; i++) {
    //cout << "Running test " << i << ":" << endl;
    time += test_concentration(i);
    //time += test_electric(i);
    //time += test_motion(i);
  }
  double a = time / CLOCKS_PER_SEC;
  cout << "Parallel: " << a << endl;

  cout << "acelerated " << b / a << endl;

  return 0;
}
