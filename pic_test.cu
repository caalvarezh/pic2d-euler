#include "pic.cpp"
#include "pic.cu"

using namespace std;

bool test_concentration(int test_num, bool parallel = true) {
  char name[40];
  bool flag = true;
  sprintf(name, "test_concentration%d.data", test_num);
  freopen (name, "r", stdin);
  cout << "Reading from " << name << endl;
  int pos_size;
  cin >> pos_size;
  double *pos;
  pos = (double *) malloc(pos_size * sizeof(double));
  for(int i = 0; i < pos_size; i++)
    cin >> pos[i];
  int n_size;
  cin >> n_size;
  double *n_in, *n_out;
  n_in = (double *) malloc(n_size * sizeof(double));
  n_out = (double *) malloc(n_size * sizeof(double));
  for(int i = 0; i < n_size; i++)
    cin >> n_in[i];
  int NSP;
  double hx;
  cin >> NSP >> hx;
  cout << "..." << endl;
  if(parallel) {
    pic_cuda::H_Concentration(pos, n_out, NSP, hx);

    for(int i = 0; i < n_size; i++) {
      if(n_in[i] - n_out[i] > 10e-3) {
        cout << "Error at " << i << " : " << n_in[i] << " " << n_out[i] << endl;
        flag = false;
      }
    }
  } else {
    pic::Concentration(pos,n_out,NSP,hx);
  }
  fclose(stdin);
  return flag;
}

bool test_motion(int test_num, bool parallel = true) {
  char name[40];
  bool flag = true;
  sprintf(name, "test_Motion%d.data", test_num);
  freopen (name, "r", stdin);
  cout << "Reading from " << name << endl;

  int pos_size;
  cin >> pos_size;
  double *pos;
  pos = (double *) malloc(pos_size * sizeof(double));
  for(int i = 0; i < pos_size; i++)
    cin >> pos[i];

  int vel_size;
  cin >> vel_size;
  double *vel;
  vel = (double *) malloc(vel_size * sizeof(double));
  for(int i = 0; i < vel_size; i++)
    cin >> vel[i];

  int NSP, especie;
  cin >> NSP >> especie;

  int E_X_size;
  cin >> E_X_size;
  double *E_X;
  E_X = (double *) malloc(E_X_size * sizeof(double));
  for(int i = 0; i < E_X_size; i++)
    cin >> E_X[i];

  int E_Y_size;
  cin >> E_Y_size;
  double *E_Y;
  E_Y = (double *) malloc(E_Y_size * sizeof(double));
  for(int i = 0; i < E_Y_size; i++)
    cin >> E_Y[i];

  double hx, mv2perdidas;
  int total_perdidos;
  cin >> hx >> total_perdidos >> mv2perdidas;

  fclose(stdin);
  cout << "..." << endl;
  if(parallel) {
    pic_cuda::H_Motion(pos, vel, NSP, especie, E_X, E_Y, hx, total_perdidos, mv2perdidas);
    /*sprintf(name, "test_out_Motion%d.data", test_num);
    freopen (name, "r", stdin);
    //cout << "Reading from " << name << endl;

    cin >> pos_size;
    double out_pos;
    for(int i = 0; i < pos_size; i++) {
      cin >> out_pos;
      if(abs(pos[i] - out_pos) > 10e-3)
        cout << "Error in pos at " << i << " with " << pos[i] << " != " << out_pos << endl;
    }

    cin >> vel_size;
    double out_vel;
    for(int i = 0; i < vel_size; i++) {
      cin >> out_vel;
      if(abs(vel[i] - out_vel) > 10e-3)
        cout << "Error in vel at " << i << " with " << vel[i] << " != " << out_vel << endl;
    }

    int n_NSP;
    cin >> n_NSP;
    if(NSP != n_NSP)
      cout << "Error in NSP " << NSP << " " << n_NSP << endl;

    double out_mv2perdidas;
    int out_total_perdidos;
    cin >> out_total_perdidos >> out_mv2perdidas;
    if(out_total_perdidos != total_perdidos)
      cout << "Error in total_perdidos " << total_perdidos << " " << out_total_perdidos << endl;
    if(abs(mv2perdidas - out_mv2perdidas) > 10e-3)
      cout << "Error in mv2perdidas " << mv2perdidas << " " << out_mv2perdidas << endl;
    fclose(stdin);*/

  } else {
    pic::Motion(pos, vel, NSP, especie, E_X, E_Y, hx, total_perdidos, mv2perdidas);
  }
  free(pos);
  free(vel);
  free(E_X);
  free(E_Y);
  return flag;
}


int main() {
  int max_test = 64;
  clock_t t = clock();
  for(int i = 0; i < max_test; i++) {
    cout << "Running test " << i << ":" << endl;
    test_motion(i);
      //cout << "ERROR IN PARALLEL" << endl;
      //cout << "Test " << i << " passed." << endl;
  }
  t = clock() - t;
  double a = ((float)t)/CLOCKS_PER_SEC;
  cout << "Parallel: " << a << endl;
  t = clock();
  for(int i = 0; i < max_test; i++) {
    cout << "Running test Sec " << i << ":" << endl;
    test_motion(i, false);
      //cout << "Test " << i << " passed." << endl;
  }
  t = clock() - t;
  a = ((float)t)/CLOCKS_PER_SEC;
  cout << "Secuential: " << a << endl;

  return 0;
}
