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

int main() {
  int max_test = 64;
  clock_t t = clock();
  for(int i = 0; i < max_test; i++) {
    cout << "Running test " << i << ":" << endl;
    if(test_concentration(i))
      //cout << "ERROR IN PARALLEL" << endl;
      cout << "Test " << i << " passed." << endl;
  }
  /*t = clock() - t;
  double a = ((float)t)/CLOCKS_PER_SEC;
  cout << "Parallel: " << a << endl;
  t = clock();
  for(int i = 0; i < max_test; i++) {
    cout << "Running test Sec " << i << ":" << endl;
    test_concentration(i, false);
      //cout << "Test " << i << " passed." << endl;
  }
  t = clock() - t;
  a = ((float)t)/CLOCKS_PER_SEC;
  cout << "Secuential: " << a << endl;
*/
  return 0;
}
