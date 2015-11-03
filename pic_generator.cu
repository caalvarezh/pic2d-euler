#include <bits/stdc++.h>
#include "pic_cuda.cu"

using namespace std;
using namespace pic_cuda;

void write_test_concentration() {
  static int num = 0;
  char name[40];
  sprintf(name, "test_concentration%d.data", num++);
  freopen (name, "w", stdout);
  int NSP = MAX_SPE;
  double hx = 0.983064;
  cout << NSP << endl;
  for(int i = 0; i < NSP; i++)
    cout << rand() % int(L_MAX_X) << " ";
  for(int i = 0; i < NSP; i++)
    cout << rand() % int(L_MAX_Y) << " ";
  cout << endl << hx << endl;
  fclose(stdout);
}

void write_test_electric() {
  static int num = 0;
  char name[40];
  sprintf(name, "test_electric%d.data", num++);
  freopen (name, "w", stdout);
  double hx = 0.983064;
  cout << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rand() % 8234 << " ";
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rand() % 2342 << " ";
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rand() % 2756 << " ";
  cout << endl << hx << endl;
  fclose(stdout);
}

void write_test_Motion(int especie) {
  static int num = 0;
  char name[40];
  sprintf(name, "test_Motion%d.data", num++);
  freopen (name, "w", stdout);

  int NSP = MAX_SPE;
  double hx = 0.983064;
  cout << NSP << endl;
  for(int i = 0; i < NSP; i++)
    cout << rand() % int(L_MAX_X) << " ";
  cout << endl;
  for(int i = 0; i < NSP; i++)
    cout << rand() % int(L_MAX_Y) << " ";
  cout << endl;
  for(int i = 0; i < NSP; i++)
    cout << (rand() % 100) / 2.0 << " ";
  cout << endl;
  for(int i = 0; i < NSP; i++)
    cout << (rand() % 100) / 2.0 << " ";

  cout << endl << especie << endl << J_X * J_Y << endl;
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rand() % 2342 << " ";
  for(int i = 0; i < J_X * J_Y; i++)
    cout << rand() % 2756 << " ";
  cout << endl << hx << endl;
  fclose(stdout);
}

int main(){
  int test = 26;
  for(int i = 0; i < test; i++) {
    write_test_concentration();
  }
  return 0;
}

