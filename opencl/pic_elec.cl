__kernel void D_electric_field (__global double *d_phi,
                                __global double *d_E_X,
                                __global double *d_E_Y,
                                double hx, int J_X, int J_Y) {
  int j = get_global_id(0);
  int k = get_global_id(1);
  if(j < J_X && k < J_Y) {
    d_E_X[j * J_Y + k] = (d_phi[(j - 1) * J_Y + k] -
                         d_phi[(j + 1) * J_Y + k]) / (2. * hx);
    d_E_Y[j * J_Y + k] = (d_phi[j * J_Y + ((J_Y + k - 1) % J_Y)]
                         - d_phi[j * J_Y + ((k + 1) % J_Y)])
                         / (2. * hx);
  }
}
