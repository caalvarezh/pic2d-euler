__kernel void D_electric_field_border (__global double *d_phi,
                                       __global double *d_E_X,
                                       __global double *d_E_Y,
                                       double hx, int J_X,
                                       int J_Y) {
  int k = get_global_id(0);
  if(k < J_Y) {
    d_E_X[k] = 0.0; //Cero en las fronteras X
    d_E_Y[k] = 0.0;
    d_E_X[(J_X - 1) * J_Y + k] = 0.0;
    d_E_Y[(J_X - 1) * J_Y + k] = 0.0;
  }
}
