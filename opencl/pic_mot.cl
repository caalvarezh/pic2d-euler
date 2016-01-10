__kernel
void D_Motion(__global double *pos_x, __global double *pos_y,
        __global double *vel_x, __global double *vel_y, int NSP, double fact,
        __global double *E_X, __global double *E_Y, double hx, double L_MAX_X,
        double L_MAX_Y, double DT, int J_X, int J_Y) {
  int j_x, j_y;
  double temp_x, temp_y, Ep_X, Ep_Y;
  double jr_x, jr_y;
  int i = get_global_id(0);
  if ( i < NSP) {
    jr_x = pos_x[i] / hx;     // Índice (real) de la posición de la superpartícula (X)
    j_x  = ((int) jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
    temp_x = jr_x - ((double) j_x);
    jr_y = pos_y[i] / hx;     // Índice (real) de la posición de la superpartícula (Y)
    j_y  = ((int) jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
    temp_y  =  jr_y - ((double) j_y);

    Ep_X = (1 - temp_x) * (1 - temp_y) * E_X[j_x * J_Y + j_y] +
      temp_x * (1 - temp_y) * E_X[(j_x + 1) * J_Y + j_y] +
      (1 - temp_x) * temp_y * E_X[j_x * J_Y + (j_y + 1)] +
      temp_x * temp_y * E_X[(j_x + 1) * J_Y + (j_y + 1)];

    Ep_Y = (1 - temp_x) * (1 - temp_y) * E_Y[j_x * J_Y + j_y] +
      temp_x * (1 - temp_y) * E_Y[(j_x + 1) * J_Y + j_y] +
      (1 - temp_x) * temp_y * E_Y[j_x * J_Y + (j_y + 1)] +
      temp_x * temp_y * E_Y[(j_x + 1) * J_Y + (j_y + 1)];


    vel_x[i] += DT * fact * Ep_X;
    vel_y[i] += DT * fact * Ep_Y;

    pos_x[i] += vel_x[i] * DT;
    pos_y[i] += vel_y[i] * DT;

    if(pos_x[i] < 0) {//Rebote en la pared del material.
      pos_x[i] = -pos_x[i];
      vel_x[i] = -vel_x[i];
    }

    while(pos_y[i] > L_MAX_Y) //Ciclo en el eje Y.
      pos_y[i] = pos_y[i] - L_MAX_Y;

    while(pos_y[i] < 0.0) //Ciclo en el eje Y.
      pos_y[i] = L_MAX_Y + pos_y[i];
  }
}

