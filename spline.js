"use strict";

class CubicSplineX {
  constructor(control_values, control_t, N_ctr, is_periodic = false) {
    //let n;
    //if (!is_periodic)
    //  n = N_ctr - 1;
    //else
    //  n = N_ctr;

    if (is_periodic) {
        this.h = new Array(N_ctr); // длины отрезков между к.т. (по параметру)
        for (let i = 0; i < N_ctr; ++i) {
          if (i + 1 == N_ctr)
            this.h[i] = 1 - control_t[i];
          else
            this.h[i] = control_t[i + 1] - control_t[i];
        }

        this.delta = new Array(N_ctr - 1);
        for (let i = 0; i < N_ctr - 1; ++i) {
            this.delta[i] = ( this.control_values[i + 2] - this.control_values[i + 1] ) / this.h[i + 1] - 
                            ( this.control_values[i + 1] - this.control_values[i] ) / this.h[i];
        }

    } else {
        this.h = new Array(N_ctr - 1); // длины отрезков между к.т. (по параметру)
        for (let i = 0; i < N_ctr - 1; ++i) {
          this.h[i] = control_t[i + 1] - control_t[i];
        }

        this.delta = new Array(N_ctr - 2);
        for (let i = 0; i < N_ctr - 2; ++i) {
            this.delta[i] = ( this.control_values[i + 2] - this.control_values[i + 1] ) / this.h[i + 1] - 
                            ( this.control_values[i + 1] - this.control_values[i] ) / this.h[i];
        }
    }

    if (is_periodic) {
        this.parameters_p = new Array(N_ctr + 1);
        this.parameters_M = new Array(N_ctr + 1);
    } else {
        this.parameters_p = new Array(N_ctr);
        this.parameters_M = new Array(N_ctr);
    }

    for (let i = 0; i < N_ctr; ++i)
      this.parameters_p[i] = control_values[i];
    if (is_periodic)
        this.parameters_p[N_ctr] = this.parameters_p[0];
    
    if (is_periodic) {
        let matrix = new Array(N_ctr),
            free = new Array(N_ctr);
        for (let i = 0; i < N_ctr; ++i)
          matrix[i] = new Array(N_ctr);
        zeros(matrix, N_ctr, N_ctr);

        for (let i = 0; i < N_ctr; ++i) {
            if (i == 0) {
                matrix[i][0] = 2 * (this.h[0] + this.h[1]);
                matrix[i][1] = this.h[1];
                matrix[i][N_ctr - 1] = this.h[0];

                free[i] = 6 * this.delta[0];
            } else if (i == N_ctr - 1) {
                matrix[i][N_ctr - 1] = 2 * (this.h[0] + this.h[1]);
                matrix[i][N_ctr - 2] = this.h[N_ctr - 1];
                matrix[i][0] = this.h[0];

                free[i] = 6 * ( (this.control_values[1] - this.control_values[0]) / this.h[0] - (this.control_values[0] - this.control_values[N_ctr - 1]) / this.h[N_ctr - 1] );
            } else {
                matrix[i][i] = 2 * (this.h[i] + this.h[i + 1]);
                matrix[i][i - 1] = this.h[i];
                matrix[i][i + 1] = this.h[i + 1];

                free[i] = 6 * this.delta[i];
            }
        }

        let tmp = solve(matrix, free);
        for (let i = 0; i < N_ctr; ++i) {
            this.parameters_M[i] = tmp[i] / 6; // bar M
        }

        this.parameters_M[N_ctr] = this.parameters_M[0];
    } else {
        let matrix = new Array(N_ctr),
            free = new Array(N_ctr);
        for (let i = 0; i < N_ctr; ++i)
          matrix[i] = new Array(N_ctr);
        zeros(matrix, N_ctr, N_ctr);

        for (let i = 0; i < N_ctr; ++i) {
            if (i == 0) {
                matrix[i][0] = this.h[1];
                matrix[i][1] = - (this.h[0] + this.h[1]);
                matrix[i][2] = this.h[0];

                free[i] = 0;
            } else if (i == N_ctr - 1) {
                matrix[i][N_ctr - 3] = this.h[N_ctr - 1];
                matrix[i][N_ctr - 2] = - ( this.h[N_ctr - 2] + this.h[N_ctr - 1] );
                matrix[i][N_ctr - 1] = this.h[N_ctr - 2];

                free[i] = 0;
            } else {
                matrix[i][i - 1] = this.h[i - 1];
                matrix[i][i] = 2 * (this.h[i - 1] + this.h[i]);
                matrix[i][i + 1] = this.h[i];

                free[i] = 6 * this.delta[i - 1];
            }
        }

        let tmp = solve(matrix, free);
        for (let i = 0; i < N_ctr; ++i) {
            this.parameters_M[i] = tmp[i] / 6; // bar M
        }
    }
  }

  calc_value(omega, i) {
    let omega1 = 1 - omega, // psi
      omega2 = omega,
      omega3 = omega * (omega - 1) * (2 - omega),
      omega4 = omega * (omega * omega - 1);
    let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
    return coefs[0] * omega1 + coefs[1] * omega2 +
      coefs[2] * this.h[i] * this.h[i] * omega3 +
      coefs[3] * this.h[i] * this.h[i] * omega4;
  }

  calc_diff(omega, i) {
    let omega1 = - 1. / this.h[i],
      omega2 = 1. / this.h[i],
      omega3 = (6 * omega - 3 * omega * omega - 2) / this.h[i],
      omega4 = (3 * omega * omega - 1) / this.h[i];
    let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
    return coefs[0] * omega1 + coefs[1] * omega2 +
      coefs[2] * this.h[i] * this.h[i] * omega3 +
      coefs[3] * this.h[i] * this.h[i] * omega4;
  }

  calc_second_diff(omega, i) {
    let omega1 = 0.,
      omega2 = 0.,
      omega3 = (6 - 6 * omega) / (this.h[i] * this.h[i]),
      omega4 = (6 * omega) / (this.h[i] * this.h[i]);
    let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
    return coefs[0] * omega1 +
      coefs[1] * omega2 +
      coefs[2] * this.h[i] * this.h[i] * omega3 +
      coefs[3] * this.h[i] * this.h[i] * omega4;
  }
}


class BicubicSplineX {
  constructor(control_values, control_u, control_v, N_ctr, M_ctr, derivative) {
    const n = N_ctr - 1;
    const m = M_ctr; // Так как мы хотим замкнуть цилиндр вдоль полярной координаты

    this.u_splines = new Array(m); // S(t, tau_l) = S(u, v_j)
    for (let j = 0; j < M_ctr; ++j) {
      let tmp_control_values = new Array(N_ctr),
        tmp_control_u = new Array(N_ctr);
      for (let i = 0; i < N_ctr; ++i) {
        tmp_control_values[i] = control_values[i][j];
        tmp_control_u[i] = control_u[i][j];
      }
      this.u_splines[j] = new CubicSplineX(tmp_control_values, tmp_control_u, N_ctr);
    }

    let v_splines = new Array(N_ctr);
    for (let i = 0; i < N_ctr; ++i) {
      let tmp_control_values = new Array(M_ctr),
        tmp_control_v = new Array(M_ctr);
      for (let j = 0; j < M_ctr; ++j) {
        tmp_control_values[j] = control_values[i][j];
        tmp_control_v[j] = control_v[i][j];
      }
      v_splines[i] = new CubicSplineX(tmp_control_values, tmp_control_v, M_ctr, true);
    }

    this.u_vv_splines = new Array(m); // bar M^{0, 2} (t, tau_l) = ... (u, v_j)
    for (let j = 0; j < m; ++j) {
      let tmp_control_values = new Array(N_ctr),
        tmp_control_u = new Array(N_ctr);
      for (let i = 0; i < n; ++i) {
        tmp_control_values[i] = v_splines[i].calc_second_diff(0, i);
        tmp_control_u[i] = control_u[i][j];
      }
      tmp_control_values[N_ctr - 1] = v_splines[N_ctr - 1].calc_second_diff(1, n - 1);
      tmp_control_u[N_ctr - 1] = control_u[N_ctr - 1][j];
      this.u_vv_splines[j] = new CubicSplineX(tmp_control_values, tmp_control_u, N_ctr);
    }

    this.d = new Array(m);
    for (let i = 0; i < m; ++i)
      this.d[i] = control_v[0][i + 1] - control_v[0][i];
  }

    calc_value(omega, xi, ii, jj) {
        let psi1 = 1 - xi,
            psi2 = xi,
            psi3 = xi * (xi - 1) * (2 - xi),
            psi4 = xi * (xi * xi - 1);
        let coefs = [this.u_splines[jj].calc_value(omega, ii),
            this.u_splines[jj].calc_value(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii)];
        return psi1 * coefs[0] +
            psi2 * coefs[1] +
            psi3 * this.d[ii] * this.d[ii] * coefs[2] +
            psi4 * this.d[ii] * this.d[ii] * coefs[3];
    }

    calc_tangent_u(omega, xi, ii, jj) {
        let psi1 = 1 - xi,
            psi2 = xi,
            psi3 = xi * (xi - 1) * (2 - xi),
            psi4 = xi * (xi * xi - 1);
        let coefs = [this.u_splines[jj].calc_value(omega, ii),
            this.u_splines[jj].calc_value(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii)];
        return coefs[0] * omega1 + coefs[1] * omega2 +
          coefs[2] * this.h[i] * this.h[i] * omega3 +
          coefs[3] * this.h[i] * this.h[i] * omega4;
    }
    calc_tangent_v(omega, xi, ii, jj) {
        let psi1 = - 1. / this.h[i],
            psi2 = 1. / this.h[i],
            psi3 = (6 * xi - 3 * xi * xi - 2) / this.h[i],
            psi4 = (3 * xi * xi - 1) / this.h[i];
        let coefs = [this.u_splines[jj].calc_value(omega, ii),
            this.u_splines[jj].calc_value(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii),
            this.u_vv_splines[jj].calc_second_diff(omega, ii)];
        return coefs[0] * omega1 + coefs[1] * omega2 +
          coefs[2] * this.h[i] * this.h[i] * omega3 +
          coefs[3] * this.h[i] * this.h[i] * omega4;
    }
}

class Spline {
  constructor(control_points, N_ctr, M_ctr) {
    let control_values_x = new Array(N_ctr),
      control_values_y = new Array(N_ctr),
      control_values_z = new Array(N_ctr),
      control_values_u = new Array(N_ctr),
      control_values_v = new Array(N_ctr);
    for (let i = 0; i < N_ctr; ++i) {
      control_values_x[i] = new Array(M_ctr);
      control_values_y[i] = new Array(M_ctr);
      control_values_z[i] = new Array(M_ctr);
      control_values_u[i] = new Array(M_ctr);
      control_values_v[i] = new Array(M_ctr);
      for (let j = 0; j < M_ctr; ++j) {
        control_values_x[i][j] = control_points[i][j].x;
        control_values_y[i][j] = control_points[i][j].y;
        control_values_z[i][j] = control_points[i][j].z;
        control_values_u[i][j] = control_points[i][j].u;
        control_values_v[i][j] = control_points[i][j].v;
      }
    }

    this.x_spline = new BicubicSplineX(control_values_x, control_values_u, control_values_v, N_ctr, M_ctr, 1);
    this.y_spline = new BicubicSplineX(control_values_y, control_values_u, control_values_v, N_ctr, M_ctr, 0);
    this.z_spline = new BicubicSplineX(control_values_z, control_values_u, control_values_v, N_ctr, M_ctr, 0);
  }

  calc_value(omega, xi, ii, jj) {
    return [this.x_spline.calc_value(omega, xi, ii, jj),
    this.y_spline.calc_value(omega, xi, ii, jj),
    this.z_spline.calc_value(omega, xi, ii, jj)]
  }

  calc_tangent_u(omega, xi, ii, jj) {
    return [this.x_spline.calc_tangent_u(omega, xi, ii, jj),
    this.y_spline.calc_tangent_u(omega, xi, ii, jj),
    this.z_spline.calc_tangent_u(omega, xi, ii, jj)]
  }

  calc_tangent_v(omega, xi, ii, jj) {
    return [this.x_spline.calc_tangent_v(omega, xi, ii, jj),
    this.y_spline.calc_tangent_v(omega, xi, ii, jj),
    this.z_spline.calc_tangent_v(omega, xi, ii, jj)]
  }
}
