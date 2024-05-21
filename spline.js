"use strict";

class CubicSplineX {
  constructor(control_values, control_t, N_ctr, is_periodic = false) {
    let n;
    if (!is_periodic)
      n = N_ctr - 1;
    else
      n = N_ctr;

    let h = new Array(n); // длины отрезков между к.т. (по параметру)
    for (let i = 0; i < n; ++i) {
      if (i + 1 == N_ctr && is_periodic)
        h[i] = 1 - control_t[i];
      else
        h[i] = control_t[i + 1] - control_t[i];
    }
    this.h = h;

    this.parameters = new Array(n);
    for (let i = 0; i < n; ++i) {
      this.parameters[i] = new Array(4); // четыре параметра для каждого отрезка
    }

    for (let i = 0; i < n; ++i) {
      this.parameters[i][0] = control_values[i];
      if (i + 1 == N_ctr && is_periodic)
        this.parameters[i][1] = control_values[0];
      else
        this.parameters[i][1] = control_values[i + 1];
    }

    let matrix = new Array(2 * n),
      free = new Array(2 * n);
    for (let i = 0; i < 2 * n; ++i)
      matrix[i] = new Array(2 * n);
    zeros(matrix, 2 * n, 2 * n);
    for (let i = 0; i < n - 1; ++i) {
      matrix[2 * i][2 * i] = - 2.;
      matrix[2 * i][2 * i + 1] = 6.;
      matrix[2 * i][2 * (i + 1)] = - 4.;
      matrix[2 * i][2 * (i + 1) + 1] = - 2.;

      matrix[2 * i + 1][2 * i] = - 1. * h[i];
      matrix[2 * i + 1][2 * i + 1] = 2. * h[i];
      matrix[2 * i + 1][2 * (i + 1)] = 2. * h[i + 1];
      matrix[2 * i + 1][2 * (i + 1) + 1] = 1. * h[i + 1];

      free[2 * i] = 0;
      free[2 * i + 1] = (control_values[i + 2] - control_values[i + 1]) / h[i + 1]
        - (control_values[i + 1] - control_values[i]) / h[i];
    }
    if (is_periodic) {
      matrix[2 * n - 2][2 * (n - 1)] = - 2.;
      matrix[2 * n - 2][2 * (n - 1) + 1] = 6.;
      matrix[2 * n - 2][0] = - 4.;
      matrix[2 * n - 2][1] = - 2.;

      matrix[2 * n - 1][2 * (n - 1)] = - 1. * h[n - 1];
      matrix[2 * n - 1][2 * (n - 1) + 1] = 2. * h[n - 1];
      matrix[2 * n - 1][0] = 2. * h[n - 1];
      matrix[2 * n - 1][1] = 1. * h[n - 1];

      free[2 * n - 2] = 0.;
      free[2 * n - 1] = (control_values[0] - control_values[n - 1]) / h[n - 1]
        - (control_values[n - 1] - control_values[n - 2]) / h[n - 1];
    } else {
      // Граничные условия первого типа:

    }

    let tmp = solve(matrix, free);
    for (let i = 0; i < n; ++i) {
      for (let ii = 0; ii < 2; ++ii)
        this.parameters[i][2 + ii] = tmp[2 * i + ii];
    }
  }

  calc_value(omega, i) {
    //let out_value = 0;
    //for (let ii = 0; ii < 4; ++ii)
    //    out_value += this.parameters[i][ii] * Math.pow(omega, ii);
    //return out_value;
    let omega1 = 1 - omega, // psi
      omega2 = omega,
      omega3 = omega * (omega - 1) * (2 - omega),
      omega4 = omega * (omega * omega - 1);
    let coefs = this.parameters[i];
    return coefs[0] * omega1 + coefs[1] * omega2 +
      coefs[2] * this.h[i] * this.h[i] * omega3 +
      coefs[3] * this.h[i] * this.h[i] * omega4;
  }

  calc_diff(omega, i) {
    let omega1 = - 1. / this.h[i],
      omega2 = 1. / this.h[i],
      omega3 = (6 * omega - 3 * omega * omega - 2) / this.h[i],
      omega4 = (3 * omega * omega - 1) / this.h[i];
    let coefs = this.parameters[i];
    return coefs[0] * omega1 + coefs[1] * omega2 +
      coefs[2] * this.h[i] * this.h[i] * omega3 +
      coefs[3] * this.h[i] * this.h[i] * omega4;
  }

  calc_second_diff(omega, i) {
    let omega1 = 0.,
      omega2 = 0.,
      omega3 = (6 - 6 * omega) / (this.h[i] * this.h[i]),
      omega4 = (6 * omega) / (this.h[i] * this.h[i]);
    let coefs = this.parameters[i];
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

    this.u_splines = new Array(m);
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

    this.u_vv_splines = new Array(m);
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
    let psi1 = 1 - omega,
      psi2 = omega,
      psi3 = omega * (omega - 1) * (2 - omega),
      psi4 = omega * (omega * omega - 1);
    let coefs = [this.u_splines[jj].calc_value(0, ii),
    this.u_splines[jj].calc_value(1, ii),
    this.u_vv_splines[jj].calc_second_diff(xi, ii),
    this.u_vv_splines[jj].calc_second_diff(xi, ii)];
    return psi1 * coefs[0] +
      psi2 * coefs[1] +
      psi3 * this.d[ii] * this.d[ii] * coefs[2] +
      psi4 * this.d[ii] * this.d[ii] * coefs[3];
  }

  calc_tangent_u(omega, xi, ii, jj) {
    //let psi1 = 1 - omega,
    //    psi2 = omega,
    //    psi3 = omega * (omega - 1) * (2 - omega),
    //    psi4 = omega * (omega*omega - 1);
    //let coefs = [ this.u_splines[jj].calc_value(0, ii),
    //              this.u_splines[jj].calc_value(1, ii),
    //              this.u_splines[jj].calc_second_diff(xi, ii),
    //              this.u_splines[jj].calc_second_diff(xi, ii) ];
    //return psi1 * coefs[0] +
    //    psi2 * coefs[1] +
    //    psi3 * this.d[ii]*this.d[ii] * coefs[2] +
    //    psi4 * this.d[ii]*this.d[ii] * coefs[3];
  }
  calc_tangent_v(omega, xi, ii, jj) {
    //let out_value = 0;
    //for (let i = 0; i < 4; ++i) 
    //    for (let j = 1; j < 4; ++j)
    //        out_value += j * this.parameters[ii][jj][i][j] * Math.pow(omega, i) * Math.pow(xi, j-1);
    //return out_value;
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
