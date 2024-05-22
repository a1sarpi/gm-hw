"use strict";

class CubicSplineX {
    constructor(control_values, control_t, N_ctr, is_periodic = false) {
        console.assert(!isNaN(control_values[0]), 'Автор долбаёб!' + control_values);
        if (is_periodic) {
            control_values.push(control_values[0]);
            this.h = new Array(N_ctr); // длины отрезков между к.т. (по параметру)
            for (let i = 0; i < N_ctr; ++i) {
                if (i + 1 == N_ctr)
                    this.h[i] = 1 - control_t[i];
                else
                    this.h[i] = control_t[i + 1] - control_t[i];
            }

            this.delta = new Array(N_ctr);
            for (let i = 0; i < N_ctr; ++i) {
                this.delta[i] = ( control_values[ (i + 2) % (N_ctr)] - control_values[ (i + 1) % (N_ctr)] ) / this.h[ (i + 1) % N_ctr] - 
                    ( control_values[ (i + 1) % (N_ctr) ] - control_values[i] ) / this.h[i];
            }
        } else {
            this.h = new Array(N_ctr - 1); // длины отрезков между к.т. (по параметру)
            for (let i = 0; i < N_ctr - 1; ++i) {
                this.h[i] = control_t[i + 1] - control_t[i];
            }

            this.delta = new Array(N_ctr - 2);
            for (let i = 0; i < N_ctr - 2; ++i) {
                this.delta[i] = ( control_values[i + 2] - control_values[i + 1] ) / this.h[i + 1] - 
                    ( control_values[i + 1] - control_values[i] ) / this.h[i];
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
                /* if (i == 0) {
                    matrix[i][0] = 2 * (this.h[0] + this.h[1]);
                    matrix[i][1] = this.h[1];
                    matrix[i][N_ctr - 1] = this.h[0];

                    free[i] = 6 * this.delta[0];
                } else if (i == N_ctr - 1) {
                    matrix[i][N_ctr - 1] = 2 * (this.h[0] + this.h[N_ctr - 1]);
                    matrix[i][N_ctr - 2] = this.h[N_ctr - 1];
                    matrix[i][0] = this.h[0];

                    free[i] = 6 * ( (control_values[1] - control_values[0]) / this.h[0] - (control_values[0] - control_values[N_ctr - 1]) / this.h[N_ctr - 1] );
                } else {
                    matrix[i][i] = 2 * (this.h[i] + this.h[i + 1]);
                    matrix[i][i - 1] = this.h[i];
                    matrix[i][i + 1] = this.h[i + 1];

                    free[i] = 6 * this.delta[i];
                } */

                matrix[i][(i - 1 + N_ctr) % (N_ctr)] = this.h[(i - 1 + N_ctr) % N_ctr];
                matrix[i][i] = 2 * (this.h[(i - 1 + N_ctr) % N_ctr] + this.h[(i) % (N_ctr)]);
                matrix[i][(i + 1) % (N_ctr)] = this.h[(i) % N_ctr];

                free[i] = 6 * this.delta[(i - 1 + N_ctr) % N_ctr];
            }

            console.log("control_values = ", control_values);
            console.log("h = ", this.h);
            console.log("delta = ", this.delta);
            console.log("Ax = b:", matrix, free);

            let tmp = solve(matrix, free);
            for (let i = 0; i < N_ctr; ++i) {
                this.parameters_M[i] = tmp[i] / 6; // bar M
            }

            this.parameters_M[N_ctr] = this.parameters_M[0];
            //this.parameters_M[0] = this.parameters_M[N_ctr];
            console.log("M_i = ", this.parameters_M);
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
                    matrix[i][N_ctr - 3] = this.h[N_ctr - 2];
                    matrix[i][N_ctr - 2] = - ( this.h[N_ctr - 3] + this.h[N_ctr - 2] );
                    matrix[i][N_ctr - 1] = this.h[N_ctr - 3];

                    free[i] = 0;
                } else {
                    matrix[i][i - 1] = this.h[i - 1];
                    matrix[i][i] = 2 * (this.h[i - 1] + this.h[i]);
                    matrix[i][i + 1] = this.h[i];

                    free[i] = 6 * this.delta[i - 1];
                }
            }

            //console.log(this.h);
            //console.log(matrix, free);

            let tmp = solve(matrix, free);
            for (let i = 0; i < N_ctr; ++i) {
                this.parameters_M[i] = tmp[i] / 6; // bar M
            }
            //console.log(this.parameters_M);
        }

        console.assert(!isNaN(this.parameters_M[0]), 'Автор дважды долбаёб' + this.parameters_M);
    }

    calc_value(omega, i) {
        let omega1 = 1 - omega, // psi
            omega2 = omega,
            omega3 = omega * (omega - 1) * (2 - omega),
            omega4 = omega * (omega * omega - 1);
        let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
        //let coefs = [this.parameters_p[i], this.parameters_p[i+1], 0, 0 ];
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
        //let coefs = [this.parameters_p[i], this.parameters_p[i+1], 0, 0 ];
        return coefs[0] * omega1 + coefs[1] * omega2 +
            coefs[2] * this.h[i] * this.h[i] * omega3 +
            coefs[3] * this.h[i] * this.h[i] * omega4;
    }

    calc_second_diff(omega, i) {
        let omega1 = 0.,
            omega2 = 0.,
            omega3 = (6 - 6 * omega),
            omega4 = (6 * omega);
        let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
        //let coefs = [this.parameters_p[i], this.parameters_p[i+1], 0, 0 ];
        return coefs[0] * omega1 +
            coefs[1] * omega2 +
            coefs[2] * omega3 +
            coefs[3] * omega4;
    }
}


class BicubicSplineX {
    constructor(control_values, control_u, control_v, N_ctr, M_ctr, derivative) {
        this.N_ctr = N_ctr;
        this.M_ctr = M_ctr;

        const n = N_ctr - 1;
        const m = M_ctr; // Так как мы хотим замкнуть цилиндр вдоль полярной координаты

        console.log("Строим u_splines");
        this.u_splines = new Array(M_ctr); // S(t, tau_l) = S(u, v_j)
        for (let j = 0; j < M_ctr; ++j) {
            let tmp_control_values = new Array(N_ctr),
                tmp_control_u = new Array(N_ctr);
            for (let i = 0; i < N_ctr; ++i) {
                tmp_control_values[i] = control_values[i][j];
                tmp_control_u[i] = control_u[i][j];
            }
            this.u_splines[j] = new CubicSplineX(tmp_control_values, tmp_control_u, N_ctr);
        }


        console.log("Строим v_splines");
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

        console.log("Строим u_vv_splines");
        this.u_vv_splines = new Array(M_ctr); // bar M^{0, 2} (t, tau_l) = ... (u, v_j)
        for (let j = 0; j < M_ctr; ++j) {
            let tmp_control_values = new Array(N_ctr),
                tmp_control_u = new Array(N_ctr);
            for (let i = 0; i < N_ctr - 1; ++i) {
                tmp_control_values[i] = v_splines[i].calc_second_diff(0, j);
                tmp_control_u[i] = control_u[i][j];
            }
            tmp_control_values[N_ctr - 1] = v_splines[N_ctr - 1].calc_second_diff(1, M_ctr - 1);
            tmp_control_u[N_ctr - 1] = control_u[N_ctr - 1][j];
            //console.log(tmp_control_values);
            this.u_vv_splines[j] = new CubicSplineX(tmp_control_values, tmp_control_u, N_ctr);
        }

        this.d = new Array(M_ctr);
        for (let i = 0; i < m; ++i)
            this.d[i] = control_v[0][i + 1] - control_v[0][i];
        this.d[M_ctr - 1] = 1 - control_v[0][m - 1];
    }

    get_coefs(omega, xi, ii, jj) {
        if (jj == this.M_ctr - 1) {
            return [this.u_splines[jj].calc_value(omega, ii),
                this.u_splines[0].calc_value(omega, ii),
                this.u_vv_splines[jj].calc_value(omega, ii) / 6,
                this.u_vv_splines[0].calc_value(omega, ii) / 6 ];
        } else {
            return [this.u_splines[jj].calc_value(omega, ii),
                this.u_splines[jj + 1].calc_value(omega, ii),
                this.u_vv_splines[jj].calc_value(omega, ii) / 6,
                this.u_vv_splines[jj + 1].calc_value(omega, ii) / 6 ];
                //0, 0];
        }

    } 

    calc_value(omega, xi, ii, jj) {
        let psi1 = 1 - xi,
            psi2 = xi,
            psi3 = xi * (xi - 1) * (2 - xi),
            psi4 = xi * (xi * xi - 1);
        let coefs = this.get_coefs(omega, xi, ii, jj);
        return psi1 * coefs[0] +
            psi2 * coefs[1] +
            psi3 * this.d[jj] * this.d[jj] * coefs[2] +
            psi4 * this.d[jj] * this.d[jj] * coefs[3];
    }

    calc_tangent_u(omega, xi, ii, jj) {
        let psi1 = 1 - xi,
            psi2 = xi,
            psi3 = xi * (xi - 1) * (2 - xi),
            psi4 = xi * (xi * xi - 1);
        let coefs = this.get_coefs(omega, xi, ii, jj);
        if (jj == this.M_ctr - 1) {
            coefs = [this.u_splines[jj].calc_diff(omega, ii),
                this.u_splines[0].calc_diff(omega, ii),
                this.u_vv_splines[jj].calc_diff(omega, ii) / 6,
                this.u_vv_splines[0].calc_diff(omega, ii) / 6 ];
        } else {
            coefs = [this.u_splines[jj].calc_diff(omega, ii),
                this.u_splines[jj + 1].calc_diff(omega, ii),
                this.u_vv_splines[jj].calc_diff(omega, ii) / 6,
                this.u_vv_splines[jj + 1].calc_diff(omega, ii) / 6 ];
                //0, 0];
        }
        return psi1 * coefs[0] / this.d[jj] +
                psi2 * coefs[1] / this.d[jj]+
                psi3 * this.d[jj] * coefs[2] +
                psi4 * this.d[jj] * coefs[3];
    }
    calc_tangent_v(omega, xi, ii, jj) {
        let psi1 = - 1. / this.d[ii],
            psi2 = 1. / this.d[ii],
            psi3 = (6 * xi - 3 * xi * xi - 2) / this.d[ii],
            psi4 = (3 * xi * xi - 1) / this.d[ii];
        let coefs = this.get_coefs(omega, xi, ii, jj);
        return psi1 * coefs[0] +
                psi2 * coefs[1] +
                psi3 * this.d[jj] * this.d[jj] * coefs[2] +
                psi4 * this.d[jj] * this.d[jj] * coefs[3];
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
        console.log("Строим сплайн x(u, v):");
        this.x_spline = new BicubicSplineX(control_values_x, control_values_u, control_values_v, N_ctr, M_ctr, 1);
        console.log("Строим сплайн y(u, v):");
        this.y_spline = new BicubicSplineX(control_values_y, control_values_u, control_values_v, N_ctr, M_ctr, 0);
        console.log("Строим сплайн z(u, v):");
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
