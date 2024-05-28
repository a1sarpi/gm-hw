"use strict";

class CubicSplineX {
    constructor(control_values, control_t, N_ctr, is_periodic = false, granich = { num: '4' }) {
        console.assert(!isNaN(control_values[0]), 'Автор неуч, control_values[0] == NaN!  ' + control_values);
        //console.assert(!isNaN(control_values[N_ctr - 1]), 'Автор неуч, control_values[N_ctr - 1] == NaN!!  ' + control_values);
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
                matrix[i][(i - 1 + N_ctr) % (N_ctr)] = this.h[(i - 1 + N_ctr) % N_ctr];
                matrix[i][i] = 2 * (this.h[(i - 1 + N_ctr) % N_ctr] + this.h[(i) % (N_ctr)]);
                matrix[i][(i + 1) % (N_ctr)] = this.h[(i) % N_ctr];

                free[i] = 6 * this.delta[(i - 1 + N_ctr) % N_ctr];
            }

            let tmp = period_progon(matrix, free);
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
            
            if (granich.num == '1') {
                console.log(granich);
                for (let i = 0; i < N_ctr; ++i) {
                    if (i == 0) {
                        matrix[i][0] = 2*this.h[0];
                        matrix[i][1] = this.h[0];

                        free[i] = 6 * ( (control_values[1] - control_values[0]) / this.h[0] - granich.p0_prime );
                    } else if (i == N_ctr - 1) {
                        matrix[i][N_ctr - 2] = this.h[N_ctr - 2];
                        matrix[i][N_ctr - 1] = 2*this.h[N_ctr - 3];

                        free[i] = 6 * ( granich.p1_prime - (control_values[N_ctr-1] - control_values[N_ctr-2]) / this.h[N_ctr - 2] );
                    } else {
                        matrix[i][i - 1] = this.h[i - 1];
                        matrix[i][i] = 2 * (this.h[i - 1] + this.h[i]);
                        matrix[i][i + 1] = this.h[i];

                        free[i] = 6 * this.delta[i - 1];
                    }
                }
            } else {
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
            }

            let tmp = period_progon(matrix, free);
            for (let i = 0; i < N_ctr; ++i) {
                this.parameters_M[i] = tmp[i] / 6; // bar M
            }
        }

        console.assert(!isNaN(this.parameters_M[0]), 'Автор дважды неуч' + this.parameters_M);
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
            omega3 = (6 - 6 * omega),
            omega4 = (6 * omega);
        let coefs = [this.parameters_p[i], this.parameters_p[i+1], this.parameters_M[i], this.parameters_M[i+1] ];
        return coefs[0] * omega1 +
            coefs[1] * omega2 +
            coefs[2] * omega3 +
            coefs[3] * omega4;
    }
}


class BicubicSplineX {
    constructor(control_values, control_u, control_v, N_ctr, M_ctr, granich = { num: '4'}) {
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
            this.u_splines[j] = new CubicSplineX(tmp_control_values, tmp_control_u, N_ctr, false, granich);
        }


        console.log("Строим v_splines");
        this.v_splines = new Array(N_ctr);
        for (let i = 0; i < N_ctr; ++i) {
            let tmp_control_values = new Array(M_ctr),
                tmp_control_v = new Array(M_ctr);
            for (let j = 0; j < M_ctr; ++j) {
                tmp_control_values[j] = control_values[i][j];
                tmp_control_v[j] = control_v[i][j];
            }
            console.log('trash', i);
            this.v_splines[i] = new CubicSplineX(tmp_control_values, tmp_control_v, M_ctr, true);
        }

        console.log("Строим v_uu_splines");
        this.v_uu_splines = new Array(N_ctr); // bar M^{0, 2} (t, tau_l) = ... (u, v_j)
        for (let i = 0; i < N_ctr; ++i) {
            let tmp_control_values = new Array(M_ctr),
                tmp_control_v = new Array(M_ctr);
            if (i == N_ctr - 1) {
                for (let j = 0; j < M_ctr; ++j) {
                    tmp_control_values[j] = this.u_splines[j].calc_second_diff(1, N_ctr - 2);
                    tmp_control_v[j] = control_v[i][j];
                }
            } else {
                for (let j = 0; j < M_ctr; ++j) {
                    tmp_control_values[j] = this.u_splines[j].calc_second_diff(0, i);
                    tmp_control_v[j] = control_v[i][j];
                }
            }
            this.v_uu_splines[i] = new CubicSplineX(tmp_control_values, tmp_control_v, M_ctr, true);
        }

        this.d = new Array(N_ctr - 1);
        for (let i = 0; i < N_ctr - 1; ++i)
            this.d[i] = control_u[i + 1][0] - control_u[i][0];
        this.d[N_ctr - 1] = 1 - control_u[N_ctr - 1][0];
    }

    calc_value(omega, xi, ii, jj) {
        let psi1 = 1 - omega,
            psi2 = omega,
            psi3 = omega * (omega - 1) * (2 - omega),
            psi4 = omega * (omega * omega - 1);
        let coefs = [this.v_splines[ii].calc_value(xi, jj),
            this.v_splines[ii + 1].calc_value(xi, jj),
            this.v_uu_splines[ii].calc_value(xi, jj) / 6,
            this.v_uu_splines[ii + 1].calc_value(xi, jj) / 6 ];
        return psi1 * coefs[0] +
            psi2 * coefs[1] +
            psi3 * this.d[ii] * this.d[ii] * coefs[2] +
            psi4 * this.d[ii] * this.d[ii] * coefs[3];
    }

    calc_tangent_v(omega, xi, ii, jj) {
        let psi1 = 1 - omega,
            psi2 = omega,
            psi3 = omega * (omega - 1) * (2 - omega),
            psi4 = omega * (omega * omega - 1);
        let coefs = [this.v_splines[ii].calc_diff(xi, jj),
                this.v_splines[ii + 1].calc_diff(xi, jj),
                this.v_uu_splines[ii].calc_diff(xi, jj) / 6,
                this.v_uu_splines[ii + 1].calc_diff(xi, jj) / 6 ];
        let outval = psi1 * coefs[0] +
                psi2 * coefs[1] +
                psi3 * this.d[ii] * this.d[ii] * coefs[2] +
                psi4 * this.d[ii] * this.d[ii] * coefs[3];
        console.assert(!isNaN(outval), 'Автор неуч, outval == NaN! в calc_tangent_v ');
        return outval;
    }
    calc_tangent_u(omega, xi, ii, jj) {
        let psi1 = - 1. / this.d[ii],
            psi2 = 1. / this.d[ii],
            psi3 = (6 * omega - 3 * omega * omega - 2) / this.d[ii],
            psi4 = (3 * omega * omega - 1) / this.d[ii];
        let coefs = [this.v_splines[ii].calc_value(xi, jj),
            this.v_splines[ii + 1].calc_value(xi, jj),
            this.v_uu_splines[ii].calc_value(xi, jj) / 6,
            this.v_uu_splines[ii + 1].calc_value(xi, jj) / 6 ];
        console.assert(!isNaN(coefs[0]), 'Автор неуч, coefs[0] == NaN! в calc_tangent_u ');
        console.assert(!isNaN(coefs[1]), 'Автор неуч, coefs[1] == NaN! в calc_tangent_u ');
        console.assert(!isNaN(coefs[2]), 'Автор неуч, coefs[2] == NaN! в calc_tangent_u ');
        console.assert(!isNaN(coefs[3]), 'Автор неуч, coefs[3] == NaN! в calc_tangent_u ');
        let outval = psi1 * coefs[0] +
                psi2 * coefs[1] +
                psi3 * this.d[ii] * this.d[ii] * coefs[2] +
                psi4 * this.d[ii] * this.d[ii] * coefs[3];
        console.assert(!isNaN(outval), 'Автор неуч, outval == NaN! в calc_tangent_u ');
        return outval;
    }
}

class Spline {
    constructor(control_points, N_ctr, M_ctr, granich = '4') {
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
        if(granich == '4') {
            this.x_spline = new BicubicSplineX(control_values_x, control_values_u, control_values_v, N_ctr, M_ctr, { num: '4' });
            console.log("Строим сплайн y(u, v):");
            this.y_spline = new BicubicSplineX(control_values_y, control_values_u, control_values_v, N_ctr, M_ctr, { num: granich });
            console.log("Строим сплайн z(u, v):");
            this.z_spline = new BicubicSplineX(control_values_z, control_values_u, control_values_v, N_ctr, M_ctr, { num: granich });
        } else {
            this.x_spline = new BicubicSplineX(control_values_x, control_values_u, control_values_v, N_ctr, M_ctr, { num: granich,
                p0_prime: (control_values_x[N_ctr - 1][0] - control_values_x[0][0]),
                p1_prime: (control_values_x[N_ctr - 1][0] - control_values_x[0][0]) });
            console.log("Строим сплайн y(u, v):");
            this.y_spline = new BicubicSplineX(control_values_y, control_values_u, control_values_v, N_ctr, M_ctr, { num: granich, p0_prime: 0, p1_prime: 0 });
            console.log("Строим сплайн z(u, v):");
            this.z_spline = new BicubicSplineX(control_values_z, control_values_u, control_values_v, N_ctr, M_ctr, { num: granich, p0_prime: 0, p1_prime: 0 });
        }
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
