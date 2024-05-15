"use strict";

class BicubicSplineX {
    constructor(control_values, N_ctr, M_ctr, derivative) {
        const n = N_ctr - 1;
        const m = M_ctr; // Так как мы хотим замкнуть цилиндр вдоль полярной координаты

        this.parameters = new Array(n);
        for (let i = 0; i < n; i++) {
            this.parameters[i] = new Array(m);
            for (let j = 0; j < m; j++) {
                this.parameters[i][j] = new Array(4);
                for(let k = 0; k < 4; ++k) 
                    this.parameters[i][j][k] = new Array(4);
            }
        }

        for (let i = 0; i < n; ++i) {
            for (let j = 0; j < m; ++j) {
                for (let ii = 0; ii < 4; ++ii)
                    for (let jj = 0; jj < 4; ++jj)
                        this.parameters[i][j][ii][jj] = 0;
            }
        }


        // Посчитать кубические сплайны вдоль полярной координаты цилиндра
        let matrix = new Array(m*4),
            free = new Array(m*4);
        for (let i = 0; i < 4*m; ++i) {
            matrix[i] = new Array(m*4);
        }
        zeros(matrix, 4*m, 4*m);
        for (let j = 0; j < m; ++j) {
            // f_ij (0, 0) = x_{i,j}
            matrix[4*j][4*j] = 1.;

            if(j < m-1) {
                // f_ij (0, 1) = x_{i, j+1}
                matrix[4*j+1][4*j] = 1.;
                matrix[4*j+1][4*j+1] = 1.;
                matrix[4*j+1][4*j+2] = 1.;
                matrix[4*j+1][4*j+3] = 1.;
                matrix[4*j+1][4*(j+1)] = -1.;
                free[4*j + 1] = 0.;

                // f_ij_v (0, 1) = f_{i, j+1}_v (0, 0)
                matrix[4*j+2][4*j  ] = 0.;
                matrix[4*j+2][4*j+1] = 1.;
                matrix[4*j+2][4*j+2] = 2.;
                matrix[4*j+2][4*j+3] = 3.;
                matrix[4*j+2][4*(j+1)+1] = -1.;
                free[4*j+2] = 0.;

                // f_ij_vv (0, 1) = f_{i, j+1}_vv (0, 0)
                matrix[4*j+3][4*j  ] = 0.;
                matrix[4*j+3][4*j+1] = 0.;
                matrix[4*j+3][4*j+2] = 2.;
                matrix[4*j+3][4*j+3] = 6.;
                matrix[4*j+3][4*(j+1)+2] = -2.;
                free[4*j+3] = 0.;
            }
        }

        // Периодические граничные условия
        // f_im (0, 1) = f_i0 (0, 0)
        matrix[4*m - 3][4*(m-1)    ] = 1.;
        matrix[4*m - 3][4*(m-1) + 1] = 1.;
        matrix[4*m - 3][4*(m-1) + 2] = 1.;
        matrix[4*m - 3][4*(m-1) + 3] = 1.;
        matrix[4*m - 3][0] = -1.;
        free[4*m-3] = 0.;

        // f_im_v (0, 1) = f_i0_v (0, 0)
        matrix[4*m - 2][4*(m-1)    ] = 0.;
        matrix[4*m - 2][4*(m-1) + 1] = 1.;
        matrix[4*m - 2][4*(m-1) + 2] = 2.;
        matrix[4*m - 2][4*(m-1) + 3] = 3.;
        matrix[4*m - 2][1] = -1.;
        free[4*m - 2] = 0.;

        // f_im_vv (0, 1) = f_i0_vv (0, 0)
        matrix[4*m - 1][4*(m-1)    ] = 0.;
        matrix[4*m - 1][4*(m-1) + 1] = 0.;
        matrix[4*m - 1][4*(m-1) + 2] = 2.;
        matrix[4*m - 1][4*(m-1) + 3] = 6.;
        matrix[4*m - 1][2] = -2.;
        free[4*m - 1] = 0.;

        for (let i = 0; i < n; ++i) {    
            for (let j = 0; j < m; ++j)
                free[4*j] = control_values[i][j];

            let tmp = solve(matrix, free);
            for (let j = 0; j < m; ++j) {
                for (let ii = 0; ii < 4; ++ii)
                    this.parameters[i][j][0][ii] = tmp[4*j + ii];
            }
        }

        
        // /*
        // Посчитать кубические сплайны вдоль оси цилиндра
        matrix = new Array(4*n);
        free = new Array(4*n);
        for (let i = 0; i < 4*n; ++i) {
            matrix[i] = new Array(4*n);
        }
        zeros(matrix, 4*n, 4*n);
        for (let j = 0; j < n; ++j) {
            // f_ij (0, 0) = x_{i,j}
            matrix[4*j][4*j] = 1.;

            if(j < n-1) {
                // f_ij (0, 1) = x_{i, j+1}
                matrix[4*j+1][4*j] = 1.;
                matrix[4*j+1][4*j+1] = 1.;
                matrix[4*j+1][4*j+2] = 1.;
                matrix[4*j+1][4*j+3] = 1.;
                matrix[4*j+1][4*(j+1)] = -1.;
                free[4*j + 1] = 0.;

                // f_ij_u (0, 1) = f_{i, j+1}_u (0, 0)
                matrix[4*j+2][4*j  ] = 0.;
                matrix[4*j+2][4*j+1] = 1.;
                matrix[4*j+2][4*j+2] = 2.;
                matrix[4*j+2][4*j+3] = 3.;
                matrix[4*j+2][4*(j+1)+1] = -1.;
                free[4*j+2] = 0.;

                // f_ij_uu (0, 1) = f_{i, j+1}_uu (0, 0)
                matrix[4*j+3][4*j  ] = 0.;
                matrix[4*j+3][4*j+1] = 0.;
                matrix[4*j+3][4*j+2] = 2.;
                matrix[4*j+3][4*j+3] = 6.;
                matrix[4*j+3][4*(j+1)+2] = -2.;
                free[4*j+3] = 0.;
            }
        }

        // Граничные условия вдоль оси цилиндра
        matrix[4*n - 3][4*(n-1)    ] = 1.;
        matrix[4*n - 3][4*(n-1) + 1] = 1.;
        matrix[4*n - 3][4*(n-1) + 2] = 1.;
        matrix[4*n - 3][4*(n-1) + 3] = 1.;

        // f_u (u_max, v) = derivative
        matrix[4*n - 2][4*(n-1)    ] = 0.;
        matrix[4*n - 2][4*(n-1) + 1] = 1.;
        matrix[4*n - 2][4*(n-1) + 2] = 2.;
        matrix[4*n - 2][4*(n-1) + 3] = 3.;
        matrix[4*n - 2][1] = -1.;
        free[4*n - 2] = derivative;

        matrix[4*n - 1][0] = 0.;
        matrix[4*n - 1][1] = 0.;
        matrix[4*n - 1][2] = 2.;
        matrix[4*n - 1][3] = 6.;
        free[4*n - 1] = derivative;

        for (let j = 0; j < m; ++j) {    
            for (let i = 0; i < n; ++i)
                free[4*i] = control_values[i][j];
            free[4*n-3] = control_values[n][j];

            let tmp = solve(matrix, free);
            for (let i = 0; i < n; ++i) {
                for (let ii = 0; ii < 4; ++ii)
                    this.parameters[i][j][ii][0] = tmp[4*i + ii];
            }
        }

        // ==================================================================
        // до этой строчки точно правильно
    
        // склеивание f_v вдоль u
        matrix = new Array(n*3);
        free = new Array(n*3);
        for (let i = 0; i < 3*n; ++i) {
            matrix[i] = new Array(n*3);
        }
        zeros(matrix, 3*n, 3*n);
        for (let i = 0; i < n - 1; ++i) {
            // f_v_ij (1, 0)
            matrix[3*i][3*i] = 1.;
            matrix[3*i][3*i+1] = 1.;
            matrix[3*i][3*i+2] = 1.;

            // f_ij_vu (1, 0) = f_{i, j+1}_vu (0, 0)
            matrix[3*i+1][3*i  ] = 1.;
            matrix[3*i+1][3*i+1] = 2.;
            matrix[3*i+1][3*i+2] = 3.;
            matrix[3*i+1][3*(i+1)] = -1.;
            free[3*i+1] = 0.;

            // f_ij_vuu (1, 0) = f_{i, j+1}_vuu (0, 0)
            matrix[3*i+2][3*i  ] = 0.;
            matrix[3*i+2][3*i+1] = 2.;
            matrix[3*i+2][3*i+2] = 6.;
            matrix[3*i+2][3*(i+1)+1] = -2.;
            free[3*i+2] = 0.;
        }

        // Граничные условия вдоль оси цилиндра
        matrix[3*n - 3][3*(n-1)    ] = 1.;
        matrix[3*n - 3][3*(n-1) + 1] = 1.;
        matrix[3*n - 3][3*(n-1) + 2] = 1.;

        matrix[3*n - 2][4*(n-1)    ] = 1.;
        matrix[3*n - 2][4*(n-1) + 1] = 2.;
        matrix[3*n - 2][4*(n-1) + 2] = 3.;
        free[3*n - 2] = 0.;

        matrix[3*n - 1][0] = 0.;
        matrix[3*n - 1][1] = 2.;
        matrix[3*n - 1][2] = 6.;
        free[3*n - 1] = 0.;

        for (let j = 0; j < m; ++j) {
            for (let i = 0; i < n - 1; ++i)
                free[3*i] = this.parameters[i+1][j][0][1] - this.parameters[i][j][0][1];
            free[3*n-3] = - this.parameters[n-1][j][0][1];

            let tmp = solve(matrix, free);
            for (let i = 0; i < n; ++i) {
                for (let ii = 0; ii < 3; ++ii)
                    this.parameters[i][j][ii + 1][1] = tmp[3*i + ii];
            }
        }


        // склеивание f_u вдоль v
        matrix = new Array(m*3);
        free = new Array(m*3);
        for (let i = 0; i < 3*m; ++i) {
            matrix[i] = new Array(m*3);
        }
        zeros(matrix, 3*m, 3*m);
        for (let j = 0; j < m - 1; ++j) {
            // f_u_ij (0, 1)
            matrix[3*j][3*j] = 1.;
            matrix[3*j][3*j+1] = 1.;
            matrix[3*j][3*j+2] = 1.;

            // f_ij_vu (0, 1) = f_{i, j+1}_vu (0, 0)
            matrix[3*j+1][3*j  ] = 1.;
            matrix[3*j+1][3*j+1] = 2.;
            matrix[3*j+1][3*j+2] = 3.;
            matrix[3*j+1][3*(j+1)] = -1.;
            free[3*j+1] = 0.;

            // f_ij_vuu (0, 1) = f_{i, j+1}_vuu (0, 0)
            matrix[3*j+2][3*j  ] = 0.;
            matrix[3*j+2][3*j+1] = 2.;
            matrix[3*j+2][3*j+2] = 6.;
            matrix[3*j+2][3*(j+1)+1] = -2.;
            free[3*j+2] = 0.;
        }

        // Периодические граничные условия
        matrix[3*m - 3][3*(m-1)    ] = 1.;
        matrix[3*m - 3][3*(m-1) + 1] = 1.;
        matrix[3*m - 3][3*(m-1) + 2] = 1.;

        matrix[3*m - 2][3*(m-1)    ] = 1.;
        matrix[3*m - 2][3*(m-1) + 1] = 2.;
        matrix[3*m - 2][3*(m-1) + 2] = 3.;
        matrix[3*m - 2][0] = -1.;
        free[3*m - 2] = 0.;

        matrix[3*m - 1][3*(m-1)    ] = 0.;
        matrix[3*m - 1][3*(m-1) + 1] = 2.;
        matrix[3*m - 1][3*(m-1) + 2] = 6.;
        matrix[3*m - 1][1] = -2.;
        free[3*m - 1] = 0.;

        for (let i = 0; i < n; ++i) {    
            for (let j = 0; j < m - 1; ++j)
                free[3*j] = this.parameters[i][j+1][1][0] - this.parameters[i][j][1][0];
            free[3*m - 3] = this.parameters[i][0][1][0] - this.parameters[i][m-1][1][0];

            let tmp = solve(matrix, free);
            for (let j = 0; j < m; ++j) {
                for (let ii = 0; ii < 3; ++ii)
                    this.parameters[i][j][1][ii + 1] = tmp[3*j + ii];
            }
        }

        // ==================================================================
        // вроде до сюда тоже всё правильно, кроме возможно граничных условий


        // склеивание f_vv вдоль u
        matrix = new Array(n*2);
        free = new Array(n*2);
        for (let i = 0; i < 2*n; ++i) {
            matrix[i] = new Array(n*2);
        }
        zeros(matrix, 2*n, 2*n);
        for (let i = 0; i < n - 1; ++i) {
            // f_vv_ij (1, 0) = f_vv_{i+1, j} (0, 0)
            matrix[2*i][2*i] = 1.;
            matrix[2*i][2*i+1] = 1.;

            // f_ij_vvu (1, 0) = f_{i+1, j}_vvu (0, 0)
            matrix[2*i+1][2*i  ] = 2.;
            matrix[2*i+1][2*i+1] = 3.;
        }

        // Граничные условия вдоль оси цилиндра
        matrix[2*n - 2][2*(n-1)    ] = 1.;
        matrix[2*n - 2][2*(n-1) + 1] = 1.;

        matrix[2*n - 1][2*(n-1)    ] = 2.;
        matrix[2*n - 1][2*(n-1) + 1] = 3.;

        for (let j = 0; j < m; ++j) {
            for (let i = 0; i < n - 1; ++i) {
                free[2*i] = this.parameters[i+1][j][0][2] - this.parameters[i][j][0][2] - this.parameters[i][j][1][2];
                free[2*i + 1] = this.parameters[i+1][j][1][2] - this.parameters[i][j][1][2];
            }
            free[2*n-2] =  - this.parameters[n-1][j][0][2] - this.parameters[n-1][j][1][2];
            free[2*n - 1] = - this.parameters[n-1][j][1][2];

            let tmp = solve(matrix, free);
            for (let i = 0; i < n; ++i) {
                for (let ii = 0; ii < 2; ++ii)
                    this.parameters[i][j][ii + 2][2] = tmp[2*i + ii];
            }
        }

        // до сюда вроде нормально


        // склеивание f_uu вдоль v
        matrix = new Array(m*2);
        free = new Array(m*2);
        for (let i = 0; i < 2*m; ++i) {
            matrix[i] = new Array(m*2);
        }
        zeros(matrix, 2*m, 2*m);
        for (let j = 0; j < m - 1; ++j) {
            // f_uu_ij (0, 1) = f_uu_{i, j+1} (0, 0)
            matrix[2*j][2*j] = 1.;
            matrix[2*j][2*j+1] = 1.;

            // f_ij_uuv (0, 1) = f_{i, j+1}_uuv (0, 0)
            matrix[2*j+1][2*j  ] = 2.;
            matrix[2*j+1][2*j+1] = 3.;
        }

        // Периодические граничные условия
        matrix[2*m - 2][2*(m-1)    ] = 1.;
        matrix[2*m - 2][2*(m-1) + 1] = 1.;

        matrix[2*m - 1][2*(m-1)    ] = 2.;
        matrix[2*m - 1][2*(m-1) + 1] = 6.;

        for (let i = 0; i < n; ++i) {    
            for (let j = 0; j < m - 1; ++j) {
                free[2*j] = this.parameters[i][j+1][2][0] - this.parameters[i][j][2][0] - this.parameters[i][j][2][1];
                free[2*j + 1] = this.parameters[i][j+1][2][1] - this.parameters[i][j][2][1];
            }
            free[2*m-2] = this.parameters[i][0][2][0] - this.parameters[i][m-1][2][0] - this.parameters[i][m-1][2][1];
            free[2*m-1] = this.parameters[i][0][2][1] - this.parameters[i][m-1][2][1];

            let tmp = solve(matrix, free);
            for (let j = 0; j < m; ++j) {
                for (let ii = 0; ii < 2; ++ii)
                    this.parameters[i][j][2][ii + 2] = tmp[2*j + ii];
            }
        }

        

        for (let i = 0; i < n; ++i) {
            for (let j = 0; j < m - 1; ++j)
                this.parameters[i][j][3][3] = this.parameters[i][j+1][3][0] -
                    this.parameters[i][j][3][0] - 
                    this.parameters[i][j][3][1] -
                    this.parameters[i][j][3][2];
            this.parameters[i][m-1][3][3] = this.parameters[i][0][3][0] -
                    this.parameters[i][m-1][3][0] - 
                    this.parameters[i][m-1][3][1] -
                    this.parameters[i][m-1][3][2];
        } // */
    }

    calc_value(omega, xi, ii, jj) {
        console.log(xi);
        let out_value = 0;
        for (let i = 0; i < 4; ++i) 
            for (let j = 0; j < 4; ++j)
                out_value += this.parameters[ii][jj][i][j] * Math.pow(omega, i) * Math.pow(xi, j);
        return out_value;
    }

    calc_tangent_u(omega, xi, ii, jj) {
        let out_value = 0;
        for (let i = 1; i < 4; ++i) 
            for (let j = 0; j < 4; ++j)
                out_value += i * this.parameters[ii][jj][i][j] * Math.pow(omega, i-1) * Math.pow(xi, j);
        return out_value;
    }
    calc_tangent_v(omega, xi, ii, jj) {
        let out_value = 0;
        for (let i = 0; i < 4; ++i) 
            for (let j = 1; j < 4; ++j)
                out_value += j * this.parameters[ii][jj][i][j] * Math.pow(omega, i) * Math.pow(xi, j-1);
        return out_value;
    }
}

class Spline {
    constructor(control_points, N_ctr, M_ctr) {
        let control_values_x = new Array(N_ctr),
            control_values_y = new Array(N_ctr),
            control_values_z = new Array(N_ctr);
            //control_values_u = new Array(N_ctr),
            //control_values_v = new Array(N_ctr);
        for (let i = 0; i < N_ctr; ++i) {
            control_values_x[i] = new Array(M_ctr);
            control_values_y[i] = new Array(M_ctr);
            control_values_z[i] = new Array(M_ctr);
            //control_values_u[i] = new Array(M_ctr);
            //control_values_v[i] = new Array(M_ctr);
            for (let j = 0; j < M_ctr; ++j) {
                control_values_x[i][j] = control_points[i][j].x;
                control_values_y[i][j] = control_points[i][j].y;
                control_values_z[i][j] = control_points[i][j].z;
                //control_values_u[i][j] = control_points[i][j].u;
                //control_values_v[i][j] = control_points[i][j].v;
            }
        }

        this.x_spline = new BicubicSplineX(control_values_x, N_ctr, M_ctr, 1);
        this.y_spline = new BicubicSplineX(control_values_y, N_ctr, M_ctr, 0);
        this.z_spline = new BicubicSplineX(control_values_z, N_ctr, M_ctr, 0);
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
