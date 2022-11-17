#ifndef temp_h_INCLUDE
#define temp_h_INCLUDE

#include "models/heisenberg.hpp"
#include "random/mersenne.hpp"

heis_spin SimpleEnergySumm(const heis_spin *const *const *const simple, 
                        const heis_spin *const *const *const xface, 
                        const heis_spin *const *const *const yface, 
                        const heis_spin *const *const *const zface,  
                        int i, int j, int k, int L, int N,
                        double J2, double J_second) {
    heis_spin res;
    res.x = res.y = res.z = 0.0;

    int iNext, jNext, kNext;
    int iPrev, jPrev, kPrev;
    double JInterlayerUp, JInterlayerDown;

    /*-------------- Neighbours for xFace ----------------*/
    jNext = (j == (L - 2) ? 0 : j);
    kNext = (k == (N - 2) ? -1 : k);
    jPrev = (j == (0) ? L - 2 : j - 1);
    kPrev = (k == (0) ? -1 : k - 1);
    JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
    JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

    if (kNext != -1) {
        res += JInterlayerUp * xface[i][jNext][kNext];
        res += JInterlayerUp * xface[i][jPrev][kNext];
    }

    if (kPrev != -1) {
        res += JInterlayerDown * xface[i][jPrev][kPrev];
        res += JInterlayerDown * xface[i][jNext][kPrev];
    }

    /*-------------- Neighbors for yFace ----------------*/
    iNext = (i == (L - 2) ? 0 : i);
    kNext = (k == (N - 2) ? -1 : k);
    iPrev = (i == (0) ? L - 2 : i - 1);
    kPrev = (k == (0) ? -1 : k - 1);

    if (kNext != -1) {
        res += JInterlayerUp * yface[iNext][j][kNext];  
        res += JInterlayerUp * yface[iPrev][j][kNext];
    }

    if (kPrev != -1) {
        res += JInterlayerDown * yface[iPrev][j][kPrev];
        res += JInterlayerDown * yface[iNext][j][kPrev];
    }

    /*-------------- Neighbors for zFace ----------------*/
    iNext = (i == (L - 2) ? 0 : i);
    jNext = (j == (L - 2) ? 0 : j);
    iPrev = (i == (0) ? L - 2 : i - 1);
    jPrev = (j == (0) ? L - 2 : j - 1);

    res += zface[iNext][jNext][k];
    res += zface[iPrev][jNext][k];
    res += zface[iPrev][jPrev][k];
    res += zface[iNext][jPrev][k];

    /*------------ Neighbours of second order ---------------*/
    if (i == 0) {
        res += J_second * simple[L - 1][j][k];
    } else {
        res += J_second * simple[i - 1][j][k];
    }

    if (i == L - 1) {
        res += J_second * simple[0][j][k];
    } else {
        res += J_second * simple[i + 1][j][k];
    }

    if (j == 0) {
        res += J_second * simple[i][L - 1][k];      
    } else {
        res += J_second * simple[i][j - 1][k];
    }

    if (j == L - 1) {
        res += J_second * simple[i][0][k];
    } else {
        res += J_second * simple[i][j + 1][k];
    }

    if (k != 0) {
        if (k == N / 2) {
            res += J2 * simple[i][j][k - 1];
        } else {
            res += J_second * simple[i][j][k - 1];
        }
    }

    if (k != N - 2) {
        if (k == N / 2 - 1) {
            res += J2 * simple[i][j][k + 1];
        } else {
            res += J_second * simple[i][j][k + 1];
        }
    }

    return res;
}

heis_spin FaceEnergySumm(const heis_spin *const *const *const simple, 
                      const heis_spin *const *const *const xface, 
                      const heis_spin *const *const *const yface, 
                      const heis_spin *const *const *const zface, 
                      int i, int j, int k, int L, int N,
                      double J2, double J_second, char type) {
    heis_spin res;
    res.x = res.y = res.z = 0.0;

    int iNext, jNext, kNext;
    int iPrev, jPrev, kPrev;
    double JInterlayerUp, JInterlayerDown;

    switch (type) {
        case 'x':
            // ----------------------- first order ----------------------------------//
            if (i != L - 2) {
                res += yface[i][j][k];
                res += yface[i][j + 1][k];
                res += zface[i][j][k];
                res += zface[i][j][k + 1];
            } else {
                res += yface[0][j][k];
                res += yface[0][j + 1][k];
                res += zface[0][j][k];
                res += zface[0][j][k + 1];
            }

            if (i != 0) {
                res += yface[i - 1][j][k];
                res += yface[i - 1][j + 1][k];
                res += zface[i - 1][j][k];
                res += zface[i - 1][j][k + 1];
            } else {
                res += yface[L - 2][j][k];
                res += yface[L - 2][j + 1][k];
                res += zface[L - 2][j][k];
                res += zface[L - 2][j][k + 1];
            }

            res += simple[i][j][k];
            res += simple[i][j + 1][k];
            res += simple[i][j][k + 1];
            res += simple[i][j + 1][k + 1];

            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * xface[iNext][j][k];
            res += J_second * xface[iPrev][j][k];
            res += J_second * xface[i][jNext][k];
            res += J_second * xface[i][jPrev][k];

            if (kNext != -1) {
                res += J_second * JInterlayerUp * xface[i][j][kNext];
            }

            if (kPrev != -1) {
                res += J_second * JInterlayerDown * xface[i][j][kPrev];
            }

            break;


        case 'y':
            // ----------------------- first order ----------------------------------//
            if (j != L - 2) {
                res += xface[i][j][k];
                res += xface[i + 1][j][k];
                res += zface[i][j][k];
                res += zface[i][j][k + 1];
            } else {
                res += xface[i][0][k];
                res += xface[i + 1][0][k];
                res += zface[i][0][k];
                res += zface[i][0][k + 1];
            }

            if (j != 0) {
                res += xface[i][j - 1][k];
                res += xface[i + 1][j - 1][k];
                res += zface[i][j - 1][k];
                res += zface[i][j - 1][k + 1];
            } else {
                res += xface[i][L - 2][k];
                res += xface[i + 1][L - 2][k];
                res += zface[i][L - 2][k];
                res += zface[i][L - 2][k + 1];
            }

            res += simple[i][j][k];
            res += simple[i + 1][j][k];
            res += simple[i][j][k + 1];
            res += simple[i][j + 1][k + 1];

            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * yface[iNext][j][k];
            res += J_second * yface[iPrev][j][k];
            res += J_second * yface[i][jNext][k];
            res += J_second * yface[i][jPrev][k];

            if (kNext != -1) {
                res += J_second * JInterlayerUp * yface[i][j][kNext];
            }

            if (kPrev != -1) {
                res += J_second * JInterlayerDown * yface[i][j][kPrev];
            }

            break;


        case 'z':
            // ----------------------- first order ----------------------------------//

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            if (k != N - 2) {
                res += JInterlayerUp * xface[i][j][k];
                res += JInterlayerUp * xface[i + 1][j][k];
                res += JInterlayerUp * yface[i][j][k];
                res += JInterlayerUp * yface[i][j + 1][k];
            }

            if (k != 0) {
                res += JInterlayerDown * xface[i][j][k - 1];
                res += JInterlayerDown * xface[i + 1][j][k - 1];
                res += JInterlayerDown * yface[i][j][k - 1];
                res += JInterlayerDown * yface[i][j + 1][k - 1];
            }

            res += simple[i][j][k];
            res += simple[i + 1][j][k];
            res += simple[i][j + 1][k];
            res += simple[i + 1][j + 1][k];

            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * zface[iNext][j][k];
            res += J_second * zface[iPrev][j][k];
            res += J_second * zface[i][jNext][k];
            res += J_second * zface[i][jPrev][k];

            if (kNext != -1) {
                res += J_second * JInterlayerUp * zface[i][j][kNext];
            }

            if (kPrev != -1) {
                res += J_second * JInterlayerDown * zface[i][j][kPrev];
            }

            break;

        default:
            break;
    }

    return res;
}

heis_spin Magnetic1(const heis_spin *const *const *const simple, 
                 const heis_spin *const *const *const xface, 
                 const heis_spin *const *const *const yface, 
                 const heis_spin *const *const *const zface, 
                 int L, int N) {
    heis_spin res;
    res.x = res.y = res.z = 0.0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            for (int k = 0; k < N / 2; ++k) {
                res += simple[i][j][k];
                res += xface[i][j][k];
                res += yface[i][j][k];
                res += zface[i][j][k];
            }
        }
    }

    res.x /= 2 * L * L * N;
    res.y /= 2 * L * L * N;
    res.z /= 2 * L * L * N;

    return res;
}

heis_spin Magnetic2(const heis_spin *const *const *const simple, 
                 const heis_spin *const *const *const xface, 
                 const heis_spin *const *const *const yface, 
                 const heis_spin *const *const *const zface, 
                 int L, int N) {
    heis_spin res;
    res.x = res.y = res.z = 0.0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            for (int k = N / 2; k < N; ++k) {
                res += simple[i][j][k];
                res += xface[i][j][k];
                res += yface[i][j][k];
                res += zface[i][j][k];
            }
        }
    }
    
    res.x /= 2 * L * L * N;
    res.y /= 2 * L * L * N;
    res.z /= 2 * L * L * N;

    return res;
}

double S(double x) {
    if (x < 0.0) {
        return -1.0;
    }
    return 1.0;
}

double GMRSimpleEnergySumm(const heis_spin *const *const *const simple, double ***simple_up_0, double ***simple_down_0, 
                           const heis_spin *const *const *const xface, double ***xface_up_0, double ***xface_down_0, 
                           const heis_spin *const *const *const yface, double ***yface_up_0, double ***yface_down_0,  
                           const heis_spin *const *const *const zface, double ***zface_up_0, double ***zface_down_0, 
                           int i, int j, int k, int L, int N, double J2, double J_second, int move) {
    double res = 0.0;

    int iNext, jNext, kNext;
    int iPrev, jPrev, kPrev;
    double JInterlayerUp, JInterlayerDown;

    /*-------------- Neighbors for xFace ----------------*/
    jNext = (j == (L - 2) ? 0 : j);
    kNext = (k == (N - 2) ? -1 : k);
    jPrev = (j == (0) ? L - 2 : j - 1);
    kPrev = (k == (0) ? -1 : k - 1);
    JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
    JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

    if (kNext != -1) {
        res += JInterlayerUp * S(xface[i][jNext][kNext].z) 
            * (xface_up_0[i][jNext][kNext] - xface_down_0[i][jNext][kNext]);
        res += JInterlayerUp * S(xface[i][jPrev][kNext].z) 
            * (xface_up_0[i][jPrev][kNext] - xface_down_0[i][jPrev][kNext]);
    }

    if (kPrev != -1) {
        res += JInterlayerDown * S(xface[i][jPrev][kPrev].z) 
            * (xface_up_0[i][jPrev][kPrev] - xface_down_0[i][jPrev][kPrev]);
        res += JInterlayerDown * S(xface[i][jNext][kPrev].z) 
            * (xface_up_0[i][jNext][kPrev] - xface_down_0[i][jNext][kPrev]);
    }

    /*-------------- Neighbors for yFace ----------------*/
    iNext = (i == (L - 2) ? 0 : i);
    kNext = (k == (N - 2) ? -1 : k);
    iPrev = (i == (0) ? L - 2 : i - 1);
    kPrev = (k == (0) ? -1 : k - 1);

    if (kNext != -1) {
        res += JInterlayerUp * S(yface[iNext][j][kNext].z) 
            * (yface_up_0[iNext][j][kNext] - yface_down_0[iNext][j][kNext]);
        res += JInterlayerUp * S(yface[iPrev][j][kNext].z) 
            * (yface_up_0[iPrev][j][kNext] - yface_down_0[iPrev][j][kNext]);
    }

    if (kPrev != -1) {
        res += JInterlayerDown * S(yface[iPrev][j][kPrev].z) 
            * (yface_up_0[iPrev][j][kPrev] - yface_down_0[iPrev][j][kPrev]);
        res += JInterlayerDown * S(yface[iNext][j][kPrev].z) 
            * (yface_up_0[iNext][j][kPrev] - yface_down_0[iNext][j][kPrev]);
    }

    /*-------------- Neighbours for zFace ----------------*/
    iNext = (i == (L - 2) ? 0 : i);
    jNext = (j == (L - 2) ? 0 : j);
    iPrev = (i == (0) ? L - 2 : i - 1);
    jPrev = (j == (0) ? L - 2 : j - 1);

    res += S(zface[iNext][jNext][k].z) 
        * (zface_up_0[iNext][jNext][k] - zface_down_0[iNext][jNext][k]);
    res += S(zface[iPrev][jNext][k].z) 
        * (zface_up_0[iPrev][jNext][k] - zface_down_0[iPrev][jNext][k]);

    res += S(zface[iPrev][jPrev][k].z) 
        * (zface_up_0[iPrev][jPrev][k] - zface_down_0[iPrev][jPrev][k]);
    res += S(zface[iNext][jPrev][k].z) 
        * (zface_up_0[iNext][jPrev][k] - zface_down_0[iNext][jPrev][k]);

    /*------------ Neighbours of second order ---------------*/
    if (i == 0) {
        res += J_second * S(simple[L - 1][j][k].z) 
            * (simple_up_0[L - 1][j][k] - simple_down_0[L - 1][j][k]);
    } else {
        res += J_second * S(simple[i - 1][j][k].z) 
            * (simple_up_0[i - 1][j][k] - simple_down_0[i - 1][j][k]);
    }

    if (i == L - 1) {
        res += J_second * S(simple[0][j][k].z) 
            * (simple_up_0[0][j][k] - simple_down_0[0][j][k]);
    } else {
        res += J_second * S(simple[i + 1][j][k].z) 
            * (simple_up_0[i + 1][j][k] - simple_down_0[i + 1][j][k]);
    }

    if (j == 0) {
        res += J_second * S(simple[i][L - 1][k].z) 
            * (simple_up_0[i][L - 1][k] - simple_down_0[i][L - 1][k]);
    } else {
        res += J_second * S(simple[i][j - 1][k].z) 
            * (simple_up_0[i][j - 1][k] - simple_down_0[i][j - 1][k]);
    }

    if (j == L - 1) {
        res += J_second * S(simple[i][0][k].z) 
            * (simple_up_0[i][0][k] - simple_down_0[i][0][k]);
    } else {
        res += J_second * S(simple[i][j + 1][k].z) 
            * (simple_up_0[i][j + 1][k] - simple_down_0[i][j + 1][k]);
    }

    if (k != 0 && move != 1) {
        if (k == N / 2) {
            res += J2 * S(simple[i][j][k - 1].z) 
                * (simple_up_0[i][j][k - 1] - simple_down_0[i][j][k - 1]);
        } else {
            res += S(simple[i][j][k - 1].z) 
                * (simple_up_0[i][j][k - 1] - simple_down_0[i][j][k - 1]);
        }
    }

    if (k != N - 2) {
        if (k == N / 2 - 1) {
            res += J2 * S(simple[i][j][k + 1].z) 
                * (simple_up_0[i][j][k + 1] - simple_down_0[i][j][k + 1]);
        } else {
            res += S(simple[i][j][k + 1].z) 
                * (simple_up_0[i][j][k + 1] - simple_down_0[i][j][k + 1]);
        }
    }

    return res;
}

double GMRFaceEnergySumm(const heis_spin *const *const *const simple, double ***simple_up, double ***simple_down, 
                         const heis_spin *const *const *const xface, double ***xface_up, double ***xface_down, 
                         const heis_spin *const *const *const yface, double ***yface_up, double ***yface_down, 
                         const heis_spin *const *const *const zface, double ***zface_up, double ***zface_down, 
                         int i, int j, int k, int L, int N, double J2, double J_second, int move, char type) {
    double res = 0.0;

    int iNext, jNext, kNext;
    int iPrev, jPrev, kPrev;
    double JInterlayerUp, JInterlayerDown;

    switch (type) {
        case 'x':
            // ----------------------- first order ----------------------------------//
            iNext = (i != L - 2) ? i : 0;

            res += S(yface[iNext][j][k].z) 
                * (yface_up[iNext][j][k] - yface_down[iNext][j][k]);
            res += S(yface[iNext][j + 1][k].z) 
                * (yface_up[iNext][j + 1][k] - yface_down[iNext][j + 1][k]);

            res += S(zface[iNext][j][k].z) 
                * (zface_up[iNext][j][k] - zface_down[iNext][j][k]);
            res += S(zface[iNext][j + 1][k].z) 
                * (zface_up[iNext][j + 1][k] - zface_down[iNext][j + 1][k]);

            iPrev = (i != 0) ? i - 1 : L - 2;

            res += S(yface[iPrev][j][k].z) 
                * (yface_up[iPrev][j][k] - yface_down[iPrev][j][k]);
            res += S(yface[iPrev][j + 1][k].z) 
                * (yface_up[iPrev][j + 1][k] - yface_down[iPrev][j + 1][k]);

            res += S(zface[iPrev][j][k].z) 
                * (zface_up[iPrev][j][k] - zface_down[iPrev][j][k]);
            res += S(zface[iPrev][j + 1][k].z) 
                * (zface_up[iPrev][j + 1][k] - zface_down[iPrev][j + 1][k]);

            res += S(simple[i][j][k].z) 
                * (simple_up[i][j][k] - simple_down[i][j][k]);
            res += S(simple[i][j + 1][k].z) 
                * (simple_up[i][j + 1][k] - simple_down[i][j + 1][k]);
            res += S(simple[i][j][k + 1].z) 
                * (simple_up[i][j][k + 1] - simple_down[i][j][k + 1]);
            res += S(simple[i][j + 1][k + 1].z) 
                * (simple_up[i][j + 1][k + 1] - simple_down[i][j + 1][k + 1]);


            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * S(xface[iNext][j][k].z) 
                * (xface_up[iNext][j][k] - xface_down[iNext][j][k]);
            res += J_second * S(xface[iPrev][j][k].z) 
                * (xface_up[iPrev][j][k] - xface_down[iPrev][j][k]);
            res += J_second * S(xface[i][jNext][k].z) 
                * (xface_up[i][jNext][k] - xface_down[i][jNext][k]);
            res += J_second * S(xface[i][jPrev][k].z) 
                * (xface_up[i][jPrev][k] - xface_down[i][jPrev][k]);

            if (kNext != -1) {
                res += JInterlayerUp * J_second * S(xface[i][j][kNext].z) 
                    * (xface_up[i][j][kNext] - xface_down[i][j][kNext]);
            }

            if (kPrev != -1 && move != -1) {
                res += JInterlayerDown * J_second * S(xface[i][j][kPrev].z) 
                    * (xface_up[i][j][kPrev] - xface_down[i][j][kPrev]);
            }

            break;


        case 'y':
            // ----------------------- first order ----------------------------------//

            jNext = (j != L - 2) ? j : 0;

            res += S(xface[i][jNext][k].z) 
                * (xface_up[i][jNext][k] - xface_down[i][jNext][k]);
            res += S(xface[i + 1][jNext][k].z) 
                * (xface_up[i + 1][jNext][k] - xface_down[i][jNext][k]);

            res += S(zface[i][jNext][k].z) 
                * (zface_up[i][jNext][k] - zface_down[i][jNext][k]);
            res += S(zface[i + 1][jNext][k].z) 
                * (zface_up[i + 1][jNext][k] - zface_down[i + 1][jNext][k]);

            jPrev = (j != 0) ? j - 1 : L - 2;

            res += S(xface[i][jPrev][k].z) 
                * (xface_up[i][jPrev][k] - xface_down[i][jPrev][k]);
            res += S(xface[i + 1][jPrev][k].z) 
                * (xface_up[i + 1][jPrev][k] - xface_down[i][jPrev][k]);

            res += S(zface[i][jPrev][k].z) 
                * (zface_up[i][jPrev][k] - zface_down[i][jPrev][k]);
            res += S(zface[i + 1][jPrev][k].z) 
                * (zface_up[i + 1][jPrev][k] - zface_down[i + 1][jPrev][k]);

            res += S(simple[i][j][k].z) 
                * (simple_up[i][j][k] - simple_down[i][j][k]);
            res += S(simple[i + 1][j][k].z) 
                * (simple_up[i + 1][j][k] - simple_down[i + 1][j][k]);
            res += S(simple[i][j][k + 1].z) 
                * (simple_up[i][j][k + 1] - simple_down[i][j][k + 1]);
            res += S(simple[i + 1][j][k + 1].z) 
                * (simple_up[i + 1][j][k + 1] - simple_down[i + 1][j][k + 1]);


            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * S(yface[iNext][j][k].z) 
                * (yface_up[iNext][j][k] - yface_down[iNext][j][k]);
            res += J_second * S(yface[iPrev][j][k].z) 
                * (yface_up[iPrev][j][k] - yface_down[iPrev][j][k]);
            res += J_second * S(yface[i][jNext][k].z) 
                * (yface_up[i][jNext][k] - yface_down[i][jNext][k]);
            res += J_second * S(yface[i][jPrev][k].z) 
                * (yface_up[i][jPrev][k] - yface_down[i][jPrev][k]);

            if (kNext != -1 && move != 1) {
                res += JInterlayerUp * J_second * S(yface[i][j][kNext].z) 
                    * (yface_up[i][j][kNext] - yface_down[i][j][kNext]);
            }

            if (kPrev != -1) {
                res += JInterlayerDown * J_second * S(yface[i][j][kPrev].z) 
                    * (yface_up[i][j][kPrev] - yface_down[i][j][kPrev]);
            }

            break;

        case 'z':
            // ----------------------- first order ----------------------------------//

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            if (k != N - 2) {
                res += JInterlayerUp * S(xface[i][j][k].z) 
                    * (xface_up[i][j][k] - xface_down[i][j][k]);
                res += JInterlayerUp * S(xface[i + 1][j][k].z) 
                    * (xface_up[i + 1][j][k] - xface_down[i + 1][j][k]);
                res += JInterlayerUp * S(yface[i][j][k].z) 
                    * (yface_up[i][j][k] - yface_down[i][j][k]);
                res += JInterlayerUp * S(yface[i][j + 1][k].z) 
                    * (yface_up[i][j + 1][k] - yface_down[i][j + 1][k]);
            }

            if (k != 0) {
                res += JInterlayerUp * S(xface[i][j][k - 1].z) 
                    * (xface_up[i][j][k - 1] - xface_down[i][j][k - 1]);
                res += JInterlayerUp * S(xface[i + 1][j][k - 1].z) 
                    * (xface_up[i + 1][j][k - 1] - xface_down[i + 1][j][k - 1]);
                res += JInterlayerUp * S(yface[i][j][k - 1].z) 
                    * (yface_up[i][j][k - 1] - yface_down[i][j][k - 1]);
                res += JInterlayerUp * S(yface[i][j + 1][k - 1].z) 
                    * (yface_up[i][j + 1][k - 1] - yface_down[i][j + 1][k - 1]);
            }

            res += S(simple[i][j][k].z) 
                * (simple_up[i][j][k] - simple_down[i][j][k]);
            res += S(simple[i + 1][j][k].z) 
                * (simple_up[i + 1][j][k] - simple_down[i + 1][j][k]);
            res += S(simple[i][j + 1][k].z) 
                * (simple_up[i][j + 1][k] - simple_down[i][j + 1][k]);
            res += S(simple[i + 1][j + 1][k].z) 
                * (simple_up[i + 1][j + 1][k] - simple_down[i + 1][j + 1][k]);

            /*------------------------second order -------------------------*/

            iNext = (i == (L - 2) ? 0 : i + 1);
            iPrev = (i == (0) ? L - 2 : i - 1);

            jNext = (j == (L - 2) ? 0 : j + 1);
            jPrev = (j == (0) ? L - 2 : j - 1);

            kNext = (k == (N - 2) ? -1 : k + 1);
            kPrev = (k == (0) ? -1 : k - 1);

            JInterlayerUp = (k == N / 2 - 2 ? J2 : 1.0);
            JInterlayerDown = (k == N / 2 - 1 ? J2 : 1.0);

            res += J_second * S(zface[iNext][j][k].z) 
                * (zface_up[iNext][j][k] - zface_down[iNext][j][k]);
            res += J_second * S(zface[iPrev][j][k].z) 
                * (zface_up[iPrev][j][k] - zface_down[iPrev][j][k]);
            res += J_second * S(zface[i][jNext][k].z) 
                * (zface_up[i][jNext][k] - zface_down[i][jNext][k]);
            res += J_second * S(zface[i][jPrev][k].z) 
                * (zface_up[i][jPrev][k] - zface_down[i][jPrev][k]);

            if (kNext != -1) {
                res += JInterlayerUp * J_second * S(zface[i][j][kNext].z) 
                    * (zface_up[i][j][kNext] - zface_down[i][j][kNext]);
            }

            if (kPrev != -1 && move != 1) {
                res += JInterlayerDown * J_second * S(zface[i][j][kPrev].z) 
                    * (zface_up[i][j][kPrev] - zface_down[i][j][kPrev]);
            }

            break;

        default:
            break;
    }

    return res;
}

void spin_transport(const heis_spin *const *const *const simple, double ***simple_up, double ***simple_down, 
                    const heis_spin *const *const *const xface, double ***xface_up, double ***xface_down, 
                    const heis_spin *const *const *const yface, double ***yface_up, double ***yface_down, 
                    const heis_spin *const *const *const zface, double ***zface_up, double ***zface_down, 
                    int L, int N, double J2, double J_second, double T_freezing, 
                    double norm, double *j_up, double *j_down) {
    double E1, E2;
    int i, j, k;
    int a;
    double delta_e, r;

    random_t gen_rand{spec_random::get_seed()};

    for (a = 0; a < 4 * L * L * N; ++a) {
        i = gen_rand(0, L);
        j = gen_rand(0, L);
        k = gen_rand(0, 4 * N - 1);
        if (k < N) {
            E1 = E2 = 0.0;
            E1 = GMRSimpleEnergySumm(simple, simple_up, simple_down, 
                                     xface, xface_up, xface_down, 
                                     yface, yface_up, yface_down, 
                                     zface, zface_up, zface_down, 
                                     i, j, k, L, N, J2, J_second, 0) 
                - S(simple[i][j][k].z) * (simple_up[i][j][k] - simple_down[i][j][k]);
            if (k == N - 1) {
                E2 = 0.0;
            } else {
                E2 = GMRSimpleEnergySumm(simple, simple_up, simple_down, 
                                         xface, xface_up, xface_down, 
                                         yface, yface_up, yface_down, 
                                         zface, zface_up, zface_down, 
                                         i, j, k, L, N, J2, J_second, 1) 
                    - S(simple[i][j][k + 1].z) * (simple_up[i][j][k] - simple_down[i][j][k] 
                    + simple_up[i][j][k + 1] - simple_down[i][j][k + 1]);
            }

            delta_e = E2 - E1;

            if (delta_e < 0.0) {
                *j_up += simple_up[i][j][k] * norm;
                *j_down += simple_down[i][j][k] * norm;

                simple_up[i][j][k + 1] += simple_up[i][j][k];
                simple_down[i][j][k + 1] += simple_down[i][j][k];
                simple_up[i][j][k] = 0.0;
                simple_down[i][j][k] = 0.0;
            } else {
                r = gen_rand();
                if (exp(-delta_e / T_freezing) > r) {
                    *j_up += simple_up[i][j][k] * norm;
                    *j_down += simple_down[i][j][k] * norm;

                    simple_up[i][j][k + 1] += simple_up[i][j][k];
                    simple_down[i][j][k + 1] += simple_down[i][j][k];
                    simple_up[i][j][k] = 0.0;
                    simple_down[i][j][k] = 0.0;
                }
            }
        }
        if (k >= N && k < 2 * N - 1) {
            k -= N;

            if (k != 0) {
                k -= 1;
            }
            if (i != 0) {
                i -= 1;
            }
            if (j != 0) {
                j -= 1;
            }

            E1 = E2 = 0.0;

            E1 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                                   xface, xface_up, xface_down, 
                                   yface, yface_up, yface_down, 
                                   zface, zface_up, zface_down, 
                                   i, j, k, L, N, J2, J_second, 0, 'x') 
                - S(xface[i][j][k].z) * (xface_up[i][j][k] - xface_down[i][j][k]);
            if (k == N - 1) {
                E2 = 0.0;
            } else {
                E2 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                                       xface, xface_up, xface_down, 
                                       yface, yface_up, yface_down, 
                                       zface, zface_up, zface_down, 
                                       i, j, k, L, N, J2, J_second, 1, 'x') 
                    - S(xface[i][j][k + 1].z) * (xface_up[i][j][k] - xface_down[i][j][k] 
                    + xface_up[i][j][k + 1] - xface_down[i][j][k + 1]);
            }

            delta_e = E2 - E1;

            if (delta_e < 0.0) {
                *j_up += xface_up[i][j][k] * norm;
                *j_down += xface_down[i][j][k] * norm;

                xface_up[i][j][k + 1] += xface_up[i][j][k];
                xface_down[i][j][k + 1] += xface_down[i][j][k];
                xface_up[i][j][k] = 0.0;
                xface_down[i][j][k] = 0.0;
            } else {
                r = gen_rand();
                if (exp(-delta_e / T_freezing) > r) {
                    *j_up += xface_up[i][j][k] * norm;
                    *j_down += xface_down[i][j][k] * norm;

                    xface_up[i][j][k + 1] += xface_up[i][j][k];
                    xface_down[i][j][k + 1] += xface_down[i][j][k];
                    xface_up[i][j][k] = 0.0;
                    xface_down[i][j][k] = 0.0;
                }
            }
        }

        if (k >= 2 * N && k < 3 * N - 1) {
            k -= 2 * N;

            if (k != 0) {
                k -= 1;
            }
            if (i != 0) {
                i -= 1;
            }
            if (j != 0) {
                j -= 1;
            }

            E1 = E2 = 0.0;

            E1 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                                   xface, xface_up, xface_down, 
                                   yface, yface_up, yface_down, 
                                   zface, zface_up, zface_down, 
                                   i, j, k, L, N, J2, J_second, 0, 'y') 
                - S(yface[i][j][k].z) * (yface_up[i][j][k] - yface_down[i][j][k]);
            if (k == N - 1) {
                E2 = 0.0;
            } else {
                E2 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                                       xface, xface_up, xface_down, 
                                       yface, yface_up, yface_down, 
                                       zface, zface_up, zface_down, 
                                       i, j, k, L, N, J2, J_second, 1, 'y') 
                    - S(yface[i][j][k + 1].z) * (yface_up[i][j][k] - yface_down[i][j][k] 
                    + yface_up[i][j][k + 1] - yface_down[i][j][k + 1]);
            }

            delta_e = E2 - E1;

            if (delta_e < 0.0) {
                *j_up += yface_up[i][j][k] * norm;
                *j_down += yface_down[i][j][k] * norm;

                yface_up[i][j][k + 1] += yface_up[i][j][k];
                yface_down[i][j][k + 1] += yface_down[i][j][k];
                yface_up[i][j][k] = 0.0;
                yface_down[i][j][k] = 0.0;
            } else {
                r = gen_rand();
                if (exp(-delta_e / T_freezing) > r) {
                    *j_up += yface_up[i][j][k] * norm;
                    *j_down += yface_down[i][j][k] * norm;

                    yface_up[i][j][k + 1] += yface_up[i][j][k];
                    yface_down[i][j][k + 1] += yface_down[i][j][k];
                    yface_up[i][j][k] = 0.0;
                    yface_down[i][j][k] = 0.0;
                }
            }
        }
    }

    if (k >= 3 * N) {
        k -= 3 * N;

        if (k != 0) {
            k -= 1;
        }
        if (i != 0) {
            i -= 1;
        }
        if (j != 0) {
            j -= 1;
        }

        E1 = E2 = 0.0;

        E1 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                               xface, xface_up, xface_down, 
                               yface, yface_up, yface_down,
                               zface, zface_up, zface_down, 
                               i, j, k, L, N, J2, J_second, 0, 'z') 
            - S(zface[i][j][k].z) * (zface_up[i][j][k] - zface_down[i][j][k]);
        if (k == N - 1) {
            E2 = 0.0;
        } else {
            E2 = GMRFaceEnergySumm(simple, simple_up, simple_down, 
                                   xface, xface_up, xface_down, 
                                   yface, yface_up, yface_down, 
                                   zface, zface_up, zface_down, 
                                   i, j, k, L, N, J2, J_second, 1, 'z') 
                - S(zface[i][j][k + 1].z) * (zface_up[i][j][k] - zface_down[i][j][k] 
                + zface_up[i][j][k + 1] - zface_down[i][j][k + 1]);
        }

        delta_e = E2 - E1;

        if (delta_e < 0.0) {
            *j_up += zface_up[i][j][k] * norm;
            *j_down += zface_down[i][j][k] * norm;

            zface_up[i][j][k + 1] += zface_up[i][j][k];
            zface_down[i][j][k + 1] += zface_down[i][j][k];
            zface_up[i][j][k] = 0.0;
            zface_down[i][j][k] = 0.0;
        } else {
            r = gen_rand();
            if (exp(-delta_e / T_freezing) > r) {
                *j_up += zface_up[i][j][k] * norm;
                *j_down += zface_down[i][j][k] * norm;

                zface_up[i][j][k + 1] += zface_up[i][j][k];
                zface_down[i][j][k + 1] += zface_down[i][j][k];
                zface_up[i][j][k] = 0.0;
                zface_down[i][j][k] = 0.0;
            }
        }
    }
}

#endif