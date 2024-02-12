#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1996390967112320251) {
   out_1996390967112320251[0] = delta_x[0] + nom_x[0];
   out_1996390967112320251[1] = delta_x[1] + nom_x[1];
   out_1996390967112320251[2] = delta_x[2] + nom_x[2];
   out_1996390967112320251[3] = delta_x[3] + nom_x[3];
   out_1996390967112320251[4] = delta_x[4] + nom_x[4];
   out_1996390967112320251[5] = delta_x[5] + nom_x[5];
   out_1996390967112320251[6] = delta_x[6] + nom_x[6];
   out_1996390967112320251[7] = delta_x[7] + nom_x[7];
   out_1996390967112320251[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2569961913780748397) {
   out_2569961913780748397[0] = -nom_x[0] + true_x[0];
   out_2569961913780748397[1] = -nom_x[1] + true_x[1];
   out_2569961913780748397[2] = -nom_x[2] + true_x[2];
   out_2569961913780748397[3] = -nom_x[3] + true_x[3];
   out_2569961913780748397[4] = -nom_x[4] + true_x[4];
   out_2569961913780748397[5] = -nom_x[5] + true_x[5];
   out_2569961913780748397[6] = -nom_x[6] + true_x[6];
   out_2569961913780748397[7] = -nom_x[7] + true_x[7];
   out_2569961913780748397[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4814550263766126742) {
   out_4814550263766126742[0] = 1.0;
   out_4814550263766126742[1] = 0;
   out_4814550263766126742[2] = 0;
   out_4814550263766126742[3] = 0;
   out_4814550263766126742[4] = 0;
   out_4814550263766126742[5] = 0;
   out_4814550263766126742[6] = 0;
   out_4814550263766126742[7] = 0;
   out_4814550263766126742[8] = 0;
   out_4814550263766126742[9] = 0;
   out_4814550263766126742[10] = 1.0;
   out_4814550263766126742[11] = 0;
   out_4814550263766126742[12] = 0;
   out_4814550263766126742[13] = 0;
   out_4814550263766126742[14] = 0;
   out_4814550263766126742[15] = 0;
   out_4814550263766126742[16] = 0;
   out_4814550263766126742[17] = 0;
   out_4814550263766126742[18] = 0;
   out_4814550263766126742[19] = 0;
   out_4814550263766126742[20] = 1.0;
   out_4814550263766126742[21] = 0;
   out_4814550263766126742[22] = 0;
   out_4814550263766126742[23] = 0;
   out_4814550263766126742[24] = 0;
   out_4814550263766126742[25] = 0;
   out_4814550263766126742[26] = 0;
   out_4814550263766126742[27] = 0;
   out_4814550263766126742[28] = 0;
   out_4814550263766126742[29] = 0;
   out_4814550263766126742[30] = 1.0;
   out_4814550263766126742[31] = 0;
   out_4814550263766126742[32] = 0;
   out_4814550263766126742[33] = 0;
   out_4814550263766126742[34] = 0;
   out_4814550263766126742[35] = 0;
   out_4814550263766126742[36] = 0;
   out_4814550263766126742[37] = 0;
   out_4814550263766126742[38] = 0;
   out_4814550263766126742[39] = 0;
   out_4814550263766126742[40] = 1.0;
   out_4814550263766126742[41] = 0;
   out_4814550263766126742[42] = 0;
   out_4814550263766126742[43] = 0;
   out_4814550263766126742[44] = 0;
   out_4814550263766126742[45] = 0;
   out_4814550263766126742[46] = 0;
   out_4814550263766126742[47] = 0;
   out_4814550263766126742[48] = 0;
   out_4814550263766126742[49] = 0;
   out_4814550263766126742[50] = 1.0;
   out_4814550263766126742[51] = 0;
   out_4814550263766126742[52] = 0;
   out_4814550263766126742[53] = 0;
   out_4814550263766126742[54] = 0;
   out_4814550263766126742[55] = 0;
   out_4814550263766126742[56] = 0;
   out_4814550263766126742[57] = 0;
   out_4814550263766126742[58] = 0;
   out_4814550263766126742[59] = 0;
   out_4814550263766126742[60] = 1.0;
   out_4814550263766126742[61] = 0;
   out_4814550263766126742[62] = 0;
   out_4814550263766126742[63] = 0;
   out_4814550263766126742[64] = 0;
   out_4814550263766126742[65] = 0;
   out_4814550263766126742[66] = 0;
   out_4814550263766126742[67] = 0;
   out_4814550263766126742[68] = 0;
   out_4814550263766126742[69] = 0;
   out_4814550263766126742[70] = 1.0;
   out_4814550263766126742[71] = 0;
   out_4814550263766126742[72] = 0;
   out_4814550263766126742[73] = 0;
   out_4814550263766126742[74] = 0;
   out_4814550263766126742[75] = 0;
   out_4814550263766126742[76] = 0;
   out_4814550263766126742[77] = 0;
   out_4814550263766126742[78] = 0;
   out_4814550263766126742[79] = 0;
   out_4814550263766126742[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4065337694278811682) {
   out_4065337694278811682[0] = state[0];
   out_4065337694278811682[1] = state[1];
   out_4065337694278811682[2] = state[2];
   out_4065337694278811682[3] = state[3];
   out_4065337694278811682[4] = state[4];
   out_4065337694278811682[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4065337694278811682[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4065337694278811682[7] = state[7];
   out_4065337694278811682[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1883874506294394197) {
   out_1883874506294394197[0] = 1;
   out_1883874506294394197[1] = 0;
   out_1883874506294394197[2] = 0;
   out_1883874506294394197[3] = 0;
   out_1883874506294394197[4] = 0;
   out_1883874506294394197[5] = 0;
   out_1883874506294394197[6] = 0;
   out_1883874506294394197[7] = 0;
   out_1883874506294394197[8] = 0;
   out_1883874506294394197[9] = 0;
   out_1883874506294394197[10] = 1;
   out_1883874506294394197[11] = 0;
   out_1883874506294394197[12] = 0;
   out_1883874506294394197[13] = 0;
   out_1883874506294394197[14] = 0;
   out_1883874506294394197[15] = 0;
   out_1883874506294394197[16] = 0;
   out_1883874506294394197[17] = 0;
   out_1883874506294394197[18] = 0;
   out_1883874506294394197[19] = 0;
   out_1883874506294394197[20] = 1;
   out_1883874506294394197[21] = 0;
   out_1883874506294394197[22] = 0;
   out_1883874506294394197[23] = 0;
   out_1883874506294394197[24] = 0;
   out_1883874506294394197[25] = 0;
   out_1883874506294394197[26] = 0;
   out_1883874506294394197[27] = 0;
   out_1883874506294394197[28] = 0;
   out_1883874506294394197[29] = 0;
   out_1883874506294394197[30] = 1;
   out_1883874506294394197[31] = 0;
   out_1883874506294394197[32] = 0;
   out_1883874506294394197[33] = 0;
   out_1883874506294394197[34] = 0;
   out_1883874506294394197[35] = 0;
   out_1883874506294394197[36] = 0;
   out_1883874506294394197[37] = 0;
   out_1883874506294394197[38] = 0;
   out_1883874506294394197[39] = 0;
   out_1883874506294394197[40] = 1;
   out_1883874506294394197[41] = 0;
   out_1883874506294394197[42] = 0;
   out_1883874506294394197[43] = 0;
   out_1883874506294394197[44] = 0;
   out_1883874506294394197[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1883874506294394197[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1883874506294394197[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1883874506294394197[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1883874506294394197[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1883874506294394197[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1883874506294394197[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1883874506294394197[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1883874506294394197[53] = -9.8000000000000007*dt;
   out_1883874506294394197[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1883874506294394197[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1883874506294394197[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1883874506294394197[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1883874506294394197[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1883874506294394197[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1883874506294394197[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1883874506294394197[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1883874506294394197[62] = 0;
   out_1883874506294394197[63] = 0;
   out_1883874506294394197[64] = 0;
   out_1883874506294394197[65] = 0;
   out_1883874506294394197[66] = 0;
   out_1883874506294394197[67] = 0;
   out_1883874506294394197[68] = 0;
   out_1883874506294394197[69] = 0;
   out_1883874506294394197[70] = 1;
   out_1883874506294394197[71] = 0;
   out_1883874506294394197[72] = 0;
   out_1883874506294394197[73] = 0;
   out_1883874506294394197[74] = 0;
   out_1883874506294394197[75] = 0;
   out_1883874506294394197[76] = 0;
   out_1883874506294394197[77] = 0;
   out_1883874506294394197[78] = 0;
   out_1883874506294394197[79] = 0;
   out_1883874506294394197[80] = 1;
}
void h_25(double *state, double *unused, double *out_2991891326874373464) {
   out_2991891326874373464[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3365004724219061885) {
   out_3365004724219061885[0] = 0;
   out_3365004724219061885[1] = 0;
   out_3365004724219061885[2] = 0;
   out_3365004724219061885[3] = 0;
   out_3365004724219061885[4] = 0;
   out_3365004724219061885[5] = 0;
   out_3365004724219061885[6] = 1;
   out_3365004724219061885[7] = 0;
   out_3365004724219061885[8] = 0;
}
void h_24(double *state, double *unused, double *out_3699121416283704588) {
   out_3699121416283704588[0] = state[4];
   out_3699121416283704588[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2446380906803715080) {
   out_2446380906803715080[0] = 0;
   out_2446380906803715080[1] = 0;
   out_2446380906803715080[2] = 0;
   out_2446380906803715080[3] = 0;
   out_2446380906803715080[4] = 1;
   out_2446380906803715080[5] = 0;
   out_2446380906803715080[6] = 0;
   out_2446380906803715080[7] = 0;
   out_2446380906803715080[8] = 0;
   out_2446380906803715080[9] = 0;
   out_2446380906803715080[10] = 0;
   out_2446380906803715080[11] = 0;
   out_2446380906803715080[12] = 0;
   out_2446380906803715080[13] = 0;
   out_2446380906803715080[14] = 1;
   out_2446380906803715080[15] = 0;
   out_2446380906803715080[16] = 0;
   out_2446380906803715080[17] = 0;
}
void h_30(double *state, double *unused, double *out_5939617317229858425) {
   out_5939617317229858425[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3494343671362301955) {
   out_3494343671362301955[0] = 0;
   out_3494343671362301955[1] = 0;
   out_3494343671362301955[2] = 0;
   out_3494343671362301955[3] = 0;
   out_3494343671362301955[4] = 1;
   out_3494343671362301955[5] = 0;
   out_3494343671362301955[6] = 0;
   out_3494343671362301955[7] = 0;
   out_3494343671362301955[8] = 0;
}
void h_26(double *state, double *unused, double *out_4957533223894560313) {
   out_4957533223894560313[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7106508043093118109) {
   out_7106508043093118109[0] = 0;
   out_7106508043093118109[1] = 0;
   out_7106508043093118109[2] = 0;
   out_7106508043093118109[3] = 0;
   out_7106508043093118109[4] = 0;
   out_7106508043093118109[5] = 0;
   out_7106508043093118109[6] = 0;
   out_7106508043093118109[7] = 1;
   out_7106508043093118109[8] = 0;
}
void h_27(double *state, double *unused, double *out_3263933395084673790) {
   out_3263933395084673790[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5669106983162726866) {
   out_5669106983162726866[0] = 0;
   out_5669106983162726866[1] = 0;
   out_5669106983162726866[2] = 0;
   out_5669106983162726866[3] = 1;
   out_5669106983162726866[4] = 0;
   out_5669106983162726866[5] = 0;
   out_5669106983162726866[6] = 0;
   out_5669106983162726866[7] = 0;
   out_5669106983162726866[8] = 0;
}
void h_29(double *state, double *unused, double *out_977237319548565193) {
   out_977237319548565193[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7382469710032277899) {
   out_7382469710032277899[0] = 0;
   out_7382469710032277899[1] = 1;
   out_7382469710032277899[2] = 0;
   out_7382469710032277899[3] = 0;
   out_7382469710032277899[4] = 0;
   out_7382469710032277899[5] = 0;
   out_7382469710032277899[6] = 0;
   out_7382469710032277899[7] = 0;
   out_7382469710032277899[8] = 0;
}
void h_28(double *state, double *unused, double *out_988404577471621656) {
   out_988404577471621656[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5418839438466951648) {
   out_5418839438466951648[0] = 1;
   out_5418839438466951648[1] = 0;
   out_5418839438466951648[2] = 0;
   out_5418839438466951648[3] = 0;
   out_5418839438466951648[4] = 0;
   out_5418839438466951648[5] = 0;
   out_5418839438466951648[6] = 0;
   out_5418839438466951648[7] = 0;
   out_5418839438466951648[8] = 0;
}
void h_31(double *state, double *unused, double *out_4142443342045657196) {
   out_4142443342045657196[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3334358762342101457) {
   out_3334358762342101457[0] = 0;
   out_3334358762342101457[1] = 0;
   out_3334358762342101457[2] = 0;
   out_3334358762342101457[3] = 0;
   out_3334358762342101457[4] = 0;
   out_3334358762342101457[5] = 0;
   out_3334358762342101457[6] = 0;
   out_3334358762342101457[7] = 0;
   out_3334358762342101457[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1996390967112320251) {
  err_fun(nom_x, delta_x, out_1996390967112320251);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2569961913780748397) {
  inv_err_fun(nom_x, true_x, out_2569961913780748397);
}
void car_H_mod_fun(double *state, double *out_4814550263766126742) {
  H_mod_fun(state, out_4814550263766126742);
}
void car_f_fun(double *state, double dt, double *out_4065337694278811682) {
  f_fun(state,  dt, out_4065337694278811682);
}
void car_F_fun(double *state, double dt, double *out_1883874506294394197) {
  F_fun(state,  dt, out_1883874506294394197);
}
void car_h_25(double *state, double *unused, double *out_2991891326874373464) {
  h_25(state, unused, out_2991891326874373464);
}
void car_H_25(double *state, double *unused, double *out_3365004724219061885) {
  H_25(state, unused, out_3365004724219061885);
}
void car_h_24(double *state, double *unused, double *out_3699121416283704588) {
  h_24(state, unused, out_3699121416283704588);
}
void car_H_24(double *state, double *unused, double *out_2446380906803715080) {
  H_24(state, unused, out_2446380906803715080);
}
void car_h_30(double *state, double *unused, double *out_5939617317229858425) {
  h_30(state, unused, out_5939617317229858425);
}
void car_H_30(double *state, double *unused, double *out_3494343671362301955) {
  H_30(state, unused, out_3494343671362301955);
}
void car_h_26(double *state, double *unused, double *out_4957533223894560313) {
  h_26(state, unused, out_4957533223894560313);
}
void car_H_26(double *state, double *unused, double *out_7106508043093118109) {
  H_26(state, unused, out_7106508043093118109);
}
void car_h_27(double *state, double *unused, double *out_3263933395084673790) {
  h_27(state, unused, out_3263933395084673790);
}
void car_H_27(double *state, double *unused, double *out_5669106983162726866) {
  H_27(state, unused, out_5669106983162726866);
}
void car_h_29(double *state, double *unused, double *out_977237319548565193) {
  h_29(state, unused, out_977237319548565193);
}
void car_H_29(double *state, double *unused, double *out_7382469710032277899) {
  H_29(state, unused, out_7382469710032277899);
}
void car_h_28(double *state, double *unused, double *out_988404577471621656) {
  h_28(state, unused, out_988404577471621656);
}
void car_H_28(double *state, double *unused, double *out_5418839438466951648) {
  H_28(state, unused, out_5418839438466951648);
}
void car_h_31(double *state, double *unused, double *out_4142443342045657196) {
  h_31(state, unused, out_4142443342045657196);
}
void car_H_31(double *state, double *unused, double *out_3334358762342101457) {
  H_31(state, unused, out_3334358762342101457);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
