#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6175132429151416991);
void live_err_fun(double *nom_x, double *delta_x, double *out_7951279187815065907);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_637224016576002731);
void live_H_mod_fun(double *state, double *out_8696275111750558888);
void live_f_fun(double *state, double dt, double *out_1735690330431513099);
void live_F_fun(double *state, double dt, double *out_1808211745531574140);
void live_h_4(double *state, double *unused, double *out_4098430510588616918);
void live_H_4(double *state, double *unused, double *out_2823963131679719256);
void live_h_9(double *state, double *unused, double *out_4534930155945164839);
void live_H_9(double *state, double *unused, double *out_3065152778309309901);
void live_h_10(double *state, double *unused, double *out_8965664941676640835);
void live_H_10(double *state, double *unused, double *out_8597485743186393706);
void live_h_12(double *state, double *unused, double *out_4431967718199545297);
void live_H_12(double *state, double *unused, double *out_7843419539711681051);
void live_h_35(double *state, double *unused, double *out_1708372239267622899);
void live_H_35(double *state, double *unused, double *out_7857761501672856856);
void live_h_32(double *state, double *unused, double *out_3490630552210583304);
void live_H_32(double *state, double *unused, double *out_5832391986255021319);
void live_h_13(double *state, double *unused, double *out_6858976564956992067);
void live_H_13(double *state, double *unused, double *out_2853644333017370113);
void live_h_14(double *state, double *unused, double *out_4534930155945164839);
void live_H_14(double *state, double *unused, double *out_3065152778309309901);
void live_h_33(double *state, double *unused, double *out_6435807191173586095);
void live_H_33(double *state, double *unused, double *out_4707204497033999252);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}