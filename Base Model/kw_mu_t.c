#include "udf.h"

#define Cmu 0.0845
#define rho_V 0.0261
#define rho_L 996.57

DEFINE_TURBULENT_VISCOSITY(mut_rev, c, t)
{
    real mu_t;
    real f_rho;
    real gamma;

    gamma = (C_R(c,t)-rho_V)/(rho_L-rho_V);
    f_rho = rho_V+pow(gamma,10)*(rho_L-rho_V);
    mu_t = f_rho*Cmu*C_K(c,t)/C_O(c,t);

    return mu_t;
}
