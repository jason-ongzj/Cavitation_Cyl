#include "udf.h"
#include "mem.h"

#define Cmu 0.0845		              /* Viscosity coefficient */
#define rho_V 0.0261      	        /* Density of water vapor */
#define rho_L 996.57                /* Density of water */
#define mu_L 0.000853               /* Dynamic viscosity of water */
#define mu_V 0.000017               /* Dynamic viscosity of water vapor */
#define D 0.048                     /* Diameter */
#define N 1.0E4                     /* Bubble number density */
#define VAPOR_PRESSURE 3560         /* Vapor pressure */
#define INLET_VELOCITY 8.92	        /* Inlet velocity */
#define PI 3.14159265               /* PI */
#define DELTA_T 1.0E-4              /* Timestep */
#define R_MIN 1.0E-6                /* Minimum bubble radius */

real void_fraction(cell_t c, Thread *mix_thread);
void compute_face_area(cell_t c, Thread *mix_thread, real cell_area);
void compute_source_term(cell_t c, Thread *mix_thread);

/* Calculation of void fraction f_g */
real void_fraction(cell_t c, Thread *mix_thread)
{
    real R_curr_timestep;
    real local_void_fraction;
    real void_fraction_lower_limit;

    R_curr_timestep = C_UDSI(c, mix_thread, 0);

    /* Invoke cavitation assumption */
    if (C_P(c,mix_thread) <= VAPOR_PRESSURE) {
        /* Void fraction calculation for local homogeneous model */
        local_void_fraction = 4 * N * PI * pow(R_curr_timestep, 2);
        void_fraction_lower_limit = 4 * N * PI * pow(R_MIN, 2);

        /* Void fraction limiter */
	      if (local_void_fraction > 0.95) {
	         local_void_fraction = 0.95;
	      }

        if (local_void_fraction < void_fraction_lower_limit) {
            local_void_fraction = void_fraction_lower_limit;
        }
    } else {
    	local_void_fraction = 0;
    }

    return local_void_fraction;
}

void compute_face_area(cell_t c, Thread *mix_thread, real cell_area)
{
    int n;
    face_t f;
    Thread *tf;
    real NV_VEC(face_area);

    c_face_loop(c, mix_thread, n){
        f = C_FACE(c, mix_thread, n);
        tf = C_FACE_THREAD(c, mix_thread, n);

        /* F_Area gives the area of the face relative to each direction in */
        /* vector form */
        F_AREA(face_area, f, tf);

        /* NV_MAG(face_area) gives the exact area of the face */
        cell_area += NV_MAG(face_area);
    }
}

void compute_source_term(cell_t c, Thread *mix_thread)
{
    real delta_P, rho_mixture;
    real rhs_denominator, source_term;
    real R_prev_timestep;

    real cell_area = 0;

    compute_face_area(c, mix_thread, cell_area);

    rho_mixture = C_R(c, mix_thread);
    R_prev_timestep = C_UDSI_M1(c, mix_thread, 0);

    if (C_P(c, mix_thread) <= VAPOR_PRESSURE){
        delta_P = VAPOR_PRESSURE - C_P(c, mix_thread);
        rhs_denominator = rho_L * (1.5 + cell_area * N * R_prev_timestep);

        /* Multiply source term in Kubota's model by rho_mixture to get source term */
        /* equivalent for the scalar transport equation in Fluent */
        source_term =  rho_mixture * sqrt(delta_P/rhs_denominator);
    } else {
        source_term = 0;
    }

    C_UDMI(c, mix_thread, 0) = source_term;
}

DEFINE_INIT(my_init_func, mix_domain)
{
    Thread *mix_thread;
    Thread *liquid_thread;
    Thread *vapor_thread;

    thread_loop_c(mix_thread, mix_domain)
    {
        cell_t c;
        liquid_thread = THREAD_SUB_THREAD(mix_thread, 0);
        vapor_thread = THREAD_SUB_THREAD(mix_thread, 1);
        begin_c_loop(c,mix_thread) {

            Thread *s_t;
            int phase_id;

            /* initialise bubble radius (R term) */
            C_UDSI(c,mix_thread,0) = 0.001;
            C_UDSI_M1(c,mix_thread,0) = 0;

            /* initalise source term for scalar transport equation */
            C_UDMI(c,mix_thread,0) = 0;

            /* loop through mixture threads */
            sub_thread_loop(s_t, mix_thread, phase_id) {
                if (phase_id == 0) {
                    C_VOF(c,liquid_thread) = 1;
                } else if (phase_id == 1) {
                    C_VOF(c,vapor_thread) = 1 - C_VOF(c,liquid_thread);
                }
            }
        } end_c_loop(c,mix_thread)
    }
}

/* Eddy viscosity model used in Fillipo's thesis */
DEFINE_TURBULENT_VISCOSITY(mut_rev, c, t)
{
    real mu_t;
    real f_rho;
    real gamma;

    gamma = (C_R(c,t) - rho_V) / (rho_L - rho_V);
    f_rho = rho_V + pow(gamma,10) * (rho_L - rho_V);
    mu_t = f_rho * Cmu *C_K(c,t) / C_O(c,t);

    return mu_t;
}

DEFINE_ADJUST(calc_vof, mix_domain)
{
    Thread *mix_thread;
    cell_t c;

    real source;

    Thread *vapor_thread;
    Thread *liquid_thread;

    thread_loop_c(mix_thread,mix_domain)
    {
        liquid_thread = THREAD_SUB_THREAD(mix_thread, 0);
        vapor_thread = THREAD_SUB_THREAD(mix_thread, 1);
    }

    mix_thread = THREAD_SUPER_THREAD(liquid_thread);

    if (Data_Valid_P()){
        begin_c_loop(c,mix_thread) {
            C_VOF(c,vapor_thread) = void_fraction(c, mix_thread);
            C_VOF(c, liquid_thread) = 1 - C_VOF(c, vapor_thread);

            compute_source_term(c, mix_thread);
        } end_c_loop(c,mix_thread)
    }
}

DEFINE_SOURCE(R_Source, c, t, dS, eqn)
{
    real source;

    source = C_UDMI(c,t,0);

    /* Explicit calculation of S */
    dS[eqn] = 0;
    return source;
}

DEFINE_DIFFUSIVITY(Zero_Diffusivity, c, t, i)
{
    /* Set to small value to avoid division by zero */
    real zero = 1E-10;
    return zero;
}

DEFINE_UDS_UNSTEADY(R_unsteady_term, c, t, i, apu, su)
{
    real rho_prev, rho_curr, R_prev_timestep;
    rho_curr = C_R(c,t);
    rho_prev = C_R_M1(c,t);
    R_prev_timestep = C_STORAGE_R(c, t, SV_UDSI_M1(i));

    *apu = -rho_curr/DELTA_T;
    *su = rho_prev*R_prev_timestep/DELTA_T;
}
