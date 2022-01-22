#ifndef _HELPER_H_
#define _HELPER_H_

#include <stdint.h>

// Simulation constants
#define R_STAR              3.0
#define EPSILON_STAR        0.2
#define L                   30.0
#define N_SYM               27
#define TOLERANCE           1.0e-7
#define DT                  1.0e-15
#define FORCE_CONVERSION    4.186e-4
#define FORCE_CONVERSION_x2 8.372e-4
#define R_CONSTANT          1.99e-3
#define M_I                 18.0
#define T_0                 300.0

// General constant
#define ALIGN 64

// Maths macros
#define abs_double(x) (x < 0.0 ? -x : x)

#define square(x) ((x) * (x))
#define cube(x)   ((x) * (x) * (x))
#define quad(x)   ((x) * (x) * (x) * (x))
#define penta(x)  ((x) * (x) * (x) * (x) * (x))
#define hexa(x)   ((x) * (x) * (x) * (x) * (x) * (x))
#define septa(x)  ((x) * (x) * (x) * (x) * (x) * (x) * (x))

// Global variable
extern uint64_t N_PARTICLES_TOTAL;
extern uint64_t N_PARTICLES_LOCAL;
extern uint64_t LOCAL_EQUAL_TOTAL;
extern uint64_t N_DL;

enum
  {
    ERR_NONE,
    ERR_USAGE,
    ERR_OPEN
  };

struct particle
{
  double x;
  double y;
  double z;
};

struct force
{
  double fx;
  double fy;
  double fz;
};

struct lennard_jones
{
  double energy;
  struct force **restrict f;
  struct force *restrict sum_i;
  struct force *restrict sum;
};

struct translation_vector
{
  double x;
  double y;
  double z;
};

struct kinetic_moment
{
  double px;
  double py;
  double pz;
};

struct ket
{
  double kinetic_energy;
  double temperature;
};

#endif // _HELPER_H_
