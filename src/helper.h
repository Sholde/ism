#ifndef _HELPER_H_
#define _HELPER_H_

#include <stdint.h>

// Simulation constants
#define R_STAR       3.0
#define EPSILON_STAR 0.2
#define L            30.0
#define N_SYM        27
#define TOLERANCE    1.0e-7

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
};

struct translation_vector
{
  double x;
  double y;
  double z;
};

#endif // _HELPER_H_
