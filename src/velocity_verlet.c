#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "helper.h"
#include "lennard_jones.h"
#include "velocity_verlet.h"

// x if y >= 0.0, -x else
#define sign_function(x, y) (y < 0.0 ? -x : x)

// Initialize seed
static inline void init_random(void)
{
  srand(time(NULL));
}

struct ket *init_ket(void)
{
  // Allocate return structure
  struct ket *restrict ket = aligned_alloc(ALIGN, sizeof(struct ket));

  ket->kinetic_energy = 0.0;
  ket->temperature = 0.0;

  return ket;
}

void compute_kinetic_energy_and_temperature(struct ket *restrict ket,
                                            const struct kinetic_moment
                                            *restrict km)
{
  // Kinetic energy
  ket->kinetic_energy = 0.0;

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      ket->kinetic_energy += square(km[i].px) + square(km[i].py) + square(km[i].pz);
    }

  ket->kinetic_energy /= (M_I * FORCE_CONVERSION_x2);

  // Temperature
  ket->temperature = ket->kinetic_energy / (N_DL * R_CONSTANT);
}

void free_ket(struct ket *restrict ket)
{
  free(ket);
}

static void first_recalibration(struct kinetic_moment *restrict km)
{
  // Compute kinetic energy and temperature for this kinetic moment
  struct ket *restrict ket = init_ket();
  compute_kinetic_energy_and_temperature(ket, km);

  // Recalibration
  double rapport = sqrt((N_DL * R_CONSTANT * T_0) / ket->kinetic_energy);

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      km[i].px *= rapport;
      km[i].py *= rapport;
      km[i].pz *= rapport;
    }

  // Release memroy
  free_ket(ket);
}

static void second_recalibration(struct kinetic_moment *restrict km)
{
  struct kinetic_moment sum_km =
    {
      .px = 0.0,
      .py = 0.0,
      .pz = 0.0
    };

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      sum_km.px += km[i].px;
      sum_km.py += km[i].py;
      sum_km.pz += km[i].pz;
    }

  sum_km.px /= N_PARTICLES_TOTAL;
  sum_km.py /= N_PARTICLES_TOTAL;
  sum_km.pz /= N_PARTICLES_TOTAL;

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      km[i].px -= sum_km.px;
      km[i].py -= sum_km.py;
      km[i].pz -= sum_km.pz;
    }
}

struct kinetic_moment *init_velocity_verlet(void)
{
  // Set the number of degree of liberty
  N_DL = 3 * N_PARTICLES_TOTAL - 3;

  // Initial kinetic moment generation
  struct kinetic_moment *restrict km =
    aligned_alloc(ALIGN, sizeof(struct kinetic_moment) * N_PARTICLES_TOTAL);

  double c = 0.0;
  double s = 0.0;

  init_random();

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      // x
      c = (double)rand() / (double)RAND_MAX;
      s = (double)rand() / (double)RAND_MAX;
      km[i].px = sign_function(1.0, 0.5 - s) * c;

      // y
      c = (double)rand() / (double)RAND_MAX;
      s = (double)rand() / (double)RAND_MAX;
      km[i].py = sign_function(1.0, 0.5 - s) * c;

      // z
      c = (double)rand() / (double)RAND_MAX;
      s = (double)rand() / (double)RAND_MAX;
      km[i].pz = sign_function(1.0, 0.5 - s) * c;
    }

  // Recalibration
  first_recalibration(km);
  second_recalibration(km);
  first_recalibration(km);

  return km;
}

void free_kinetic_moment(struct kinetic_moment *restrict km)
{
  free(km);
}

void velocity_verlet(struct particle *restrict p,
                     __attribute__ ((unused)) struct translation_vector *restrict tv,
                     struct lennard_jones *restrict plj,
                     struct kinetic_moment *restrict km,
                     __attribute__ ((unused)) const double r_cut)
{
  // Compute forces
#if CLASSICAL
  lennard_jones(plj, p);
#elif PERIODICAL
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);
#else
  lennard_jones(plj, p);
#endif

  // Update kinetic moments
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      km[i].px -= DT * FORCE_CONVERSION * plj->sum_i[i].fx * 0.5;
      km[i].py -= DT * FORCE_CONVERSION * plj->sum_i[i].fy * 0.5;
      km[i].pz -= DT * FORCE_CONVERSION * plj->sum_i[i].fz * 0.5;
    }

  // Update positions
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      p[i].x += DT * km[i].px / M_I;
      p[i].y += DT * km[i].py / M_I;
      p[i].z += DT * km[i].pz / M_I;
    }

  // Re-compute forces
#if CLASSICAL
  lennard_jones(plj, p);
#elif PERIODICAL
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);
#else
  lennard_jones(plj, p);
#endif

  // Update kinetic moments
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      km[i].px -= DT * FORCE_CONVERSION * plj->sum_i[i].fx * 0.5;
      km[i].py -= DT * FORCE_CONVERSION * plj->sum_i[i].fy * 0.5;
      km[i].pz -= DT * FORCE_CONVERSION * plj->sum_i[i].fz * 0.5;
    }
}

void berendsen_thermostat(struct kinetic_moment *restrict km,
                          struct ket *restrict ket)
{
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      km[i].px += km[i].px * GAMMA * (ket->temperature / T_0 - 1);
      km[i].py += km[i].py * GAMMA * (ket->temperature / T_0 - 1);
      km[i].pz += km[i].pz * GAMMA * (ket->temperature / T_0 - 1);
    }
}
