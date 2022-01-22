#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "helper.h"
#include "lennard_jones.h"
#include "velocity_verlet.h"

// x if y >= 0.0, -x else
#define sign_function(x, y) (y < 0.0 ? -x : x)

// Initialize seed
void init_random(void)
{
  srand(time(NULL));
}

struct ket *compute_kinetic_energy_and_temperature(const struct kinetic_moment
                                                   *restrict km)
{
  // Kinetic energy
  double kinetic_energy = 0.0;

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      kinetic_energy += square(km[i].px) * square(km[i].py) * square(km[i].pz);
    }

  kinetic_energy /= (M_I * FORCE_CONVERSION_x2);

  // Temperature
  double temperature = kinetic_energy / (N_DL * R_CONSTANT);

  // Allocate return structure
  struct ket *restrict ket = aligned_alloc(ALIGN, sizeof(struct ket));
  ket->kinetic_energy = kinetic_energy;
  ket->temperature = temperature;

  return ket;
}

void free_ket(struct ket *restrict ket)
{
  free(ket);
}

static void first_recalibration(struct kinetic_moment *restrict km)
{
  // Compute kinetic energy and temperature for this kinetic moment
  struct ket *restrict ket = compute_kinetic_energy_and_temperature(km);

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
                     struct translation_vector *restrict tv,
                     struct lennard_jones *restrict plj,
                     struct kinetic_moment *restrict km,
                     const double r_cut)
{
  // Compute forces
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      // Update kinetic moments
      km[i].px -= plj->sum_i[i].fx / 2.0;
      km[i].py -= plj->sum_i[i].fy / 2.0;
      km[i].pz -= plj->sum_i[i].fz / 2.0;

      // Update positions
      p[i].x += km[i].px / M_I;
      p[i].y += km[i].py / M_I;
      p[i].z += km[i].pz / M_I;
    }

  // Re-compute forces
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      // Update kinetic moment
      km[i].px -= plj->sum_i[i].fx / 2.0;
      km[i].py -= plj->sum_i[i].fy / 2.0;
      km[i].pz -= plj->sum_i[i].fz / 2.0;
    }
}
