#include <stdlib.h>

#include "helper.h"
#include "common.h"
#include "lennard_jones.h"

//
struct lennard_jones *init_lennard_jones(void)
{
  // Allocate memory
  struct lennard_jones *restrict lj =
    aligned_alloc(ALIGN, sizeof(struct lennard_jones));

  lj->f =
    aligned_alloc(ALIGN, sizeof(struct force *restrict) * N_PARTICLES_LOCAL);

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    lj->f[i] = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_LOCAL);

  lj->sum_i = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_LOCAL);
  lj->sum = aligned_alloc(ALIGN, sizeof(struct force));

  // Init energy to 0
  lj->energy = 0.0;

  // Init force
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          lj->f[i][j].fx = 0.0;
          lj->f[i][j].fy = 0.0;
          lj->f[i][j].fz = 0.0;
        }
    }

  // Init sum of force apply on particle i to 0
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      lj->sum_i[i].fx = 0.0;
      lj->sum_i[i].fy = 0.0;
      lj->sum_i[i].fz = 0.0;
    }

  // Init sum of force to 0
  lj->sum->fx = 0.0;
  lj->sum->fy = 0.0;
  lj->sum->fz = 0.0;

  return lj;
}

//
void free_lennard_jones(struct lennard_jones *restrict lj)
{
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    free(lj->f[i]);

  free(lj->sum_i);
  free(lj->sum);
  free(lj->f);
  free(lj);
}

//
void lennard_jones(struct lennard_jones *restrict lj,
                   const struct particle *restrict p)
{
  // Init energy to 0
  lj->energy = 0.0;

  // Init force
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          lj->f[i][j].fx = 0.0;
          lj->f[i][j].fy = 0.0;
          lj->f[i][j].fz = 0.0;
        }
    }

  // Init sum of force apply on particle i to 0
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      lj->sum_i[i].fx = 0.0;
      lj->sum_i[i].fy = 0.0;
      lj->sum_i[i].fz = 0.0;
    }

  // Init sum of force to 0
  lj->sum->fx = 0.0;
  lj->sum->fy = 0.0;
  lj->sum->fz = 0.0;

  // Compute
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      for (uint64_t j = i + 1; j < N_PARTICLES_LOCAL; j++)
        {
          const double distance = compute_square_distance_3D(p + i, p + j);

          const double R_STAR_distance = square(R_STAR) / distance;

          const double u_ij =
            (hexa(R_STAR_distance) - 2.0 * cube(R_STAR_distance));

          // Update energy
          lj->energy += u_ij;

          // Update forces
          const double du_ij =
            -48.0 * EPSILON_STAR * (septa(R_STAR_distance) - quad(R_STAR_distance));

          // Update force on particle i with j
          lj->f[i][j].fx = du_ij * (p[i].x - p[j].x);
          lj->f[i][j].fy = du_ij * (p[i].y - p[j].y);
          lj->f[i][j].fz = du_ij * (p[i].z - p[j].z);

          // Update force on particle j with i
          lj->f[j][i].fx = - lj->f[i][j].fx;
          lj->f[j][i].fy = - lj->f[i][j].fy;
          lj->f[j][i].fz = - lj->f[i][j].fz;

          // Update sum
          lj->sum_i[i].fx += lj->f[i][j].fx;
          lj->sum_i[i].fy += lj->f[i][j].fy;
          lj->sum_i[i].fz += lj->f[i][j].fz;

          lj->sum_i[j].fx += lj->f[j][i].fx;
          lj->sum_i[j].fy += lj->f[j][i].fy;
          lj->sum_i[j].fz += lj->f[j][i].fz;
        }

      // Update force on particle i with i
      lj->f[i][i].fx = 0.0;
      lj->f[i][i].fy = 0.0;
      lj->f[i][i].fz = 0.0;

      // Update sum
      lj->sum->fx += lj->sum_i[i].fx;
      lj->sum->fy += lj->sum_i[i].fy;
      lj->sum->fz += lj->sum_i[i].fz;
    }

  // Update energy
  lj->energy *= 4.0 * EPSILON_STAR;
}

//
void periodical_lennard_jones(struct lennard_jones *restrict plj,
                              const struct particle *restrict p,
                              const struct translation_vector *restrict tv,
                              const double r_cut, const uint64_t n)
{
  // Init energy to 0
  plj->energy = 0.0;

  // Init force
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          plj->f[i][j].fx = 0.0;
          plj->f[i][j].fy = 0.0;
          plj->f[i][j].fz = 0.0;
        }
    }

  // Init sum of force apply on particle i to 0
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      plj->sum_i[i].fx = 0.0;
      plj->sum_i[i].fy = 0.0;
      plj->sum_i[i].fz = 0.0;
    }

  // Init sum of force to 0
  plj->sum->fx = 0.0;
  plj->sum->fy = 0.0;
  plj->sum->fz = 0.0;

  // Compute
  for (uint64_t k = 0; k < n; k++)
    {
      for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
        {
          for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
            {
              // Test if i == j and then ignore this step
              if (i == j)
                continue;

              const struct particle tmp_j =
                {
                  .x = p[j].x + tv[k].x,
                  .y = p[j].y + tv[k].y,
                  .z = p[j].z + tv[k].z
                };

              const double distance = compute_square_distance_3D(p + i, &tmp_j);

              // Test if the distance is under r_cut and then ignore this step
              if (distance > square(r_cut))
                continue;

              const double R_STAR_distance = square(R_STAR) / distance;

              const double u_ij =
                (hexa(R_STAR_distance) - 2.0 * cube(R_STAR_distance));

              // Update energy
              plj->energy += u_ij;

              // Update forces
              const double du_ij =
                -48.0 * EPSILON_STAR * (septa(R_STAR_distance) - quad(R_STAR_distance));

              // Update force on particle i with j
              plj->f[i][j].fx += du_ij * (p[i].x - tmp_j.x);
              plj->f[i][j].fy += du_ij * (p[i].y - tmp_j.y);
              plj->f[i][j].fz += du_ij * (p[i].z - tmp_j.z);
            }
        }
    }

  // Update sum
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      // Update sum_i
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          plj->sum_i[i].fx += plj->f[i][j].fx;
          plj->sum_i[i].fy += plj->f[i][j].fy;
          plj->sum_i[i].fz += plj->f[i][j].fz;
        }

      // Update sum
      plj->sum->fx += plj->sum_i[i].fx;
      plj->sum->fy += plj->sum_i[i].fy;
      plj->sum->fz += plj->sum_i[i].fz;
    }

  // Update energy
  plj->energy *= 2.0 * EPSILON_STAR;
}
