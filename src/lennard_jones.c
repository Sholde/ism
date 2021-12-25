#include <stdlib.h>

#include "helper.h"
#include "common.h"
#include "lennard_jones.h"

struct lennard_jones *lennard_jones(const struct particle *restrict p)
{
  // Init lennard_jones
  struct lennard_jones *restrict lj =
    aligned_alloc(ALIGN, sizeof(struct lennard_jones));

  // Init forces
  lj->f = aligned_alloc(ALIGN, sizeof(struct force *restrict) * N_PARTICLES_LOCAL);

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    lj->f[i] = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_LOCAL);

  // Init energy to 0
  lj->energy = 0.0;

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
        }

      // Update force on particle i with i
      lj->f[i][i].fx = 0.0;
      lj->f[i][i].fy = 0.0;
      lj->f[i][i].fz = 0.0;
    }

  //
  lj->energy *= 4.0 * EPSILON_STAR;

  return lj;
}

void free_lennard_jones(struct lennard_jones *restrict lj)
{
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    free(lj->f[i]);

  free(lj->f);
  free(lj);
}

struct lennard_jones *periodical_lennard_jones(const struct particle
                                               *restrict p,
                                               const struct translation_vector
                                               *restrict tv,
                                               const double r_cut,
                                               const uint64_t n)
{
  // Init lennard_jones
  struct lennard_jones *restrict plj =
    aligned_alloc(ALIGN, sizeof(struct lennard_jones));

  // Init forces
  plj->f = aligned_alloc(ALIGN, sizeof(struct force *restrict) * N_PARTICLES_LOCAL);

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      plj->f[i] = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_LOCAL);

      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          plj->f[i][j].fx = 0.0;
          plj->f[i][j].fy = 0.0;
          plj->f[i][j].fz = 0.0;
        }
    }

  //
  plj->energy = 0.0;

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

  //
  plj->energy *= 2.0 * EPSILON_STAR;

  return plj;
}
