#include <stdio.h>
#include <stdlib.h>

#include "helper.h"
#include "common.h"

struct particle *get_particles(char *restrict filename)
{
  // Begin exploration of the file
  FILE *restrict  stream = fopen(filename, "r");

  if (!stream)
    {
      printf("Error when open the file %s\n", filename);
      exit(ERR_OPEN);
    }

  // Count number of line
  uint64_t count = 0;

  // Skip data
  uint64_t first   = 0;
  uint64_t second  = 0;
  uint64_t useless = 0;
  double one       = 0.0;
  double two       = 0.0;
  double three     = 0.0;

  // Skip the comment line
  fscanf(stream, "%lu %lu\n", &first, &second);

  while (fscanf(stream, "%lu %lf %lf %lf\n", &useless, &one, &two, &three) != EOF)
    {
      count++;
    }

  // Set the number of particles
  N_PARTICLES_TOTAL = count;

  if (LOCAL_EQUAL_TOTAL)
    N_PARTICLES_LOCAL = N_PARTICLES_TOTAL;

  // Close file to re-open it
  fclose(stream);

  // Init particles
  struct particle *restrict p =
    aligned_alloc(ALIGN, sizeof(struct particle) * N_PARTICLES_LOCAL);

  // Begin exploration of the file
  FILE *restrict f = fopen(filename, "r");

  if (!f)
    {
      printf("Error when open the file %s\n", filename);
      exit(ERR_OPEN);
    }

  // Skip the comment line
  fscanf(f, "%lu %lu\n", &first, &second);

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      fscanf(f, "%lu %lf %lf %lf\n", &useless, &(p[i].x), &(p[i].y), &(p[i].z));
    }

  fclose(f);

  return p;
}

void free_particles(struct particle *restrict p)
{
  free(p);
}

void print_particles(const struct particle *restrict p)
{
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      printf("%13e %13e %13e\n", p[i].x, p[i].y, p[i].z);
    }
}

void print_energy(const struct lennard_jones *restrict lj)
{
  printf("energy: %lf\n", lj->energy);
}

double compute_square_distance_3D(const struct particle *restrict a,
                                  const struct particle *restrict b)
{
  const double x_2 = square(a->x - b->x);
  const double y_2 = square(a->y - b->y);
  const double z_2 = square(a->z - b->z);

  return x_2 + y_2 + z_2;
}

uint64_t check_forces(const struct force **restrict f,
                      const double tolerance)
{
  uint64_t error = 0;

  struct force sum =
    {
      .fx = 0.0,
      .fy = 0.0,
      .fz = 0.0
    };

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      // Init sum of forces on particle i
      struct force sum_i =
        {
          .fx = 0.0,
          .fy = 0.0,
          .fz = 0.0
        };

      // Sum on particle i
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          sum_i.fx += f[i][j].fx;
          sum_i.fy += f[i][j].fy;
          sum_i.fz += f[i][j].fz;
        }

      // Update global sum
      sum.fx += sum_i.fx;
      sum.fy += sum_i.fy;
      sum.fz += sum_i.fz;
    }

  // Get absolute value of sum to compare with a tolerance
  struct force abs_sum =
    {
      .fx = abs_double(sum.fx),
      .fy = abs_double(sum.fy),
      .fz = abs_double(sum.fz)
    };

  // Compare with a tolerance
  if (abs_sum.fx > tolerance
      || abs_sum.fy > tolerance
      || abs_sum.fz > tolerance)
    {
      fprintf(stderr, "==error== %s at line %d: "
              "sum of forces apply on particles are NOT null:\n"
              "            -> fx: %13e, fy: %13e,  fz: %13e\n",
              __func__, __LINE__, sum.fx, sum.fy, sum.fz);
      error = 1;
    }
  else
    {
      printf("sum of forces apply on particles are null:\n"
             "  -> fx: %13e, fy: %13e, fz: %13e\n",
             sum.fx, sum.fy, sum.fz);
    }

  return error;
}

struct translation_vector *init_translation_vectors(const uint64_t n)
{
  //
  struct translation_vector *restrict tv =
    aligned_alloc(ALIGN, sizeof(struct translation_vector) * n);

  for (uint64_t i = 0; i < n; i++)
    {
      tv[i].x = (double)((int64_t)(i / 9)       - (int64_t)1) * L;
      tv[i].y = (double)((int64_t)((i / 3) % 3) - (int64_t)1) * L;
      tv[i].z = (double)((int64_t)(i % 3)       - (int64_t)1) * L;
    }

  return tv;
}

void print_translation_vectors(const struct translation_vector *restrict tv,
                               const uint64_t n)
{
  for (uint64_t i = 0; i < n; i++)
    printf("x: %13e, y: %13e, z: %13e\n", tv[i].x, tv[i].y, tv[i].z);
}

void free_translation_vector(struct translation_vector *restrict tv)
{
  free(tv);
}
