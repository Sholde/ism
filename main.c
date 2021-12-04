#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define R_STAR       3.0
#define EPSILON_STAR 0.2

#define ALIGN 64

#define abs_double(x) (x < 0.0 ? -x : x)

// Maths macros
#define square(x) ((x) * (x))
#define cube(x)   ((x) * (x) * (x))
#define quad(x)   ((x) * (x) * (x) * (x))
#define penta(x)  ((x) * (x) * (x) * (x) * (x))
#define hexa(x)   ((x) * (x) * (x) * (x) * (x) * (x))
#define septa(x)  ((x) * (x) * (x) * (x) * (x) * (x) * (x))

uint64_t N_PARTICLES_TOTAL = 0;
uint64_t N_PARTICLES_LOCAL = 0;
uint64_t LOCAL_EQUAL_TOTAL = 1;

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
      printf("%e %e %e\n", p[i].x, p[i].y, p[i].z);
    }
}

const double compute_square_distance_3D(const struct particle *restrict a,
                                        const struct particle *restrict b)
{
  const double x_2 = square(a->x - b->x);
  const double y_2 = square(a->y - b->y);
  const double z_2 = square(a->z - b->z);

  return x_2 + y_2 + z_2;
}

struct lennard_jones *lennard_jones(const struct particle *restrict p)
{
  // Init lennard_jones
  struct lennard_jones *restrict lj =
    aligned_alloc(ALIGN, sizeof(struct lennard_jones));

  // Init forces
  lj->f = aligned_alloc(ALIGN, sizeof(struct force *) * N_PARTICLES_LOCAL);

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    lj->f[i] = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_LOCAL);

  //
  lj->energy = 0.0;

  // Compute forces
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      for (uint64_t j = i + 1; j < N_PARTICLES_LOCAL; j++)
        {
          const double distance = compute_square_distance_3D(p + i, p + j);

          const double R_STAR_distance = square(R_STAR) / distance;

          const double u_ij =
            EPSILON_STAR * (hexa(R_STAR_distance) - 2.0 * cube(R_STAR_distance));

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

          // Update energy
          lj->energy += u_ij;
        }

      // Update force on particle i with i
      lj->f[i][i].fx = 0.0;
      lj->f[i][i].fy = 0.0;
      lj->f[i][i].fz = 0.0;
    }

  lj->energy *= 4.0;

  return lj;
}

void free_lennard_jones(struct lennard_jones *restrict lj)
{
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    free(lj->f[i]);

  free(lj->f);
  free(lj);
}

void print_energy(const struct lennard_jones *restrict lj)
{
  printf("energy: %lf\n", lj->energy);
}

uint64_t check_forces(const struct force **restrict f)
{
  uint64_t error = 0;
  uint64_t count = 0;

  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      // Init sum of forces on particle i
      struct force sum_i =
        {
          .fx = 0.0,
          .fy = 0.0,
          .fz = 0.0
        };

      // Sum
      for (uint64_t j = 0; j < N_PARTICLES_LOCAL; j++)
        {
          sum_i.fx += f[i][j].fx;
          sum_i.fy += f[i][j].fy;
          sum_i.fz += f[i][j].fz;
        }

      // Get absolute value of sum to compare with a tolerance
      struct force abs_sum_i =
        {
          .fx = abs_double(sum_i.fx),
          .fy = abs_double(sum_i.fy),
          .fz = abs_double(sum_i.fz)
        };

      // Compare with a tolerance
      double tolerance = 1.0e-8;

      if (abs_sum_i.fx > tolerance
          || abs_sum_i.fy > tolerance
          || abs_sum_i.fz > tolerance)
        {
#ifdef DEBUG
          fprintf(stderr, "==error== %s at line %d: "
                  "sum of forces apply on particle %lu are not null\n"
                  "            -> x y z: %lf %lf %lf\n\n",
                  __func__, __LINE__, i, sum_i.fx, sum_i.fy, sum_i.fz);
#endif
          count++;
          error = 1;
        }
    }

  printf("number of non-zero force sums: %lu\n", count);

  return error;
}

int main(int argc, char **argv)
{
  // Check argument
  if (argc < 2)
    {
      printf("Usage: %s [FILENAME] [NB PARTICLES LOCAL - OPTIONAL]\n", argv[0]);
      return ERR_USAGE;
    }

  // Update argument
  if (argc == 3)
    {
      N_PARTICLES_LOCAL = atol(argv[2]);
      LOCAL_EQUAL_TOTAL = 0;
    }

  // Particles
  struct particle *restrict p = get_particles(argv[1]);
  //print_particles(p);

  // Force
  struct lennard_jones *restrict lj = lennard_jones(p);
  print_energy(lj);
  uint64_t error __attribute__((unused)) = check_forces(lj->f);

  // Release memory
  free_lennard_jones(lj);
  free_particles(p);

  return 0;
}
