#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define N_PARTICLES_TOTAL 1000
#define N_PARTICLES_LOCAL 500

#define R       3.0
#define EPSILON 0.2

#define ALIGN 64

#define abs_double(x) (x < 0.0 ? -x : x)

#if N_PARTICLES_LOCAL >= N_PARTICLES_TOTAL
#error "N_PARTICLES_LOCAL >= N_PARTICLES_TOTAL"
#endif

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

char *get_field(char *line, uint64_t num)
{
  char* tok;
  for (tok = strtok(line, " ");
       tok && *tok;
       tok = strtok(NULL, " "))
    {
      if (!--num)
        return tok;
    }
  return NULL;
}

struct particle *get_particles(char *filename)
{
  // Init particles
  struct particle *p =
    aligned_alloc(ALIGN, sizeof(struct particle) * N_PARTICLES_TOTAL);

  // Begin exploration of the file
  FILE* stream = fopen(filename, "r");

  if (!stream)
    {
      printf("Error when open the file %s\n", filename);
      exit(ERR_OPEN);
    }

  char line[1024];
  uint64_t count = 0;

  // Skip the comment line
  fgets(line, 1024, stream);

  while (fgets(line, 1024, stream))
    {
      char *tmp = strdup(line);

      char *x = get_field(tmp, 2);
      x[strlen(x) - 1] = '\0';

      tmp = strdup(line);

      char *y = get_field(tmp, 3);
      y[strlen(y) - 1] = '\0';

      tmp = strdup(line);

      char *z = get_field(tmp, 4);
      z[strlen(z) - 1] = '\0';

      p[count].x = strtod(x, NULL);
      p[count].y = strtod(y, NULL);
      p[count].z = strtod(z, NULL);

      free(tmp);
      count++;
    }

  fclose(stream);

  return p;
}

void free_particles(struct particle *p)
{
  free(p);
}

void print_particles(const struct particle *p)
{
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      printf("%e %e %e\n", p[i].x, p[i].y, p[i].z);
    }
}

const double compute_distance_3D(const struct particle *a,
                                 const struct particle *b)
{
  const double x_2 = pow(a->x - b->x, 2);
  const double y_2 = pow(a->y - b->y, 2);
  const double z_2 = pow(a->z - b->z, 2);

  return sqrt(x_2 + y_2 + z_2);
}

struct force **get_forces(const struct particle *p)
{
  // Init forces
  struct force **f =
    aligned_alloc(ALIGN, sizeof(struct force *) * N_PARTICLES_TOTAL);

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    f[i] = aligned_alloc(ALIGN, sizeof(struct force) * N_PARTICLES_TOTAL);

  double energy = 0.0;

  // Compute forces
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      for (uint64_t j = i + 1; j < N_PARTICLES_TOTAL; j++)
        {
          double distance = compute_distance_3D(p + i, p + j);
          double R_distance = R / distance;

          double u_ij =
            -48.0 * EPSILON * (pow(R_distance, 12) - 2.0 * pow(R_distance, 6));

          // Update force on particle i and j
          f[i][j].fx = u_ij * (p[i].x - p[j].x);
          f[i][j].fy = u_ij * (p[i].y - p[j].y);
          f[i][j].fz = u_ij * (p[i].z - p[j].z);

          f[j][i].fx = f[i][j].fx;
          f[j][i].fy = f[i][j].fy;
          f[j][i].fz = f[i][j].fz;

          // Update energy
          energy += f[i][j].fx;
          energy += f[i][j].fy;
          energy += f[i][j].fz;
        }

      // Update force on particle i and i
      f[i][i].fx = 0.0;
      f[i][i].fy = 0.0;
      f[i][i].fz = 0.0;
    }

  energy *= 4.0;
  printf("energy: %lf\n", energy);

  return f;
}

void free_forces(struct force **f)
{
  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      free(f[i]);
    }

  free(f);
}

uint64_t check_forces(const struct force **f)
{
  uint64_t error = 0;

  for (uint64_t i = 0; i < N_PARTICLES_TOTAL; i++)
    {
      // Init sum of forces on particle i
      struct force sum_i =
        {
          .fx = 0.0,
          .fy = 0.0,
          .fz = 0.0
        };

      // Sum
      for (uint64_t j = 0; j < N_PARTICLES_TOTAL; j++)
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
      double tolerance = 1.0e-3;

      if (abs_sum_i.fx > tolerance
          || abs_sum_i.fy > tolerance
          || abs_sum_i.fz > tolerance)
        {
          fprintf(stderr, "==error== %s at line %d: "
                  "sum of forces apply on particle %lu are not null\n"
                  "            -> x y z: %lf %lf %lf\n\n",
                  __func__, __LINE__, i, sum_i.fx, sum_i.fy, sum_i.fz);

          error = 1;
        }
    }

  return error;
}

int main(int argc, char **argv)
{
  // Check argument
  if (argc != 2)
    {
      printf("Usage: %s [FILENAME]\n", argv[0]);
      return ERR_USAGE;
    }

  // Particles
  struct particle *p = get_particles(argv[1]);
  print_particles(p);

  // Force
  struct force **f = get_forces(p);
  uint64_t error = check_forces(f);

  // Release memory
  free_forces(f);
  free_particles(p);

  return 0;
}
