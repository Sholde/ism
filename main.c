#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

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

struct translation_vector
{
  double x;
  double y;
  double z;
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
      printf("%13e %13e %13e\n", p[i].x, p[i].y, p[i].z);
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

void print_energy(const struct lennard_jones *restrict lj)
{
  printf("energy: %lf\n", lj->energy);
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

  // Structure to monitoring
  struct timespec clock;
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(argv[1]);
  //print_particles(p);

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &clock);
  before = clock.tv_sec + clock.tv_nsec * 1e-9;

  // Force and Energy
  struct lennard_jones *restrict lj = lennard_jones(p);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &clock);
  after = clock.tv_sec + clock.tv_nsec * 1e-9;

  printf("Lennard Jones:\n");
  print_energy(lj);
  uint64_t error __attribute__((unused)) = check_forces(lj->f, TOLERANCE);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Generate translation vectors
  struct translation_vector *restrict tv = init_translation_vectors(N_SYM);
  //print_translation_vectors(tv, N_SYM);

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &clock);
  before = clock.tv_sec + clock.tv_nsec * 1e-9;

  // Force and Energy
  const double r_cut = 10.0;
  struct lennard_jones *restrict plj = periodical_lennard_jones(p, tv, r_cut, N_SYM);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &clock);
  after = clock.tv_sec + clock.tv_nsec * 1e-9;

  printf("Periodical Lennard Jones:\n");
  print_energy(plj);
  uint64_t plj_error __attribute__((unused)) = check_forces(plj->f, TOLERANCE);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Release memory
  free_lennard_jones(plj);
  free_translation_vector(tv);
  free_lennard_jones(lj);
  free_particles(p);

  return 0;
}
