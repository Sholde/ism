#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "helper.h"
#include "common.h"
#include "lennard_jones.h"

// Global variable
uint64_t N_PARTICLES_TOTAL = 0;
uint64_t N_PARTICLES_LOCAL = 0;
uint64_t LOCAL_EQUAL_TOTAL = 1;

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
