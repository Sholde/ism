#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "helper.h"
#include "common.h"
#include "lennard_jones.h"
#include "velocity_verlet.h"

// Global variable
uint64_t N_PARTICLES_TOTAL = 0;
uint64_t N_PARTICLES_LOCAL = 0;
uint64_t LOCAL_EQUAL_TOTAL = 1;
uint64_t N_DL = 0;

// Structure to monitoring simulation
struct timespec simulation_clock;

//
static void handle_argument(const int argc, const char **argv,
                            char *restrict filename)
{
  // Check argument
  if (argc < 2)
    {
      printf("Usage: %s [FILENAME] [NB PARTICLES LOCAL - OPTIONAL]\n", argv[0]);
      exit(ERR_USAGE);
    }

  // Update argument
  if (argc == 3)
    {
      N_PARTICLES_LOCAL = atol(argv[2]);
      LOCAL_EQUAL_TOTAL = 0;
    }

  //
  strcpy(filename, argv[1]);
}

//
static void run_lennard_jones(const char *restrict filename)
{
  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(filename);
  //print_particles(p);

  // Init lennard jones
  struct lennard_jones *restrict lj = init_lennard_jones();

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  // Run lennard jones
  lennard_jones(lj, p);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  // Print
  printf("== Lennard Jones ==\n");
  print_energy(lj);
  uint64_t error __attribute__((unused)) = check_forces(lj->f, TOLERANCE);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Release memory
  free_lennard_jones(lj);
}

//
static void run_periodical_lennard_jones(const char *restrict filename)
{
  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(filename);
  //print_particles(p);

  // Generate translation vectors
  struct translation_vector *restrict tv = init_translation_vectors(N_SYM);
  //print_translation_vectors(tv, N_SYM);

  // Init lennard jones
  struct lennard_jones *restrict plj = init_lennard_jones();

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  // Run periodical lennard jones
  const double r_cut = 10.0;
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  // Print
  printf("== Periodical Lennard Jones ==\n");
  print_energy(plj);
  uint64_t plj_error __attribute__((unused)) = check_forces(plj->f, TOLERANCE);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Release memory
  free_lennard_jones(plj);
  free_translation_vector(tv);
}

//
static void run_velocity_verlet(const char *restrict filename)
{
  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(filename);
  //print_particles(p);

  // Generate translation vectors
  struct translation_vector *restrict tv = init_translation_vectors(N_SYM);

  // Init lennard jones
  struct lennard_jones *restrict plj = init_lennard_jones();

  // r_cut
  const double r_cut = 10.0;

  // Velocity verlet
  printf("\n== Velocity Verlet ==\n");

  struct kinetic_moment *restrict km = init_velocity_verlet();
  periodical_lennard_jones(plj, p, tv, r_cut, N_SYM);

  //
  uint64_t n_step = 10;

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  printf("\n");
  printf("           "
         "%14s "
         "%15s "
         "%16s "
         "%18s "
         "%17s "
         "\n",
         "TEMPERATURE",
         "TOTAL_ENERGY",
         "KINETIC_ENERGY",
         "POTENTIAL_ENERGY",
         "FORCES_SUM_NORM");

  struct ket *restrict ket = compute_kinetic_energy_and_temperature(km);
  printf("STEP %2d -- ", 0);
  printf("%14e ", ket->temperature);
  printf("%15e ", ket->kinetic_energy + plj->energy);
  printf("%16e ", ket->kinetic_energy);
  printf("%18e ", plj->energy);
  printf("%17e ", abs_double(plj->sum->fx + plj->sum->fy + plj->sum->fz));
  printf("\n");
  free_ket(ket);

  for (uint64_t n = 1; n < n_step + 1; n++)
    {
      //
      velocity_verlet(p, tv, plj, km, r_cut);

      //
      ket = compute_kinetic_energy_and_temperature(km);

      printf("STEP %2ld -- ", n);
      printf("%14e ", ket->temperature);
      printf("%15e ", ket->kinetic_energy + plj->energy);
      printf("%16e ", ket->kinetic_energy);
      printf("%18e ", plj->energy);
      printf("%17e ", square(plj->sum->fx) + square(plj->sum->fy) + square(plj->sum->fz));
      printf("\n");

      free_ket(ket);
    }

  printf("\n");

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1e-9;

  // Print
  printf("Simulate: %e seconds\n", (double)n_step * DT);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Release memory
  free_kinetic_moment(km);
  free_lennard_jones(plj);
  free_particles(p);
}

int main(int argc, char **argv)
{
  //
  char filename[256];

  // Handle command line argument
  handle_argument(argc, argv, filename);

  // Run
  run_lennard_jones(filename);
  run_periodical_lennard_jones(filename);
  run_velocity_verlet(filename);

  return 0;
}
