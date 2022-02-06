#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "helper.h"
#include "common.h"
#include "lennard_jones.h"
#include "velocity_verlet.h"
#include "io.h"
#include "arguments.h"

// Global variable
uint64_t N_PARTICLES_TOTAL = 0;
uint64_t N_PARTICLES_LOCAL = 0;
uint64_t LOCAL_EQUAL_TOTAL = 1;
uint64_t N_DL = 0;
uint64_t N_STEP = 10000;
uint64_t M_STEP = 100;
double R_CUT = 10.0;

const char *const VERSION = "1.0.0";
char INPUT_FILE[256] = "";
char OUTPUT_FILE[256] = "output.pdb";

// Structure to monitoring simulation
struct timespec simulation_clock;

int print_usage(const char *const software)
{
  printf("[Usage] %s [Options]\n", software);
  return EXIT_SUCCESS;
}

int default_action(const char *const arg)
{
  printf("Unrecognized argument: %s\n", arg);
  return EXIT_FAILURE;
}

int print_version(void)
{
  printf("[Version] %s\n", VERSION);
  return EXIT_SUCCESS;
}

int select_input(const char *const arg)
{
  //
  const char *ptr = strchr(arg, '=');
  const char *value = ++ptr;
  strcpy(INPUT_FILE, value);
  return EXIT_SUCCESS;
}

int select_output(const char *const arg)
{
  //
  const char *ptr = strchr(arg, '=');
  const char *value = ++ptr;
  strcpy(OUTPUT_FILE, value);
  return EXIT_SUCCESS;
}

int select_n_step(const char *const arg)
{
  //
  const char *ptr = strchr(arg, '=');
  const uint64_t value = atoll(++ptr);
  N_STEP = value;
  return EXIT_SUCCESS;
}

int select_r_cut(const char *const arg)
{
  //
  const char *ptr = strchr(arg, '=');
  const uint64_t value = atoll(++ptr);
  R_CUT = value;
  return EXIT_SUCCESS;
}

//
static void handle_argument(const int argc, const char **argv)
{
  // Check argument
  if (argc < 2)
    {
      print_usage(argv[0]);
      printf("Use --help or -h option for more details.\n");
      exit(ERR_USAGE);
    }

  //
  initArguments(print_usage, default_action);

  //
  addArgument("--version", "-v", print_version, "Display the software version.");
  addArgument("--input=", "-i=", select_input, "Select input file.");
  addArgument("--output=", "-o=", select_output, "Select output file.");
  addArgument("--nstep=", NULL, select_n_step, "Select N_STEP value.");
  addArgument("--rcut=", NULL, select_r_cut, "Select R_CUT value.");

  //
  parseArguments(argc, argv);

  if (strcmp(INPUT_FILE, "") == 0)
    exit(EXIT_SUCCESS);
}

//
static void print_column_name(void)
{
  printf("              "
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
}

//
static void print_step(const uint64_t step,
                       const double temperature,
                       const double tot_energy,
                       const double kinetic_energy,
                       const double potential_energy,
                       const double norm_sum_forces)
{
#if DEBUG
  printf("STEP %5ld -- ", step);
  printf("%14e ", temperature);
  printf("%15e ", tot_energy);
  printf("%16e ", kinetic_energy);
  printf("%18e ", potential_energy);
  printf("%17e ", norm_sum_forces);
  printf("\n");
#endif
}

//
static void run_lennard_jones(void)
{
  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(INPUT_FILE);
  //print_particles(p);

  // Init lennard jones
  struct lennard_jones *restrict lj = init_lennard_jones();

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

  // Run lennard jones
  lennard_jones(lj, p);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

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
static void run_periodical_lennard_jones(void)
{
  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(INPUT_FILE);
  //print_particles(p);

  // Generate translation vectors
  struct translation_vector *restrict tv = init_translation_vectors(N_SYM);
  //print_translation_vectors(tv, N_SYM);

  // Init lennard jones
  struct lennard_jones *restrict plj = init_lennard_jones();

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

  // Run periodical lennard jones
  periodical_lennard_jones(plj, p, tv, R_CUT, N_SYM);

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

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
static void run_velocity_verlet(void)
{
  //
  reset_file(OUTPUT_FILE);

  //
  double before;
  double after;

  // Particles
  struct particle *restrict p = get_particles(INPUT_FILE);
  //print_particles(p);

  // Generate translation vectors
  struct translation_vector *restrict tv = init_translation_vectors(N_SYM);

  // Init lennard jones
  struct lennard_jones *restrict plj = init_lennard_jones();

  // Velocity verlet
  printf("\n== Velocity Verlet ==\n");

  struct kinetic_moment *restrict km = init_velocity_verlet();
  periodical_lennard_jones(plj, p, tv, R_CUT, N_SYM);

  //
  print_column_name();

  //
  struct ket *restrict ket = init_ket();

  // Step 0
  compute_kinetic_energy_and_temperature(ket, km);
  print_step(0, ket->temperature, ket->kinetic_energy + plj->energy,
             ket->kinetic_energy, plj->energy,
             norm_3d(plj->sum->fx, plj->sum->fz, plj->sum->fz));
  store_particles(OUTPUT_FILE, p, 0);

  // Take time before
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  before = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

  // Launch velocity verlet
  for (uint64_t step = 1; step < N_STEP + 1; step++)
    {
      //
      velocity_verlet(p, tv, plj, km, R_CUT);

      //
      compute_kinetic_energy_and_temperature(ket, km);

      print_step(step, ket->temperature, ket->kinetic_energy + plj->energy,
                 ket->kinetic_energy, plj->energy,
                 norm_3d(plj->sum->fx, plj->sum->fz, plj->sum->fz));

      //
      store_particles(OUTPUT_FILE, p, step);

      //
      if (step % M_STEP == 0)
        berendsen_thermostat(km, ket);
    }

  // Take time after
  clock_gettime(CLOCK_MONOTONIC, &simulation_clock);
  after = simulation_clock.tv_sec + simulation_clock.tv_nsec * 1.0e-9;

  // Print
  printf("\n");
  printf("Simulate: %lf fento-seconds\n", (double)N_STEP * DT);
  printf("Take: %lf seconds\n", after - before);
  printf("\n");

  // Release memory
  free_ket(ket);
  free_kinetic_moment(km);
  free_lennard_jones(plj);
  free_particles(p);
}

int main(int argc, char **argv)
{
  // Handle command line argument
  handle_argument(argc, argv);

  // Run
  run_lennard_jones();
  run_periodical_lennard_jones();
  run_velocity_verlet();

  return 0;
}
