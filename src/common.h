#ifndef _COMMON_H_
#define _COMMON_H_

// Particles
struct particle *get_particles(const char *restrict filename);
void free_particles(struct particle *restrict p);
void print_particles(const struct particle *restrict p);
void print_energy(const struct lennard_jones *restrict lj);

// Compute
double compute_square_distance_3D(const struct particle *restrict a,
                                  const struct particle *restrict b);
uint64_t check_forces(const struct force **restrict f,
                      const double tolerance);

// Translation vectors
struct translation_vector *init_translation_vectors(const uint64_t n);
void print_translation_vectors(const struct translation_vector *restrict tv,
                               const uint64_t n);
void free_translation_vector(struct translation_vector *restrict tv);

#endif // _COMMON_H_
