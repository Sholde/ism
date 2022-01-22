#ifndef _VELOCITY_VERLET_H_
#define _VELOCITY_VERLET_H_

//
struct ket *compute_kinetic_energy_and_temperature(const struct kinetic_moment
                                                   *restrict km);
void free_ket(struct ket *restrict ket);

//
struct kinetic_moment *init_velocity_verlet(void);
void free_kinetic_moment(struct kinetic_moment *restrict km);

//
void velocity_verlet(struct particle *restrict p,
                     struct translation_vector *restrict tv,
                     struct lennard_jones *restrict plj,
                     struct kinetic_moment *restrict km,
                     const double r_cut);

#endif // _VELOCITY_VERLET_H_
