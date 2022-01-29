#ifndef _VELOCITY_VERLET_H_
#define _VELOCITY_VERLET_H_

//
struct ket *init_ket(void);

//
void compute_kinetic_energy_and_temperature(struct ket *restrict ket,
                                            const struct kinetic_moment
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

//
void berendsen_thermostat(struct kinetic_moment *restrict km,
                          struct ket *restrict ket);

#endif // _VELOCITY_VERLET_H_
