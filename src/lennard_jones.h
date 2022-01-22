#ifndef _LENNARD_JONES_H_
#define _LENNARD_JONES_H_

//
struct lennard_jones *init_lennard_jones(void);

//
void free_lennard_jones(struct lennard_jones *restrict lj);

//
void lennard_jones(struct lennard_jones *restrict lj,
                   const struct particle *restrict p);

//
void periodical_lennard_jones(struct lennard_jones *restrict plj,
                              const struct particle *restrict p,
                              const struct translation_vector *restrict tv,
                              const double r_cut, const uint64_t n);

#endif // _LENNARD_JONES_H_
