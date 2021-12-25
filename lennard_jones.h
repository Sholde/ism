#ifndef _LENNARD_JONES_H_
#define _LENNARD_JONES_H_

struct lennard_jones *lennard_jones(const struct particle *restrict p);
void free_lennard_jones(struct lennard_jones *restrict lj);
struct lennard_jones *periodical_lennard_jones(const struct particle
                                               *restrict p,
                                               const struct translation_vector
                                               *restrict tv,
                                               const double r_cut,
                                               const uint64_t n);

#endif // _LENNARD_JONES_H_
