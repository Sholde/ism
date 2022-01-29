#ifndef _IO_H_
#define _IO_H_

/**
 * reset_file - Reset file
 */
void reset_file(const char *filename);

/**
 * store_particles - Store particles in file nammed filename in PDB format
 * @param filename: file name
 * @param p       : sturct that contain position of particles
 * @param ite     : iteration number
 * @return
 */
void store_particles(const char *filename, const struct particle *restrict p,
                     const double temperature, const uint64_t ite);

#endif // _IO_H_
