#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "helper.h"
#include "io.h"

void reset_file(const char *filename)
{
  FILE *restrict f = fopen(filename, "w");

  if (!f)
    {
      printf("Error when open the file %s\n", filename);
      exit(ERR_OPEN);
    }

  fclose(f);

  return;
}

void store_particles(const char *filename, const struct particle *restrict p,
                     const uint64_t ite)
{
  // Begin exploration of the file
  FILE *restrict f = fopen(filename, "a");

  if (!f)
    {
      printf("Error when open the file %s\n", filename);
      exit(ERR_OPEN);
    }

  fseek(f, 0, SEEK_END);

  // Print first lines
  fprintf(f, "CRYST1  %.2lf  %.2lf  %.2lf  90.00  90.00  90.00  P  1\n", L, L, L);
  fprintf(f, "MODEL  %ld\n", ite);

  // Print positions
  for (uint64_t i = 0; i < N_PARTICLES_LOCAL; i++)
    {
      fprintf(f, "ATOM  %5ld  C   0  %10.3lf  %10.3lf  %10.3lf  MRES\n",
              i + 1, p[i].x, p[i].y, p[i].z);
    }

  // Print last lines
  fprintf(f, "TER\n");
  fprintf(f, "ENDMDL\n");

  // Close file
  fclose(f);

  return;
}
