#ifndef COMPUTE_H
#define COMPUTE_H

#include <memory>
#include <vector>

#include "Particle.h"

void compute_neighbors_orthorhombic(std::vector<std::shared_ptr<Particle> >& particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, unsigned int nneigh_track);
void compute_q6(std::vector<std::shared_ptr<Particle> >& particles);
void compute_distance_surface(std::vector<std::shared_ptr<Particle> >& particles, unsigned int one_surf_idx, double zmin, double zmax);

#endif
