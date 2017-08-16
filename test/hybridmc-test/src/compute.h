#ifndef COMPUTE_H
#define COMPUTE_H

#include <memory>
#include <vector>

#include "Particle.h"

struct computed_particlewise_vecdouble
{
	unsigned int m_id = -1;
	std::vector<double> m_value;
};

struct computed_particlewise_double
{
	unsigned int m_id = -1;
	double m_value = 0.00;
};

void compute_neighbors_orthorhombic(std::vector<std::shared_ptr<Particle> >& particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, unsigned int nneigh_track);
void compute_q6(std::vector<std::shared_ptr<Particle> >& particles, std::vector<computed_particlewise_vecdouble>& q6_all);

#endif
