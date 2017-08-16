#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

#include "Neighbor.h"

class Particle
{
	public:
		unsigned int m_id = -1,
					 m_type = -1;
		std::vector<double> m_x;
		std::vector<std::shared_ptr<Neighbor> > m_neighbor_ptrs;

		Particle(unsigned int arg_id, unsigned int arg_type, std::vector<double> arg_x): m_id(arg_id), m_type(arg_type), m_x(arg_x) {};
};

#endif

