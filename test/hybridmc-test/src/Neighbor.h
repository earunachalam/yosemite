#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <memory>
#include <vector>

class Particle;

class Neighbor
{
	public:
		std::shared_ptr<Particle> 	m_ptr;
		double 						m_dist = - 1.00;
		std::vector<double>			m_tnvec;

		Neighbor(const std::shared_ptr<Particle>& particle_ptr, const double dist, const std::vector<double>& tnvec):
			m_ptr(particle_ptr), m_dist(dist), m_tnvec(tnvec) {};
};


#endif
