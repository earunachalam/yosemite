#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <armadillo>
#include <memory>
#include <vector>

class Particle;

class Neighbor
{
	public:
		std::shared_ptr<Particle> 	selfptr;
		double 						dist = - 1.00;
		arma::vec					tnvec;

		Neighbor(const std::shared_ptr<Particle>& particle_ptr, const double dist, const arma::vec& tnvec):
			selfptr(particle_ptr), dist(dist), tnvec(tnvec) {};
};


#endif
