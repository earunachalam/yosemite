#ifndef PARTICLE_H
#define PARTICLE_H

#include <armadillo>

#include "Neighbor.h"

#define WATERTYPE 1

class Particle
{
	public:
		unsigned int id = -1;										// unique identifier for each particle
		unsigned int type = -1;										// 1 for water, 2 for surface
		arma::vec x;												// position
		std::vector<std::shared_ptr<Neighbor> > neighbor_ptrs;		// list of pointers to nearest neighbors

		double zdist; 												// distance to surface (z axis)

		arma::vec q6m = arma::zeros(13);							// values for each of the 2l+1 components
		double q6 = 0;												// prop to sum of squares of each component

		Particle(unsigned int arg_id, unsigned int arg_type, double x, double y, double z): id(arg_id), type(arg_type)
		{
			this->x = {x, y, z};
		};
};

#endif

