#include "compute.h"

#include <armadillo>
#include <cassert>
#include <cmath>
#include <csignal>
#include <limits>
#include <omp.h>

#include "Neighbor.h"
#include "Particle.h"


#define WATERTYPE 1

void compute_neighbors_orthorhombic(std::vector<std::shared_ptr<Particle> >& particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, unsigned int nneigh_track)
{
	unsigned int ntotal = particles.size();
	std::vector<std::vector<std::shared_ptr<Neighbor> > > all_neighbors(ntotal);
	arma::vec L = {fabs(xmax - xmin), fabs(ymax - ymin), fabs(zmax - zmin)};

	for (unsigned int ip = 0; ip < particles.size(); ++ip)
	{
		all_neighbors.at(ip).reserve(ntotal);

		std::shared_ptr<Particle> iptr = particles.at(ip);
		if (iptr->type != WATERTYPE) continue;
		arma::vec ipos = iptr->x;

		for (unsigned int jp = 0; jp < ip; ++jp)
		{
			std::shared_ptr<Particle> jptr = particles.at(jp);
			if (jptr->type != WATERTYPE) continue;
			arma::vec jpos = jptr->x;

			arma::vec dx = jpos - ipos;
			arma::vec tnvec = dx - arma::round(dx/L) % L;

			double distance = arma::norm(tnvec);

			all_neighbors.at(ip).push_back(std::make_shared<Neighbor>(jptr, distance, tnvec));
			all_neighbors.at(jp).push_back(std::make_shared<Neighbor>(iptr, distance, -1.00*tnvec));
		}
	}

	for (unsigned int ip = 0; ip < particles.size(); ++ip)
	{
		if (particles.at(ip)->type != WATERTYPE) continue;
		std::sort(all_neighbors.at(ip).begin(), all_neighbors.at(ip).end(), [](std::shared_ptr<Neighbor> n1, std::shared_ptr<Neighbor> n2){return n1->dist < n2->dist;});
		particles.at(ip)->neighbor_ptrs = std::vector<std::shared_ptr<Neighbor> >(all_neighbors.at(ip).begin(), all_neighbors.at(ip).begin()+nneigh_track);
	}
}
