#include "compute.h"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <omp.h>

#include "Neighbor.h"
#include "Particle.h"
#include "real_spherical_harmonics.h"

#define THRESHOLD_q6 3.7

// first 13 elements of vector to contain q6m values and the final element to contain the q6 value
void compute_q6(std::vector<std::shared_ptr<Particle> >& particles)
{
	for (unsigned int ip = 0; ip < particles.size(); ++ip)
	{
		if (particles.at(ip)->type == WATERTYPE)
		{
			// get number of neighbors which are less than THRESHOLD_q6 away from current particle
			// note that this relies on pre-sorting of neighbor list
			unsigned int currNNeigh = 0;
			for (auto& jNeigh : particles.at(ip)->neighbor_ptrs)
				if ((jNeigh->dist < THRESHOLD_q6) && (jNeigh->selfptr->type == WATERTYPE))
					++currNNeigh;
			
			if (currNNeigh > 0)
			{
				// reset to zero
				particles.at(ip)->q6m.zeros();

				double weight_ij = 1.0;
				double totweight_i = static_cast<double>(currNNeigh);

				for (unsigned int j_neigh_idx = 0; j_neigh_idx < currNNeigh; ++j_neigh_idx)
				{
					if (ip != j_neigh_idx)
					{
						double costheta = particles.at(ip)->neighbor_ptrs.at(j_neigh_idx)->tnvec.at(2)/particles.at(ip)->neighbor_ptrs.at(j_neigh_idx)->dist;
						double phi = atan2(particles.at(ip)->neighbor_ptrs.at(j_neigh_idx)->tnvec.at(1), particles.at(ip)->neighbor_ptrs.at(j_neigh_idx)->tnvec.at(0));

						particles.at(ip)->q6m(0) += S6neg6(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(1) += S6neg5(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(2) += S6neg4(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(3) += S6neg3(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(4) += S6neg2(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(5) += S6neg1(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(6) += S6zer0(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(7) += S6pos1(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(8) += S6pos2(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(9) += S6pos3(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(10) += S6pos4(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(11) += S6pos5(costheta,phi)*weight_ij;
						particles.at(ip)->q6m(12) += S6pos6(costheta,phi)*weight_ij;
					}
				}

				// Normalize each q6m component by the number of neighbors, and calculate the statistics over the entire system
				
				double totsq = 0.00;
				for (unsigned int m = 0; m < 13; ++m)
				{
					particles.at(ip)->q6m(m) /= totweight_i;
					totsq += pow(particles.at(ip)->q6m(m), 2.0);
				}

				particles.at(ip)->q6 = sqrt(totsq*4.0*M_PI/13.0);
				
				// check for anomalies
				if (particles.at(ip)->q6m(6) < -0.4)
				{
					printf("ID = %d, q60 = %lf, currNNeigh = %d\n", particles.at(ip)->id, particles.at(ip)->q6m(6), currNNeigh);
					printf("%lf %lf %lf \n", particles.at(ip)->x(0), particles.at(ip)->x(1), particles.at(ip)->x(2));
					for (uint n = 0; n < currNNeigh; ++n)
					{
						printf("dist = %lf\n", particles.at(ip)->neighbor_ptrs.at(n)->dist);
						std::cout << particles.at(ip)->neighbor_ptrs.at(n)->selfptr->x << std::endl;
					}
				}
			}
			else
			{
				printf("No neighbors:  ID = %d\n", particles.at(ip)->id);
			}
		}
	}
}
