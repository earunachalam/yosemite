#include "compute.h"

#include <cmath>
#include <iostream>

#include "Neighbor.h"
#include "Particle.h"
#include "real_spherical_harmonics.h"

#define WATERTYPE 1
#define THRESHOLD_q6 3.5


// first 13 elements of vector to contain q6m values and the final element to contain the q6 value
void compute_q6(std::vector<std::shared_ptr<Particle> >& particles, std::vector<computed_particlewise_vecdouble>& q6_all_allparticles)
{
	for (unsigned int i_particle_idx = 0; i_particle_idx < particles.size(); ++i_particle_idx)
	{
		if (particles.at(i_particle_idx)->m_type == WATERTYPE)
		{
			// get number of neighbors which are less than THRESHOLD_q6 away from current particle
			// note that this relies on pre-sorting of neighbor list
			unsigned int currNNeigh = 0;
			for (auto& jNeigh : particles.at(i_particle_idx)->m_neighbor_ptrs)
				if ((jNeigh->m_dist < THRESHOLD_q6) && (jNeigh->m_ptr->m_type == WATERTYPE))
					++currNNeigh;
			
			if (currNNeigh > 0)
			{
				// reset to zero
				q6_all_allparticles.at(i_particle_idx).m_value.resize(14, 0.0);

				double weight_ij = 1.0;
				double totweight_i = static_cast<double>(currNNeigh);

				for (unsigned int j_neigh_idx = 0; j_neigh_idx < currNNeigh; ++j_neigh_idx)
				{
					if (i_particle_idx != j_neigh_idx)
					{
						double costheta = particles.at(i_particle_idx)->m_neighbor_ptrs.at(j_neigh_idx)->m_tnvec.at(2)/particles.at(i_particle_idx)->m_neighbor_ptrs.at(j_neigh_idx)->m_dist;
						double phi = atan2(particles.at(i_particle_idx)->m_neighbor_ptrs.at(j_neigh_idx)->m_tnvec.at(1), particles.at(i_particle_idx)->m_neighbor_ptrs.at(j_neigh_idx)->m_tnvec.at(0));

						q6_all_allparticles.at(i_particle_idx).m_value.at(0) += S6neg6(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(1) += S6neg5(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(2) += S6neg4(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(3) += S6neg3(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(4) += S6neg2(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(5) += S6neg1(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(6) += S6zer0(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(7) += S6pos1(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(8) += S6pos2(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(9) += S6pos3(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(10) += S6pos4(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(11) += S6pos5(costheta,phi)*weight_ij;
						q6_all_allparticles.at(i_particle_idx).m_value.at(12) += S6pos6(costheta,phi)*weight_ij;
					}
				}

				// Normalize each q6m component by the number of neighbors, and calculate the statistics over the entire system
				
				double totsq = 0.00;
				for (unsigned int m = 0; m < 13; ++m)
				{
					q6_all_allparticles.at(i_particle_idx).m_value.at(m) /= totweight_i;
					totsq += pow(q6_all_allparticles.at(i_particle_idx).m_value.at(m), 2.0);
				}

				q6_all_allparticles.at(i_particle_idx).m_value.back() = sqrt(totsq*4.0*M_PI/13.0);
				
				// check for anomalies
				if (q6_all_allparticles.at(i_particle_idx).m_value.at(6) < -0.4)
				{
					printf("ID = %d, q60 = %lf, currNNeigh = %d\n", particles.at(i_particle_idx)->m_id, q6_all_allparticles.at(i_particle_idx).m_value.at(6), currNNeigh);
				}
			}
			else
			{
				printf("No neighbors:  ID = %d\n", particles.at(i_particle_idx)->m_id);
			}
		}
	}
}
