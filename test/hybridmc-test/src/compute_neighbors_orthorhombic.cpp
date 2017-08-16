#include "compute.h"

#include <cassert>
#include <cmath>
#include <csignal>
#include <limits>

#include "Neighbor.h"
#include "Particle.h"


#define WATERTYPE 1

void compute_neighbors_orthorhombic(std::vector<std::shared_ptr<Particle> >& particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, unsigned int nneigh_track)
{
	double Lx = fabs(xmax - xmin),
		   Ly = fabs(ymax - ymin),
		   Lz = fabs(zmax - zmin);

	for (unsigned int i_particle_idx = 0; i_particle_idx < particles.size(); ++i_particle_idx)
	{
		std::shared_ptr<Particle> i_part_ptr = particles.at(i_particle_idx);
		std::vector<double> i_pos = i_part_ptr->m_x;
		if (i_part_ptr->m_type != WATERTYPE) continue;

		double i_x = i_pos.at(0),
			   i_y = i_pos.at(1),
			   i_z = i_pos.at(2);

		std::vector<std::vector<double> > tnvecs(nneigh_track);
		std::vector<double> distances(nneigh_track, std::numeric_limits<double>::max());
		std::vector<std::shared_ptr<Particle> > ptcl_neigh_ptrs(nneigh_track, nullptr);

		for (unsigned j_particle_idx = 0; j_particle_idx < particles.size(); ++j_particle_idx)
		{
			if (i_particle_idx == j_particle_idx) continue;

			std::shared_ptr<Particle> j_part_ptr = particles.at(j_particle_idx);
			std::vector<double> j_pos = j_part_ptr->m_x;
			if (j_part_ptr->m_type != WATERTYPE) continue;

			double j_x = j_pos.at(0),
				   j_y = j_pos.at(1),
				   j_z = j_pos.at(2);

			double dx = j_x - i_x,
				   dy = j_y - i_y,
				   dz = j_z - i_z;

			dx -= std::round(dx/Lx)*Lx;
			dy -= std::round(dy/Ly)*Ly;
			dz -= std::round(dz/Lz)*Lz;

			double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

			for (unsigned int k_neigh_idx = 0; k_neigh_idx < nneigh_track; ++k_neigh_idx)
			{
				if (distance < distances.at(k_neigh_idx))
				{
					std::vector<double> tnvec = {dx, dy, dz};

					for (unsigned int l_neigh_idx = nneigh_track - 2; l_neigh_idx == k_neigh_idx; l_neigh_idx--)
					{
						tnvecs.at(l_neigh_idx + 1) = tnvec;
						distances.at(l_neigh_idx + 1) = distances.at(l_neigh_idx);
						ptcl_neigh_ptrs.at(l_neigh_idx + 1) = ptcl_neigh_ptrs.at(l_neigh_idx);
					}

					tnvecs.at(k_neigh_idx) = tnvec;
					distances.at(k_neigh_idx) = distance;
					ptcl_neigh_ptrs.at(k_neigh_idx) = j_part_ptr;
					break;
				}
			}
		}
		
		for (unsigned int close_neigh_idx = 0;  close_neigh_idx < nneigh_track; ++close_neigh_idx)
		{
			particles.at(i_particle_idx)->m_neighbor_ptrs.push_back(std::make_shared<Neighbor>(ptcl_neigh_ptrs.at(close_neigh_idx), distances.at(close_neigh_idx), tnvecs.at(close_neigh_idx)));
		}
	}
}
