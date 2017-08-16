#include "compute.h"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <omp.h>

#include "Neighbor.h"
#include "Particle.h"
#include "real_spherical_harmonics.h"

#define GC(A) {int status; char * demangled = abi::__cxa_demangle(typeid(A).name(),0,0,&status); std::cout << __LINE__ << ": " #A << "\t" << demangled <<"\n";}                                                                                           
#define P2S(a) std::cout << __LINE__ << ": " << #a << ": " << (a) << std::endl

// calculate z distance of particle from surface
void compute_distance_surface(std::vector<std::shared_ptr<Particle> >& particles, unsigned int one_surf_idx, double zmin, double zmax)
{
	double Lz = zmax - zmin;

	for (unsigned int ip = 0; ip < particles.size(); ++ip)
	{
		if (particles.at(ip)->type == WATERTYPE)
		{
			// z distance
			particles.at(ip)->zdist = fabs(particles.at(ip)->x(2) - particles.at(one_surf_idx)->x(2));
			particles.at(ip)->zdist -= std::floor(particles.at(ip)->zdist/Lz + 0.5)*Lz;
			particles.at(ip)->zdist = fabs(particles.at(ip)->zdist);
		}
	}
}
