#ifndef FRAME_H
#define FRAME_H

#include <vector>

#include "Particle.h"

class Frame
{
	public:
		unsigned int tstep;
		std::vector<std::shared_ptr<Particle> > particle_ptrs;
		double xmin,
			   xmax,
			   ymin,
			   ymax,
			   zmin,
			   zmax;

		Frame(unsigned int tstep, std::vector<std::shared_ptr<Particle> >& particle_ptrs, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax): tstep(tstep), particle_ptrs(particle_ptrs), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax) {};
};


#endif
