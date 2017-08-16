#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <memory>

class Particle;

class Neighbor
{
	std::shared_ptr<Particle> m_ptr;
	double dist = -1.00;
};

#endif
