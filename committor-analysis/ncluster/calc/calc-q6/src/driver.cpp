#include <csignal>
#include <iostream>
#include <omp.h>

#include "Frame.h"
#include "InputParser.h"
#include "Particle.h"

#include "compute.h"

int main(int argc, char* argv[])
{
	FILE* fp_traj = fopen("traj.lammpstrj", "r");
	if (fp_traj == nullptr) printf("Error opening file: check cwd.\n");

	// file for q6 profile
	// particle fmt: zdist-from-surf q6val
	// frames not separated
	FILE* fp_zdist_q6 = fopen("zdist_q6.dat", "w");
	if (fp_zdist_q6 == nullptr) printf("Error opening file: check cwd.\n");

	// file for neighbors (used in clustering)
	// particle fmt: x-coord y-coord z-coord q6val neigh1id neigh1q6val neigh2id neigh2q6val ...
	// frames separated by: FRAME %framenum
	FILE* fp_x_y_z_q6 = fopen("x_y_z_q6.dat", "w");
	if (fp_x_y_z_q6 == nullptr) printf("Error opening file: check cwd.\n");

	unsigned int startframe = 0,
				 endframe = 1000,
				 skipframe = 1,
				 nframe,
				 nwat = 2000,
				 nother = 416,
				 ntotal;
	bool scaledcoords = false; // default

	InputParser ip(argc, argv);
	
	startframe = ip.cmdOptionExists("startframe") ? stoi(ip.getCmdOption("startframe")) : 0;
	endframe = ip.cmdOptionExists("endframe") ? stoi(ip.getCmdOption("endframe")) : 0;
	skipframe = ip.cmdOptionExists("skipframe") ? stoi(ip.getCmdOption("skipframe")) : 0;
	nframe = (endframe - startframe) / skipframe;
	nwat = ip.cmdOptionExists("nwat") ? stoi(ip.getCmdOption("nwat")) : 0;
	nother = ip.cmdOptionExists("nother") ? stoi(ip.getCmdOption("nother")) : 0;
	ntotal = nwat + nother;
	scaledcoords = ip.cmdOptionExists("scaled");

	std::vector<std::shared_ptr<Frame> > frame_ptrs(nframe);

	for (unsigned int ifr = 0; ifr < endframe; ++ifr)
	{
		unsigned int frameTimestep, nParticles;
		double xmin, xmax, ymin, ymax, zmin, zmax;

		// if at selected frame
		bool useframe = (ifr >= startframe) && ((ifr-startframe)%skipframe == 0);

		char ignore[1000];

		if (fgets(ignore, 1000, fp_traj) == NULL)
		{printf("Error reading header at frame %d\n", ifr); abort();}

		if (fscanf(fp_traj, "%d ", &frameTimestep) != 1)
		{printf("Error reading header at frame %d\n", ifr); abort();}

		if (fgets(ignore, 1000, fp_traj) == NULL)
		{printf("Error reading header at frame %d\n", ifr); abort();}

		if (fscanf(fp_traj, "%d ", &nParticles) != 1)
		{printf("Error reading number of pptrs at frame %d\n", ifr); abort();}
		else if (nParticles != (nwat + nother))
		{printf("Error: Npptrs in lammpstrj != Npptrs at frame %d\n", ifr); abort();}

		if (fgets(ignore, 1000, fp_traj) == NULL)
		{printf("Error reading header at frame %d\n", ifr); abort();}

		if (fscanf(fp_traj, "%lf %lf ", &xmin, &xmax) != 2)
		{printf("Error reading x cell vectors at frame %d\n", ifr); abort();}

		if (fscanf(fp_traj, "%lf %lf ", &ymin, &ymax) != 2)
		{printf("Error reading y cell vectors at frame %d\n", ifr); abort();}

		if (fscanf(fp_traj, "%lf %lf ", &zmin, &zmax) != 2)
		{printf("Error reading z cell vectors at frame %d\n", ifr); abort();}

		if (fgets(ignore, 1000, fp_traj) == NULL)
		{printf("Error reading header at frame %d\n", ifr); abort();}

		// box dims
		double Lx = xmax - xmin, Ly = ymax - ymin, Lz = zmax - zmin;

		// particle properties
		unsigned int id, type;
		double x, y, z;

		if (useframe)
		{
			std::vector<std::shared_ptr<Particle> > pptrs(ntotal);

			for (unsigned int pidx = 0; pidx < ntotal; ++pidx)
			{
				if (fscanf(fp_traj, "%d %d %lf %lf %lf ", &id, &type, &x, &y, &z) != 5)
				{printf("Error reading particle %d at frame %d\n", pidx, ifr); abort();}

				if (scaledcoords) { x *= Lx; y *= Ly; z *= Lz; }

				pptrs.at(id-1) = std::make_shared<Particle>(id, type, x, y, z);
			}

			frame_ptrs.at(ifr-startframe) = std::make_shared<Frame>(ifr, pptrs, xmin, xmax, ymin, ymax, zmin, zmax);
		}
		else
		{
			for (unsigned int pidx = 0; pidx < ntotal; ++pidx)
			{
				if (fscanf(fp_traj, "%d %d %lf %lf %lf ", &id, &type, &x, &y, &z) != 5)
				{printf("Error reading particle %d at frame %d\n", pidx, ifr); abort();}
			}
		}
	}
	
	#pragma omp parallel for
	for (unsigned int ifr = 0; ifr < nframe; ++ifr)
	{
		std::shared_ptr<Frame> f = frame_ptrs.at(ifr);
		compute_distance_surface(f->particle_ptrs, 2001, f->zmin, f->zmax);
		compute_neighbors_orthorhombic(f->particle_ptrs, f->xmin, f->xmax, f->ymin, f->ymax, f->zmin, f->zmax, 6);
		compute_q6(f->particle_ptrs);

		printf("Frame %d complete\n", ifr);
	}

	for (unsigned int ifr = 0; ifr < nframe; ++ifr)
	{
		std::shared_ptr<Frame> f = frame_ptrs.at(ifr);
		fprintf(fp_x_y_z_q6, "FRAME %u %lf %lf %lf %lf %lf %lf\n", ifr, f->xmin, f->xmax, f->ymin, f->ymax, f->zmin, f->zmax);

		for (unsigned int jp = 0; jp < ntotal; ++jp)
		{
			std::shared_ptr<Particle> p = f->particle_ptrs.at(jp);
			if (p->type == WATERTYPE)
			{
				// file for profile
				fprintf(fp_zdist_q6, "%lf %lf\n", p->zdist, p->q6);
				
				// file for neighbors (used in clustering)
				fprintf(fp_x_y_z_q6, "%u %lf %lf %lf %lf ", p->id, p->x(0), p->x(1), p->x(2), p->q6);
				for (unsigned int kn = 0; kn < p->neighbor_ptrs.size(); ++kn)
				{
					std::shared_ptr<Neighbor> kn_ptr = p->neighbor_ptrs.at(kn);
					fprintf(fp_x_y_z_q6, "%u %lf %lf ", kn_ptr->selfptr->id, kn_ptr->dist, kn_ptr->selfptr->q6);
				}
				fprintf(fp_x_y_z_q6, "\n");
			}
		}
	}
	return 0;
}
