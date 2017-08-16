// default include files from sample program
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

// LAMMPS headers
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"

// additional include files
#include <iostream>
#include <memory>
#include <vector>

// additional headers
#include "compute.h"
#include "Particle.h"

// ignore warnings regarding strings in sample program
#pragma GCC diagnostic ignored "-Wwrite-strings"



#define NNEIGH_TRACK 30

using namespace LAMMPS_NS;

int main(int narg, char **arg)
{
	// setup MPI and various communicators
	// driver runs on all procs in MPI_COMM_WORLD
	// comm_lammps only has 1st P procs (could be all or any subset)

	MPI_Init(&narg,&arg);

	if (narg != 3)
	{
		printf("Syntax: simpleCC P in.lammps\n");
		exit(1);
	}

	int me,nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int nprocs_lammps = atoi(arg[1]);
	if (nprocs_lammps > nprocs)
	{
		if (me == 0)
		{
			printf("ERROR: LAMMPS cannot use more procs than available\n");
		}

		MPI_Abort(MPI_COMM_WORLD,1);
	}

	int lammps;
	if (me < nprocs_lammps)		lammps = 1;
	else 						lammps = MPI_UNDEFINED;

	MPI_Comm comm_lammps;
	MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);

	// open LAMMPS input script

	FILE *fp = NULL;
	if (me == 0)
	{
		fp = fopen(arg[2],"r");
		if (fp == NULL)
		{
			printf("ERROR: Could not open LAMMPS input script\n");
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	// run the input script through LAMMPS one line at a time until end-of-file
	// driver proc 0 reads a line, Bcasts it to all procs
	// (could just send it to proc 0 of comm_lammps and let it Bcast)
	// all LAMMPS procs call input->one() on the line

	LAMMPS *lmp = NULL;
	if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);

	int n;
	char line[1024];
	while (1)
	{
		if (me == 0)
		{
			if (fgets(line,1024,fp) == NULL) 	n = 0;
			else 								n = strlen(line) + 1;

			if (n == 0) fclose(fp);
		}

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (n == 0) break;
		MPI_Bcast(line, n, MPI_CHAR, 0, MPI_COMM_WORLD);
		if (lammps == 1) lammps_command(lmp, line);
	}


	// construct particle objects from lammps data

	unsigned int natoms = static_cast<unsigned int>(lmp->atom->natoms);
	int* id_raw = new int[natoms];
	int* type_raw = new int[natoms];
	double* x_raw = new double[3*natoms];
	
	lammps_gather_atoms(lmp, "id", 0, 1, id_raw);
	std::vector<int> id(id_raw, id_raw+natoms);
	
	lammps_gather_atoms(lmp, "type", 0, 1, type_raw);
	std::vector<int> type(type_raw, type_raw+natoms);
	
	lammps_gather_atoms(lmp, "x", 1, 3, x_raw);
	std::vector<double> x(x_raw, x_raw+3*natoms);

	std::vector<std::shared_ptr<Particle> > particles;
	for (unsigned int i = 0; i < natoms; ++i)
	{
		std::vector<double> curr_x = {x.at(i*3), x.at(i*3 + 1), x.at(i*3 + 2)};
		particles.push_back(std::make_shared<Particle>(static_cast<unsigned int>(id.at(i)), static_cast<unsigned int>(type.at(i)), curr_x));
	}

	double xmin = *((double*) lammps_extract_global(lmp, "boxxlo"));
	double xmax = *((double*) lammps_extract_global(lmp, "boxxhi"));
	double ymin = *((double*) lammps_extract_global(lmp, "boxylo"));
	double ymax = *((double*) lammps_extract_global(lmp, "boxyhi"));
	double zmin = *((double*) lammps_extract_global(lmp, "boxzlo"));
	double zmax = *((double*) lammps_extract_global(lmp, "boxzhi"));

	compute_neighbors_orthorhombic(particles, xmin, xmax, ymin, ymax, zmin, zmax, NNEIGH_TRACK);

	// calculate q6 profile for system

	std::vector<computed_particlewise_vecdouble> q6_all_allparticles(particles.size());
	// first 13 elements of vector to contain q6m values and the final element to contain the q6 value
	
	compute_q6(particles, q6_all_allparticles);
	

	//for (auto& elem: particles)
	//{
		//if (elem->m_id != 43) continue;

		//printf("id = %u\n", elem->m_id);
		//printf("pos = %lf, %lf, %lf\n", elem->m_x[0], elem->m_x[1], elem->m_x[2]);
		//for (auto& elem2: elem->m_neighbor_ptrs)
		//{
			//printf("id = %u, type = %u, dist = %lf\n", elem2->m_ptr->m_id, elem2->m_ptr->m_type, elem2->m_dist);
		//}
		//return 0;
	//}

	//if (lammps == 1) {
	//lmp->input->one("run 10");

	//double epsilon = 0.1;
	//x[0] += epsilon;
	//lammps_scatter_atoms(lmp,"x",1,3,x);

	//// these 2 lines are the same

	//// lammps_command(lmp,"run 1");
	//lmp->input->one("run 1");
	//}

	//// extract force on single atom two different ways

	//if (lammps == 1) {
	//double **f = (double **) lammps_extract_atom(lmp,"f");
	//printf("Force on 1 atom via extract_atom: %g\n",f[0][0]);

	//double *fx = (double *) lammps_extract_variable(lmp,"fx","all");
	//printf("Force on 1 atom via extract_variable: %g\n",fx[0]);
	//}

	//// use commands_string() and commands_list() to invoke more commands

	//char *strtwo = "run 10\nrun 20";
	//if (lammps == 1) lammps_commands_string(lmp,strtwo);

	//char *cmds[2];
	//cmds[0] = "run 10";
	//cmds[1] = "run 20";
	//if (lammps == 1) lammps_commands_list(lmp,2,cmds);

	//// delete all atoms
	//// create_atoms() to create new ones with old coords, vels
	//// initial thermo should be same as step 20

	//int *type = NULL;

	//if (lammps == 1) {
	//int natoms = static_cast<int> (lmp->atom->natoms);
	//type = new int[natoms];
	//for (int i = 0; i < natoms; i++) type[i] = 1;

	//lmp->input->one("delete_atoms group all");
	//lammps_create_atoms(lmp,natoms,NULL,type,x,v,NULL,0);
	//lmp->input->one("run 10");
	//}

	//delete [] x;
	//delete [] v;
	//delete [] type;

	//// close down LAMMPS

	//delete lmp;

	//// close down MPI

	//if (lammps == 1) MPI_Comm_free(&comm_lammps);
	//MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Finalize();
}
