#include <algorithm>
#include <iostream>
#include <cassert>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>

#include "voro++.hh"

#include "Histogram.h"


#define NCHAR 			2000

struct molec
{
	unsigned int 				id = -1;
	double 						x, y, z;
	double 						q6 = 0.00;
	std::vector<unsigned int>	neighbors,
								neighbors_xtalline;
};

struct cluster
{
	std::vector<unsigned int> 	members,
								branchpoints;
};

struct frame
{
	unsigned int 	id;
	double			xmin,
					xmax,
					ymin,
					ymax,
					zmin,
					zmax,
					dx,
					dy,
					dz;

	std::vector<std::shared_ptr<molec> >	molecs;
	std::vector<std::shared_ptr<cluster> >	clusters;
};

int main(int argc, char* argv[])
{
	FILE* ifp; // input file
	std::string statsdir;
	double threshold_q6 = 0.00;

	if (argc < 2)
	{
		std::cerr << "Error: insufficient arguments provided." << std::endl;
		abort();
	}
	else
	{
		ifp = fopen(argv[1], "r");
		if (ifp == nullptr)
		{
			printf("File open operation failed: check cwd.\n");
			abort();
		}

		statsdir = std::string(argv[2]);
		threshold_q6 = std::atof(argv[3]);
	}

	std::vector<std::shared_ptr<frame> > frame_ptrs;
	std::shared_ptr<frame> f;
	
	while (!feof(ifp))
	{
		char currline[NCHAR];
		if(fgets(currline, NCHAR, ifp))
		{
			char*						token = strtok(currline, " ");
			std::vector<std::string> 	tokens;
			while(token != NULL)
			{
				std::string newtoken(token);
				tokens.push_back(std::string(token));
				token = strtok(NULL, " \n");
			}

			if (tokens.at(0) == "FRAME")
			{
				// create latest frame
				frame_ptrs.push_back(std::make_shared<frame>());
				f = frame_ptrs.back();
				
				f->id 	= std::stoi(tokens.at(1));
				std::cout << "Reading frame " << f->id << "\r" << std::flush;
				f->xmin	= std::stod(tokens.at(2));
				f->xmax	= std::stod(tokens.at(3));
				f->ymin	= std::stod(tokens.at(4));
				f->ymax	= std::stod(tokens.at(5));
				f->zmin	= std::stod(tokens.at(6));
				f->zmax	= std::stod(tokens.at(7));
				f->dx 	= f->xmax - f->xmin;
				f->dy 	= f->ymax - f->ymin;
				f->dz 	= f->zmax - f->zmin;
			}
			else
			{
				f->molecs.push_back(std::make_shared<molec>());
				auto m = f->molecs.back();
				m->id 	= std::stoi(tokens.at(0));
				m->x 	= std::stod(tokens.at(1));
				m->y 	= std::stod(tokens.at(2));
				m->z 	= std::stod(tokens.at(3));
				m->q6 	= std::stod(tokens.at(4));

				for (unsigned int idx = 5; idx < tokens.size(); idx += 3)
				{
					unsigned int id = std::stoi(tokens.at(idx));
					double dist 	= std::stoi(tokens.at(idx+1));
					double q6 		= std::stoi(tokens.at(idx+2));
					
					if (dist < 3.7)
					{
						m->neighbors.push_back(id);

						if (q6 >= threshold_q6)
						{
							m->neighbors_xtalline.push_back(id);
						}
					}
				}
			}
		}
	}

	size_t n_molecs = f->molecs.size();

	//#pragma omp parallel for
	for (size_t iframe = 0; iframe < frame_ptrs.size(); ++iframe)
	{
		auto f = frame_ptrs.at(iframe);			// pointer to current frame
		auto m = f->molecs;						// pointers to all molecs in current frame

		// carry out DFS of molecs in each frame to identify icelike connected components (clusters)
		
		// get index of element (of the current frame's molec vector m) which has given id
		// simple version, very fast BUT must make sure that particle IDs start at 1...
		auto id2idx = [&m](unsigned int id){ return id-1; };
		// general version, slow BUT works for any starting value for particle IDs
		//auto id2idx = [&m](unsigned int id){ for (size_t i = 0; i < m.size(); ++i){ if (m.at(i)->id == id) return i; } return SIZE_MAX;};
		
		// vector holding ids of molecules that must be considered
		std::vector<unsigned int> molecs_remaining(n_molecs);
		std::iota(molecs_remaining.begin(), molecs_remaining.end(), 1);

		// vector holding ids of molecules that are already in a cluster
		std::vector<unsigned int> already_in_cluster;

		unsigned int curr_cluster_idx = -1;		// to hold index of cluster to which molec must be added

		// while there are still molecules that have not yet been assigned to a cluster (if q6 abov threshold) or neglected due to low q6
		while (!molecs_remaining.empty())
		{
			f->clusters.push_back(std::make_shared<cluster>());	// add new cluster
			curr_cluster_idx = f->clusters.size() - 1;			// index of new cluster in list of clusters - this is the default cluster to which molecs being considered are added, i.e. as long they are not (by virtue of common neighbors) members of a previously created cluster

			// neglect all remaining molecs which are not icelike
			while (m.at(id2idx(molecs_remaining.front()))->q6 < threshold_q6)
			{
				molecs_remaining.erase(molecs_remaining.begin());
				if (molecs_remaining.empty()) break;
			}

			// if no more icelike molecs remaining, write cluster data and move onto next frame
			if (molecs_remaining.empty()) break;

			// add the new molec as a 'seed' for the new cluster. Also add to branchpoint list (needed for DFS)
			f->clusters.at(curr_cluster_idx)->members.push_back(molecs_remaining.front());
			f->clusters.at(curr_cluster_idx)->branchpoints.push_back(molecs_remaining.front());
			already_in_cluster.push_back(molecs_remaining.front());
			
			// do not consider the current molec (the new seed) in a future loop
			molecs_remaining.erase(molecs_remaining.begin());

			// true if there are no more branchpoints to consider in the current connected component
			bool curr_branchpoint_list_empty = f->clusters.at(curr_cluster_idx)->branchpoints.empty();

			// while there are still more branchpoints to consider in the current connected component, i.e. the chain of neighbors is not yet terminated
			while (!curr_branchpoint_list_empty)
			{
				// last branchpoint in vecotr, i.e. id of molec whose neighbors we must now consider
				unsigned int curr_branchpoint 				= f->clusters.at(curr_cluster_idx)->branchpoints.back();
				// index of molec within vector m of molecs that has id == last branchpoint
				unsigned int curr_branchpoint_idx_in_m 		= id2idx(curr_branchpoint);

				// loop over neighbors of current branchpoint molec
				for (auto& n: m.at(curr_branchpoint_idx_in_m)->neighbors)
				{
					// get index of molec within vector m of molecs that has id == n
					unsigned int nidx = id2idx(n);

					// if current neighbor is not ice-like, skip
					if (m.at(nidx)->q6 < threshold_q6) continue;

					unsigned int member_cluster_idx // to hold index of cluster of which the current neighbor is a part (if any)
						= curr_cluster_idx;			// default is the current cluster
					bool add_to_current = true;		// add to current cluster (if current neighbor not part of previous cluster)
					for (unsigned int ic = 0; ic < f->clusters.size(); ++ic)
					{
						for (auto& member: f->clusters.at(ic)->members)
						{
							if (n == member)
							{
								member_cluster_idx = ic;
								add_to_current = false;
							}
						}
					}

					// if current neighbor is part of a previously created cluster
					if (member_cluster_idx != curr_cluster_idx)
					{
						// merge the current cluster and the previously created cluster
						
						f->clusters.at(member_cluster_idx)->branchpoints.insert(
								f->clusters.at(member_cluster_idx)->branchpoints.end(),
								f->clusters.at(curr_cluster_idx)->branchpoints.begin(),
								f->clusters.at(curr_cluster_idx)->branchpoints.end());
								
						f->clusters.at(member_cluster_idx)->members.insert(
								f->clusters.at(member_cluster_idx)->members.end(),
								f->clusters.at(curr_cluster_idx)->members.begin(),
								f->clusters.at(curr_cluster_idx)->members.end());
						
						f->clusters.erase(f->clusters.begin() + curr_cluster_idx);
						if (curr_cluster_idx < member_cluster_idx) 	curr_cluster_idx = member_cluster_idx - 1;
						else 										curr_cluster_idx = member_cluster_idx;
					}
					
					// if current neighbor does not bridge the current cluster with a previously created cluster
					if (add_to_current) 
					{
						// add current neighbor to most recently created cluster

						f->clusters.at(curr_cluster_idx)->members.push_back(n);
						f->clusters.at(curr_cluster_idx)->branchpoints.push_back(n);
						
						// get index of molec within vector molecs_remaining that has id == n
						for (size_t i = 0; i < molecs_remaining.size(); ++i)
						{
							if (m.at(id2idx(molecs_remaining.at(i)))->id == n)
							{
								molecs_remaining.erase(molecs_remaining.begin()+i);
								break;
							}
						}
					}
				}
				
				// erase current branchpoint since we have found all of its neighbors and made them branchpoints
				unsigned int curr_branchpoint_idx_in_itself = -1;
				for (unsigned int ibp = 0; ibp < f->clusters.at(curr_cluster_idx)->branchpoints.size(); ++ibp)
				{
					if (f->clusters.at(curr_cluster_idx)->branchpoints.at(ibp) == curr_branchpoint)
					{
						curr_branchpoint_idx_in_itself = ibp;
						break;
					}
				}
				f->clusters.at(curr_cluster_idx)->branchpoints.erase(
						(f->clusters.at(curr_cluster_idx)->branchpoints.begin())+curr_branchpoint_idx_in_itself);
				
				// again, check if current branchpoint list is empty
				curr_branchpoint_list_empty = f->clusters.at(curr_cluster_idx)->branchpoints.empty();
			}
		}

		std::string curr_frame_fname = statsdir + "/f" + std::to_string(f->id) + ".dat";
		std::cout << curr_frame_fname << std::endl;
		std::ofstream curr_frame_clusters(curr_frame_fname);
		for (auto& c: f->clusters)
		{
			for (auto& member: c->members) curr_frame_clusters << member << " "; // subtract 1 to get vmd indices
			curr_frame_clusters << "\n";
		}
		curr_frame_clusters.close();
		std::cout << iframe << " completed." << std::endl;
	}

	return 0;
}
