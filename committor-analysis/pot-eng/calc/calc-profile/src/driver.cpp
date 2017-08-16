#include <cstdio>
#include <cstdlib>

#include "Histogram.h"

int main()
{
	FILE* ifp = fopen("zdist_q6.dat", "r");
	if (ifp == nullptr) { printf("File open operation failed: check cwd.\n"); abort(); }

	Histogram h(0.0, 29.0, 0.1);

	while (!feof(ifp))
	{
		double zdist, opval;
		
		// assumes format is z-distance from surface, order param value
		fscanf(ifp, "%lf %lf", &zdist, &opval);
	
		h.addData(zdist, opval);
	}

	std::vector<double>
		bins 	= h.g_bins(),
		counts 	= h.g_counts(),
		values 	= h.g_values();

	FILE* ofp = fopen("profile_q6.dat","w");

	for (unsigned int ibin = 0; ibin < bins.size(); ++ibin)
	{
		fprintf(ofp, "%lf %lf %lf\n", bins.at(ibin), counts.at(ibin), values.at(ibin));
	}

	return 0;
}
