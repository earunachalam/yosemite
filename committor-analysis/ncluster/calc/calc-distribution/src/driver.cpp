#include <cstdio>
#include <cstdlib>

#include "Histogram.h"

int main()
{
	FILE* ifp = fopen("zdist_q6.dat", "r");
	if (ifp == nullptr) { printf("File open operation failed: check cwd.\n"); abort(); }

	Histogram h(-1.0, 1.0, 0.01);

	while (!feof(ifp))
	{
		double zdist, opval;

		// assumes format is z-distance from surface, order param value
		if (!fscanf(ifp, "%lf %lf", &zdist, &opval))
		{ printf("Error reading file.\n"); abort(); }

		if (zdist < 4.500) h.addData(opval);
	}

	std::vector<double>
		bins 	= h.g_bins(),
		counts 	= h.g_counts(),
		values 	= h.g_values();

	FILE* ofp = fopen("contact_q6.dat","w");

	for (unsigned int ibin = 0; ibin < bins.size(); ++ibin)
	{
		fprintf(ofp, "%lf %lf %lf\n", bins.at(ibin), counts.at(ibin), values.at(ibin));
	}

	return 0;
}
