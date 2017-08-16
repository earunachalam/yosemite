#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int main()
{
	size_t N = 10000, numts = 10000, sampf = 100;
	std::vector<double> x(N), v(N), a(N), xrec;
	xrec.reserve((size_t) numts/sampf);

	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
	std::uniform_real_distribution<double> xd(-3,3), vd(-0.05,0.05);

	for(auto& elem: x) elem = xd(rng);
	for(auto& elem: v) elem = vd(rng);

	double dt = 0.002;
	double tol = -5.00e-15;

	unsigned int ctr = 0;
	for (double wcenter = 0; wcenter < 10.0; wcenter += 1.0)
	{
		for (size_t ts = 1; ts < numts; ++ts)
		{
			for (size_t i = 0; i < N; ++i)
			{
				if (x.at(i) < 0) a.at(i) = 0; 
				else if (x.at(i) > 10) a.at(i) = 0;
				else a.at(i) = -1.00 + 0.5*std::pow(x.at(i) - wcenter,2); //fabs(x.at(i)) < 1.00 ? 0 : x.at(i) - pow(x.at(i),3.00);
				x.at(i) += v.at(i)*dt;
				v.at(i) += a.at(i)*dt;
			}

			if (ts>1000 && (ts%sampf)==0) for (size_t i = 0; i < N; ++i) xrec.push_back(x.at(i));
		}

		std::string fname = "bias"
		FILE* fp = fopen("nobias.dat", "w");
		for (size_t i = 0; i < xrec.size(); ++i) fprintf(fp, "%f ", xrec.at(i));
		fclose(fp);
	}
	
	return 0;
}
