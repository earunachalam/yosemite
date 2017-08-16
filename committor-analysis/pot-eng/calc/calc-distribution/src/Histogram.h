#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "excepthandle.h"

class Histogram
{
protected:

	double min;
	double max;
	double range;
	double binwidth;

	int nBins;
	std::vector<double> bins;
	std::vector<double> counts;
	std::vector<double> values;

public:

	Histogram(const double iarg_min, const double iarg_max, const double iarg_binwidth);

	Histogram(const Histogram &obj);

	void addData(const double locationWithinRange, const double data = 1.0);

	std::vector<double> g_bins();

	std::vector<double> g_counts();

	std::vector<double> g_values();

};

#endif
