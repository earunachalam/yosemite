#include "Histogram.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "except_handler.h"

Histogram::Histogram(const double iarg_min, const double iarg_max, const double iarg_binwidth)
{
	try
	{
		if (iarg_min >= iarg_max)
		{
			std::string message = std::string(__FILE__) + ": " + std::to_string(__LINE__) + ": Minimum value must be less than maximum value.";
			throw std::logic_error(message);
		}

		this->min = iarg_min;
		this->max = iarg_max;
		this->binwidth = iarg_binwidth;

		this->range = this->max - this->min;
		this->nBins = static_cast<int>(range/binwidth) + 1;

		this->bins.resize(nBins);
		this->counts.resize(nBins);
		this->values.resize(nBins);

		std::iota((this->bins).begin(), (this->bins).end(), 0.0);
		std::transform((this->bins).begin(), (this->bins).end(), (this->bins).begin(),[&](double i){return binwidth*(i+0.5)+this->min;});

		std::fill((this->counts).begin(), (this->counts).end(), 0.0);
		std::fill((this->values).begin(), (this->values).end(), 0.0);
	}
	catch(const std::exception &ex)
	{
		#pragma omp critical
		{
			__userdef_fatal_error(ex);
		}
	}
}

Histogram::Histogram(const Histogram &obj)
{
	this->min = obj.min;
	this->max = obj.max;
	this->binwidth = obj.binwidth;

	this->range = obj.range;
	this->nBins = obj.nBins;

	this->bins.resize(nBins);
	this->counts.resize(nBins);
	this->values.resize(nBins);

	this->bins = obj.bins;
	std::fill((this->counts).begin(), (this->counts).end(), 0.0);
	std::fill((this->values).begin(), (this->values).end(), 0.0);
}

void Histogram::addData(const double locationWithinRange, const double data)
{
	int bin = static_cast<int>((locationWithinRange - min)/this->binwidth);
	if ((bin >= 0) && (bin < this->nBins))
	{
		this->counts.at(bin) += 1.0;
		this->values.at(bin) += data;
	}
}

std::vector<double> Histogram::g_bins()
{
	return this->bins;
}

std::vector<double> Histogram::g_counts()
{
	return this->counts;
}

std::vector<double> Histogram::g_values()
{
	return this->values;
}
