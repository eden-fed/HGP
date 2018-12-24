#pragma once

#include <complex>
#include <limits>



namespace LIMIT
{
	bool hasAllFloatingPointSpecialSupport();
	bool isNAN(double num);
	bool isInfinity(double num);
	bool isFinite(double num);
	bool isNAN(float num);
	bool isInfinity(float num);
	bool isFinite(float num);
	bool isNAN(const std::complex<double>& num);
	bool isInfinity(const std::complex<double>& num);
	bool isFinite(const std::complex<double>& num);
	bool isNAN(const std::complex<float>& num);
	bool isInfinity(const std::complex<float>& num);
	bool isFinite(const std::complex<float>& num);

	void setInfinity(double& num);
	void setInfinity(float& num);
	void setInfinity(std::complex<double>& num);
	void setInfinity(std::complex<float>& num);

}

