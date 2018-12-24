#include "Limit_Util.h"




bool LIMIT::hasAllFloatingPointSpecialSupport()
{
	if(!std::numeric_limits<double>::has_infinity)
	{
		return false;
	}
	if(!std::numeric_limits<double>::has_quiet_NaN)
	{
		return false;
	}
	if(!std::numeric_limits<double>::has_signaling_NaN)
	{
		return false;
	}

	if(!std::numeric_limits<float>::has_infinity)
	{
		return false;
	}
	if(!std::numeric_limits<float>::has_quiet_NaN)
	{
		return false;
	}
	if(!std::numeric_limits<float>::has_signaling_NaN)
	{
		return false;
	}
	return true;
}

bool LIMIT::isNAN(double num)
{
	if(num == std::numeric_limits<double>::quiet_NaN() || num == std::numeric_limits<double>::signaling_NaN())
	{
		return true;
	}
	return false;
}

bool LIMIT::isInfinity(double num)
{
	if(num == std::numeric_limits<double>::infinity() || -num == std::numeric_limits<double>::infinity())
	{
		return true;
	}
	return false;
}

bool LIMIT::isFinite(double num)
{
	if(isNAN(num) || isInfinity(num))
	{
		return false;
	}
	return true;
}

bool LIMIT::isNAN(float num)
{
	if(num == std::numeric_limits<float>::quiet_NaN() || num == std::numeric_limits<float>::signaling_NaN())
	{
		return true;
	}
	return false;
}

bool LIMIT::isInfinity(float num)
{
	if(num == std::numeric_limits<float>::infinity() || -num == std::numeric_limits<float>::infinity())
	{
		return true;
	}
	return false;
}

bool LIMIT::isFinite(float num)
{
	if(isNAN(num) || isInfinity(num))
	{
		return false;
	}
	return true;
}


bool LIMIT::isNAN(const std::complex<double>& num)
{
	if(isNAN(real(num)) || isNAN(imag(num)))
	{
		return true;
	}
	return false;
}

bool LIMIT::isInfinity(const std::complex<double>& num)
{
	if(isInfinity(real(num)) || isInfinity(imag(num)))
	{
		return true;
	}
	return false;
}

bool LIMIT::isFinite(const std::complex<double>& num)
{
	if(isFinite(real(num)) && isFinite(imag(num)))
	{
		return true;
	}
	return false;
}


bool LIMIT::isNAN(const std::complex<float>& num)
{
	if(isNAN(real(num)) || isNAN(imag(num)))
	{
		return true;
	}
	return false;
}

bool LIMIT::isInfinity(const std::complex<float>& num)
{
	if(isInfinity(real(num)) || isInfinity(imag(num)))
	{
		return true;
	}
	return false;
}

bool LIMIT::isFinite(const std::complex<float>& num)
{
	if(isFinite(real(num)) && isFinite(imag(num)))
	{
		return true;
	}
	return false;
}

void LIMIT::setInfinity(double& num)
{
	num = std::numeric_limits<double>::infinity();
}


void LIMIT::setInfinity(float& num)
{
	num = std::numeric_limits<float>::infinity();
}


void LIMIT::setInfinity(std::complex<double>& num)
{
	num = std::complex<double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
}


void LIMIT::setInfinity(std::complex<float>& num)
{
	num = std::complex<float>(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity());
}



