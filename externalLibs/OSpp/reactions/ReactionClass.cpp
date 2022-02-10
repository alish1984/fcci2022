/* 
 * File:   TransportPolicy_CHEMKIN.cpp
 * Author: acuoci
 * 
 * Created on October 2, 2012, 7:02 PM
 */

#include <iomanip>
#include "ReactionClass.h"

std::string GetUnitsOfKineticConstantsSI(const double lambda)
{
	std::string units;

		 if (lambda == 1.)	units += "[1/s]";
	else if (lambda == 2.)	units += "[m3/kmol/s]";
	else if (lambda == 3.)	units += "[m6/kmol2/s]";
	else if (lambda == 4.)	units += "[m9/kmol3/s]";
	else
	{
		std::stringstream exponent;
		std::stringstream exponent3;
		exponent << lambda - 1.;
		exponent3 << 3.*(lambda - 1.);
		units += "[m" + exponent3.str() + "/kmol" + exponent.str() + "/s]";
	}

	return units;
}

std::string GetUnitsOfKineticConstantsCGS(const double lambda)
{
	std::string units;

		 if (lambda == 1.)	units += "[1/s]";
	else if (lambda == 2.)	units += "[cm3/mol/s]";
	else if (lambda == 3.)	units += "[cm6/mol2/s]";
	else if (lambda == 4.)	units += "[cm9/mol3/s]";
	else
	{
		std::stringstream exponent;
		std::stringstream exponent3;
		exponent << lambda - 1.;
		exponent3 << 3.*(lambda - 1.);
		units += "[cm" + exponent3.str() + "/mol" + exponent.str() + "/s]";
	}

	return units;
}

void GetKineticConstantString(const double A, const double beta, const double E, const double lambda, std::string& line_kForwardStringSI, std::string& line_kForwardStringCGS)
{
	line_kForwardStringSI  = "k = ";
	line_kForwardStringCGS = "k = ";

	// Frequency factor
	{
		std::stringstream ANumber_;
		ANumber_ << std::scientific << std::setprecision(6) << A;
		line_kForwardStringSI += ANumber_.str();
	}
	{
		std::stringstream ANumber_;
		ANumber_ << std::scientific << std::setprecision(6) << A/std::pow(1.e3, 1.-lambda);
		line_kForwardStringCGS += ANumber_.str();
	}
		
	// Temperature exponent
	if (beta != 0.) 
	{
		std::stringstream betaNumber_;
		betaNumber_ << std::fixed << std::setprecision(3) << beta;
		line_kForwardStringSI  += " T^" + betaNumber_.str();
		line_kForwardStringCGS += " T^" + betaNumber_.str();
	}

	// Activation energy
	if (E != 0.)    
	{
		std::stringstream ENumber_;
		ENumber_ << std::scientific << std::setprecision(6) << -E;
		line_kForwardStringSI += " exp(" + ENumber_.str() + "/RT)";
	}
	{
		std::stringstream ENumber_;
		ENumber_ << std::fixed << std::setprecision(2) << -E/(1.e3*Conversions::J_from_cal);
		line_kForwardStringCGS += " exp(" + ENumber_.str() + "/RT)";
	}
}
