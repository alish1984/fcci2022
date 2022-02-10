/* 
 * File:   TransportPolicy_CHEMKIN.h
 * Author: acuoci
 *
 * Created on October 2, 2012, 7:02 PM
 */

#ifndef REACTIONCLASS_H
#define	REACTIONCLASS_H

#include <vector>

class ThirdBodyPolicy
{
	double efficiency(const std::vector<double>& c) { return 2.; }
};

class NoThirdBodyPolicy
{
	double efficiency(const std::vector<double>& c) { return 1.; }
};

class IrreversiblePolicy
{
	// Stoichiometric coefficients
	const std::vector<double>& reactant_nu() const { return reactant_nu_;}
	const std::vector<double>& product_nu() const { return product_nu_;}
	const std::vector<unsigned int>& reactant_nu_indices() const { return reactant_nu_indices_;}
	const std::vector<unsigned int>& product_nu_indices() const { return product_nu_indices_;}

	// Set stoichiometric coefficients
	void set_reactant_nu(const unsigned int i, const double value) { reactant_nu_[i] = value; }
	void set_product_nu(const unsigned int i, const double value)  { product_nu_[i] = value; }

	double Product(const std::vector<double>& c) { return 2.; }

	private:
};

class ReversiblePolicy
{
	double Product(const std::vector<double>& c) { return 2.; }
};

class ReactionRateSimplePolicy
{
	double Product(const std::vector<double>& c) { return 2.; }
};

class ReactionRateFallOffPolicy
{
	double Product(const std::vector<double>& c) { return 2.; }
};

class ReactionRateChmicallyActivatedBimolecularPolicy
{
	double Product(const std::vector<double>& c) { return 2.; }
};


template
<
	template <class> class Reversibility,
	template <class> class ThirdBody,
	template <class> class ReactionRate
>
class ReactionClass
{
public:

	ReactionClass() {};
	ReactionClass(const ReactionClass& orig) {};
	virtual ~ReactionClass() {};

public:



	double r(const double T, const double P_Pa, const std::vector<double>& c) const
	{
		r = ReactionRate.k(T, P_Pa,c) * ThirdBody.efficiency(c) * Reversibility.Product(c);
	}

protected:

	// Units
	PhysicalConstants::UNITS_REACTION a_units_;
	PhysicalConstants::UNITS_REACTION e_units_;
	
	// Stoichiometric coefficients
	std::vector<double> reactant_nu_;
	std::vector<double> product_nu_;
	std::vector<unsigned int>	reactant_nu_indices_;
	std::vector<unsigned int>	product_nu_indices_;

	// Exponent coefficients
	std::vector<double> reactant_lambda_;
	std::vector<unsigned int>	reactant_lambda_indices_;
	
	// Conversion factors
	double sumLambdaReactants_;
	double sumNuReactants_;
	double sumNuProducts_;
	double a_conversion_;
	double e_conversion_;
};

#endif	/* REACTIONCLASS_H */

