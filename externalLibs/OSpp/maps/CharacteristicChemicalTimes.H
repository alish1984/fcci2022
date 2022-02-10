/*----------------------------------------------------------------------------*\
|                                                                              |
|                                                                              |
|                    ┌─┐┌─┐┌─┐┬  ╦┌─┐╔═╗╔╦╗╔═╗╦╔═╔═╗                           |
|                    │  │ │├─┤│  ║└─┐╚═╗║║║║ ║╠╩╗║╣                            |
|                    └─┘└─┘┴ ┴┴─┘╩└─┘╚═╝╩ ╩╚═╝╩ ╩╚═╝                           |
|                                                                              |
|                                                                              |
|                                                                              |
|   Authors: A. Shamooni,                                                      |
|                                                                              |
|   Contacts: Ali Shamooni                                                     |
|   email: ali.shamooni@itv.uni-stuttgart.de                                   |
|   email: ali.shamooni@gmail.com                                              |
|   Institut für Technische Verbrennung (ITV),                                 |
|   Universität Stuttgart,                                                     |
|   Herdweg 51 (2nd floor), 70174  Stuttgart, Germany                          |
|                                                                              |
|------------------------------------------------------------------------------|
|                                                                              |
|   This file is part of coalIsSMOKE solver.                                   |
|                                                                              |
|	License                                                                    |
|                                                                              |
|   Copyright(C) 2019 A. Shamooni                                              |
|   coalIsSMOKE is free software: you can redistribute it and/or modify        |
|   it under the terms of the GNU General Public License as published by       |
|   the Free Software Foundation, either version 3 of the License, or          |
|   (at your option) any later version.                                        |
|                                                                              |
|   coalIsSMOKE is distributed in the hope that it will be useful,             |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|   GNU General Public License for more details.                               |
|                                                                              |
|   You should have received a copy of the GNU General Public License          |
|   along with coalIsSMOKE. If not, see <http://www.gnu.org/licenses/>.        |
|                                                                              |
\*----------------------------------------------------------------------------*/
#ifndef CharacteristicChemicalTimes_H
#define CharacteristicChemicalTimes_H

class CharacteristicChemicalTimes
{
	public:

	enum CharacteristicChemicalTimesType { CHEMICAL_TIMES_FORMATION_RATES, CHEMICAL_TIMES_REACTION_RATES, CHEMICAL_TIMES_EIGENVALUES };

	public:

	CharacteristicChemicalTimes(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&  		thermodynamicsMapXML,
					OpenSMOKE::KineticsMap_CHEMKIN<double>& 		kineticsMapXML);

	void SetThresholdChemicalTime(const double threshold);
	void SetThresholdTemperatureChemicalTime(const double threshold);
	void SetType(const CharacteristicChemicalTimesType type);
	void SetLapackMode(const bool flag);
	double CalculateCharacteristicChemicalTime(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble rf_;
	OpenSMOKE::OpenSMOKEVectorDouble tauChemSpecies_;

	OpenSMOKE::OpenSMOKEVectorDouble& getMoleFractions(){return  x_;}
	OpenSMOKE::OpenSMOKEVectorDouble& getFormationRates(){return  R_;}
	OpenSMOKE::OpenSMOKEVectorDouble& getTauChemSpecies(){return  tauChemSpecies_;}


	private:

	double FromFormationRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);
	double FromReactionRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);
	double FromEigenValueAnalysis(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

	void Equations(const OpenSMOKE::OpenSMOKEVectorDouble& c, OpenSMOKE::OpenSMOKEVectorDouble& dc_over_dt);

	void NumericalJacobian(const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd& J);

	void SearchConservativeModes(const std::vector<double>& lambda_real, const std::vector<double>& lambda_mod, std::vector<double>& lambda_real_cleaned);

	private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN<double>& 	kineticsMapXML_;

	unsigned int ne_;
	unsigned int ns_;
	unsigned int nr_;
	unsigned int n_;

	bool lapack_mode_;

	double T_;
	double P_;

	Eigen::MatrixXd Jc_;
	Eigen::VectorXd sum_nu_;

	Eigen::VectorXd vr_;
	Eigen::VectorXd vl_;



	// Local variables
	std::vector<double> lambda_real_;
	std::vector<double> lambda_imag_;
	std::vector<double> lambda_mod_;
	std::vector<double> lambda_real_cleaned_;

	double tauChem_threshold_;
	double tauChem_temperature_threshold_;

	CharacteristicChemicalTimesType type_;
};

#include "CharacteristicChemicalTimes.C"

#endif // CharacteristicChemicalTimes_H
