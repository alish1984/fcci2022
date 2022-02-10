/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_BatchReactor_Isothermal_ConstantVolume_H
#define	OpenSMOKE_BatchReactor_Isothermal_ConstantVolume_H

// Parent class
#include "BatchReactor.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "maps/SensitivityMap.h"

// Options
#include "BatchReactor_Options.h"
#include "ode/ODE_Parameters.h"


namespace OpenSMOKE
{
	//!  A class for simulating batch reactors with constant volume in isothermal conditions
	/*!
		 The purpose of this class is to simulate a batch reactor with constant volume, in isothermal conditions
		 The conservation equations of species are solved in terms of mass fractions.
	*/

	class BatchReactor_Isothermal_ConstantVolume : public BatchReactor
	{

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param kineticsMap map containing the kinetic mechanism
		*@param ode_parameters parameters governing the solution of the stiff ODE system
		*@param batch_options options governing the output
		*@param V0 initial volume [m3]
		*@param T0 initial temperature [K]
		*@param P0 initial pressure [Pa]
		*@param omega0 initial mass fractions of species
		*/
		BatchReactor_Isothermal_ConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
												OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap,
												OpenSMOKE::ODE_Parameters& ode_parameters,
												OpenSMOKE::BatchReactor_Options& batch_options,
												const double V0, const double T0, const double P0, 
												const OpenSMOKE::OpenSMOKEVectorDouble& omega0 );

		/**
		*@brief Solves the batch reactor
                *@param t0 the starting time of integration [s]
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double t0, const double tf);

		/**
		*@brief Ordinary differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

		/**
		*@brief Sparse Analystical Jacobian
		*@param t current time [s]
		*@param y current solution
		*@param J sparse analytical Jacobian
		*/
		virtual void SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J);

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	protected:

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&	thermodynamicsMap_;		//!< thermodynamic map
		OpenSMOKE::KineticsMap_CHEMKIN<double>&			kineticsMap_;			//!< kinetic map
		OpenSMOKE::ODE_Parameters&						ode_parameters_;		//!< ode parameters
		OpenSMOKE::BatchReactor_Options&				batch_options_;			//!< options
	};
}

#include "BatchReactor_Isothermal_ConstantVolume.hpp"

#endif	/* OpenSMOKE_BatchReactor_Isothermal_ConstantVolume_H */

