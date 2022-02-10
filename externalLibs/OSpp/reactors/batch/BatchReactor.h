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

#ifndef OpenSMOKE_BatchReactor_H
#define	OpenSMOKE_BatchReactor_H

#include "math/OpenSMOKEVector.h"
#include "BatchReactor_Options.h"
#include "reactors/utilities/SensitivityAnalysis_Options.h"
#include "maps/SensitivityMap.h"
#include "boost/filesystem.hpp"
#include "ode/ODE_Parameters.h"

namespace OpenSMOKE
{
	enum BatchReactor_Type { BATCH_REACTOR_NONISOTHERMAL_CONSTANTP, BATCH_REACTOR_NONISOTHERMAL_CONSTANTV,
							 BATCH_REACTOR_ISOTHERMAL_CONSTANTP, BATCH_REACTOR_ISOTHERMAL_CONSTANTV,
							 BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME };

	//!  A class for simulating batch reactors
	/*!
		 The purpose of this class is to simulate batch reactors. It provides a common interface to different
		 batch reactors
	*/

	class BatchReactor
	{

	public:

		/**
		*@brief Solves the batch reactor
		*@param tf the final time of integration [s]
		*/
		virtual void Solve(const double t0, const double tf) = 0;

		/**
		*@brief Ordinary differential Equations corresponding to the reactor under investigation
		*@param t current time [s]
		*@param y current solution
		*@param dy current unsteady terms
		*/
		virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;

		/**
		*@brief Sparse Analytical Jacobian
		*@param t current time [s]
		*@param y current solution
		*@param J sparse analytical Jacobian matrix
		*/
		virtual void SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J) = 0;

		/**
		*@brief Writes the output (called at the end of each time step)
		*@param t current time [s]
		*@param y current solution
		*/
		virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param t current time [s]
		*@param y current solution
		*/
		virtual void SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y) = 0;

		/**
		*@brief Return the total number of equations
		*@return the total number of equations
		*/
		unsigned int NumberOfEquations() const { return NE_; };
		
		/**
		*@brief Performs the sensitivity analysis at a specific time
		*@param sensitivityMap map of sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param sensitivity_options options for sensitivity analysis
		*/		
		void EnableSensitivityAnalysis(OpenSMOKE::SensitivityMap<double>& sensitivityMap, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, BatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
	
		/**
		*@brief Returns the current value of temperature [K]
		*/	
		double T() const { return T_; }

		/**
		*@brief Returns the current value of pressure [Pa]
		*/	
		double P_Pa() const { return P_; }

		/**
		*@brief Returns the current value of volume [m3]
		*/	
		double V() const { return V_; }

		/**
		*@brief Returns the current values of mass fractions [-]
		*/	
		const OpenSMOKE::OpenSMOKEVectorDouble& omega() const { return omega_; }
		
		/**
		*@brief Writes on the video the final status of the batch reactor
		*@param t0 start time of integration (usually 0)
        *@param tf final time of integration
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options governing the solution
		*/			
		void FinalStatus(const double t0, const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::BatchReactor_Options& batch_options);

		/**
		*@brief Writes on a file the final summary of the batch reactor
		*@param summary_file file where the summary will be written
        *@param t0 start time of integration (usually 0)
        *@param tf end time of integration
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options governing the solution
		*/		
		void FinalSummary(const boost::filesystem::path summary_file, const double t0, const double tf, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::BatchReactor_Options& batch_options);

		/**
		*@brief Returns the final composition of the batch reactor, together with temperature and pressure
		*@param T temperature [K]
		*@param P pressure [Pa]
		*@param omega mass fractions [-]
		*/	
		void GetFinalStatus(double& T, double& P, OpenSMOKE::OpenSMOKEVectorDouble& omega);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param psr_options options for solving the psr
		*/
		void PrepareASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, BatchReactor_Options& batch_options);

		/**
		*@brief Prepares the ASCII file
		*@param fOutput ASCII file where the output will be written
		*@param output_file_ascii name of ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param psr_options options for solving the psr
		*/
		void PrepareParametricASCIIFile(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, BatchReactor_Options& batch_options);


		/**
		*@brief Print the final status on output stream
		*@param fOutput output stream
		*@param t time
		*/
		void PrintFinalStatus(std::ostream& fOutput, const double t);

		/**
		*@brief Print the final status on output stream
		*@param fOutput output stream
		*@param t time
		*/
		void PrintParametricFinalStatus(std::ostream& fOutput, const double t);

	protected:

		/**
		*@brief Prepares the ASCII file
		*@param output_file_ascii ASCII file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareASCIIFile(const boost::filesystem::path output_file_ascii, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, BatchReactor_Options& batch_options);
	
		/**
		*@brief Prepares the XML file
		*@param output_file_xml XML file where the output will be written
		*@param thermodynamicsMap map of thermodynamic data
		*/
		void PrepareXMLFile(const boost::filesystem::path output_file_xml, OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap);
		
		/**
		*@brief Closes the XML file
		*/		
		void CloseXMLFile();

		/**
		*@brief Calculates the numerical Jacobian
		*@param t current time
		*@param y current solution
		*@param J Jacobian matrix
		*/
		void NumericalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEMatrixDouble& J);	
		
		/**
		*@brief Allocates the memory for the vectors and the matrices contained in the batch reactor
		*/			
		void MemoryAllocation();

		/**
		*@brief Prepares the XML files for the sensitivity analysis
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void PrepareSensitivityXMLFiles(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, BatchReactor_Options& batch_options, OpenSMOKE::SensitivityAnalysis_Options& sensitivity_options);
		
		/**
		*@brief Closes the XML files for the sensitivity analysis
		*/			
		void CloseSensitivityXMLFiles();

		/**
		*@brief Open all the files (ASCII and XML)
		*@param thermodynamicsMap map of thermodynamic data
		*@param batch_options options for solving the batch reactor
		*/
		void OpenAllFiles(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::BatchReactor_Options& batch_options);
		
		/**
		*@brief Closes all the files (ASCII and XML)
		*@param batch_options options for solving the batch reactor
		*/		
		void CloseAllFiles(OpenSMOKE::BatchReactor_Options& batch_options);

		/**
		*@brief Solves the Ordinary Differential Equations using one of the provided open/source ODE solvers
		*@param t0 start time of integrations (usually 0)
                *@param tf final time of integrations
		*@param ode_parameters parameters governing the ODE integration
		*/	
		void SolveOpenSourceSolvers(const double t0, const double tf, OpenSMOKE::ODE_Parameters& ode_parameters);
	
protected:

		BatchReactor_Type type_;					//!< type of batch reactor

		double T0_;									//!< initial temperature [K]
		double P0_;									//!< initial pressure [Pa]
		double V0_;									//!< initial volume [m3]
		double MW0_;								//!< initial molecular weight [kg/kmol]
		double H0_;									//!< initial enthalpy [J/kmol]
		double U0_;									//!< initial internal energy [J/kmol]
		OpenSMOKE::OpenSMOKEVectorDouble omega0_;	//!< initial composition [mass fractions]
		OpenSMOKE::OpenSMOKEVectorDouble x0_;		//!< initial composition [mole fractions]

		double rho_;	//!< density [kg/m3]
		double mass_;	//!< mass [kg]
                
		double global_thermal_exchange_coefficient_;	//!< global thermal exchange coefficient [W/m2/K]
		double exchange_area_;				//!< exchange area [m2]
		double T_environment_;				//!< environment temperature                

		unsigned int NC_;	//!< number of species [-]
		unsigned int NR_;	//!< number of reactions [-]
		unsigned int NE_;	//!< number of equations [-]

		OpenSMOKE::OpenSMOKEVectorDouble omega_;	//!< current mass fractions
		OpenSMOKE::OpenSMOKEVectorDouble x_;		//!< current mole fractions
		OpenSMOKE::OpenSMOKEVectorDouble c_;		//!< current concentrations [kmol/m3]
		OpenSMOKE::OpenSMOKEVectorDouble R_;		//!< current formation rates [kmol/m3/s]
		
		double T_;									//!< current temperature [K]
		double P_;									//!< current pressure [Pa]
		double V_;									//!< current volume [V]
		double cTot_;								//!< current concentration [kmol/m3]
		double MW_;									//!< current molecular weight [kg/kmol]
		double QR_;									//!< current heat release [W/m3]
		double CvMixMass_;							//!< current specific heat (constant volume) [J/kg/K]
		double CpMixMass_;							//!< current specific heat (constant pressure) [J/kg/K]

		unsigned int iteration_;					//!< iteration index
		unsigned int counter_file_video_;			//!< iteration counter for writing on video
		unsigned int counter_file_ASCII_;			//!< iteration counter for writing on ASCII file
		unsigned int counter_file_XML_;				//!< iteration counter for writing on XML file
		unsigned int counter_sensitivity_XML_;		//!< iteration counter for writing on XML sensitivity files

		OpenSMOKE::OpenSMOKEVectorDouble y0_;		//!< vector containing the initial values for all the variables
		OpenSMOKE::OpenSMOKEVectorDouble yf_;		//!< vector containing the final values for all the variables

		std::ofstream fXML_;						//!< XML file where the output is written
		std::ofstream fASCII_;						//!< ASCII file where the output is written
		std::ofstream  fSensitivityParentXML_;		//!< XML file where the sensitivity analysis is written (parent file)
		std::ofstream* fSensitivityChildXML_;		//!< XML files where the sensitivity analysis is written (children files)

		std::vector<unsigned int> indices_of_output_species_;		//!< indices of species for which the profiles are written on the ASCII file
		std::vector<unsigned int> indices_of_sensitivity_species_;	//!< indices of species for which the sensitivity analysis is written on file
		std::vector<unsigned int> widths_of_output_species_;        //!< width   of species for which the profiles are written on the ASCII file
		
		OpenSMOKE::SensitivityMap<double>* sensitivityMap_;			//!< sensitivity map
		OpenSMOKE::OpenSMOKEVectorDouble scaling_Jp_;				//!< vector containing datafor performing the sensitivity analysis
		OpenSMOKE::OpenSMOKEMatrixDouble Jnum_;						//!< Jacobian matrix (required for sensitivity analysis)
		Eigen::SparseMatrix<double> Jan_;							//!< Jacobian matrix (required for sensitivity analysis)
	};
}

#include "BatchReactor.hpp"

#endif	/* OpenSMOKE_BatchReactor_H */

