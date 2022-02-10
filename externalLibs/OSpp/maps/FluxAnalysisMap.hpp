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

#include "math/OpenSMOKEUtilities.h"
#include <Eigen/Sparse>
#include <iomanip>

namespace OpenSMOKE
{
	template<typename map>
	void FluxAnalysisMap<map>::AnalyzeNetFluxes(	const unsigned int index_j,
													std::vector<unsigned int>& important_indices,
													std::vector<double>& important_normal_fluxes,
													std::vector<double>& important_fluxes)
	{
		std::cout << "/* Entering AnalyzeNetFluxes */" << '\n';
		std::cout << "/* index j */" <<index_j<< '\n';
		std::cout << "/* index_atomN */" <<index_atom_ <<'\n';
		const double n_j = thermodynamicsMapXML_.atomic_composition()(index_j,index_atom_);

		std::cout << "/* species j */" << thermodynamicsMapXML_.NamesOfSpecies()[index_j]<<'\n';
		std::cout << "/* n_j */" << n_j<<'\n';

		std::vector<double> number_of_atoms(kineticsMapXML_.NumberOfReactions());
		std::vector<std::vector<unsigned int> > reactions_of_single_reactant(NC);
		std::vector<std::vector<unsigned int> > reactions_of_single_product(NC);

		std::cout << "/* before  stoichiometry map AnalyzeNetFluxes */" << '\n';
		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants(),k); it; ++it)
			{
				reactions_of_single_reactant[k].push_back(it.row());
				number_of_atoms[it.row()] += thermodynamicsMapXML_.atomic_composition()(k,index_atom_);
			}

		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_products(),k); it; ++it)
				reactions_of_single_product[k].push_back(it.row());

std::cout << "/* Species j as reactant  */" << '\n';
		// Species j as reactant
		std::vector<double> fluxes_as_reactant(NC);
		{
			for (unsigned int k=0; k<NC; ++k)
				fluxes_as_reactant[k] = 0.;

			for (unsigned int k=0; k<NC; ++k)
			{
				const double n_k = thermodynamicsMapXML_.atomic_composition()(k,index_atom_);

				if (n_k != 0 && index_j != k)
				{
					// Recognize common reactions
					std::vector<unsigned int> couples;
					for (unsigned int kk=0; kk<reactions_of_single_product[k].size(); ++kk)
						for (unsigned int jj=0; jj<reactions_of_single_reactant[index_j].size(); ++jj)
							if (reactions_of_single_product[k][kk] == reactions_of_single_reactant[index_j][jj])
							{
								couples.push_back(reactions_of_single_product[k][kk]);
								break;
							}

					if (destruction_analysis_ == true)
					{
						for (unsigned int j=0;j<couples.size();j++)
							if (r_[couples[j]+1]>0.) fluxes_as_reactant[k] += n_k*n_j * r_[couples[j]+1]/number_of_atoms[couples[j]];
					}
					else
					{
						for (unsigned int j=0;j<couples.size();j++)
							if (r_[couples[j]+1]<0.) fluxes_as_reactant[k] += n_k*n_j * r_[couples[j]+1]/number_of_atoms[couples[j]];
					}
				}
			}
		}
		std::cout << "/* Species j as product  */" << '\n';

		// Species j as product
		std::vector<double> fluxes_as_product(NC);
		{
			for (unsigned int k=0; k<NC; ++k)
				fluxes_as_product[k] = 0.;

			for (unsigned int k=0; k<NC; ++k)
			{
				const double n_k = thermodynamicsMapXML_.atomic_composition()(k,index_atom_);

				if (n_k != 0 && index_j != k)
				{
					// Recognize common reactions
					std::vector<unsigned int> couples;
					for (unsigned int kk=0; kk<reactions_of_single_reactant[k].size(); ++kk)
						for (unsigned int jj=0; jj<reactions_of_single_product[index_j].size(); ++jj)
							if (reactions_of_single_reactant[k][kk] == reactions_of_single_product[index_j][jj])
							{
								couples.push_back(reactions_of_single_reactant[k][kk]);
								break;
							}

					if (destruction_analysis_ == true)
					{
						for (unsigned int j=0;j<couples.size();j++)
							if (r_[couples[j]+1]<0.) fluxes_as_product[k] += n_k*n_j * r_[couples[j]+1]/number_of_atoms[couples[j]];
					}
					else
					{
						for (unsigned int j=0;j<couples.size();j++)
							if (r_[couples[j]+1]>0.) fluxes_as_product[k] += n_k*n_j * r_[couples[j]+1]/number_of_atoms[couples[j]];
					}
				}
			}
		}
std::cout << "/* fluxes */" << '\n';
		std::vector<double> fluxes(NC);
		for (unsigned int j=0;j<NC;j++)
			fluxes[j] = std::fabs(fluxes_as_reactant[j]+fluxes_as_product[j]);

std::cout << "/* indices */" << '\n';

		std::vector<unsigned int> indices(NC);
		for (unsigned int j=0;j<NC;j++)
			indices[j] = j;

std::cout << "/* reorder */" << '\n';

		OpenSMOKE_Utilities::ReorderPairsOfVectors(fluxes, indices);
		std::reverse(fluxes.begin(), fluxes.end());
		std::reverse(indices.begin(), indices.end());

		// Select only important fluxes
		{
			double sum = 0.;
			for (unsigned int j=0;j<NC;j++)
				sum += fluxes[j];

			std::vector<double> normal_fluxes(NC);
			for (unsigned int j=0;j<NC;j++)
				normal_fluxes[j] = fluxes[j] / sum * 100.;

			int index_max = 0;
			for (unsigned int j=0;j<NC;j++)
				if (normal_fluxes[j] < min_percentage_threshold_)
				{
					index_max = j;
					break;
				}
std::cout << "/* important_indices */" << '\n';
			important_indices.resize(std::min(max_width_,index_max));
			important_normal_fluxes.resize(std::min(max_width_,index_max));
			important_fluxes.resize(std::min(max_width_,index_max));
			for (int j=0;j<std::min(max_width_,index_max);j++)
			{
				important_indices[j] = indices[j];
				important_normal_fluxes[j] = normal_fluxes[j];
				important_fluxes[j] = normal_fluxes[j]/100. * sum;
			}
		}
	}

	template<typename map>
	void FluxAnalysisMap<map>::OpenGraphFile(const std::string file_name)
	{
		fOut.open(file_name.c_str(), std::ios::out);
		fOut << "digraph myFirst {" << std::endl;
		fOut << "node [shape=circle];  " << std::endl;
	}

	template<typename map>
	void FluxAnalysisMap<map>::CloseGraphFile()
	{
		fOut << "overlap=false" << std::endl;
		fOut << "label=\"Reaction Flux Analysis from OpenSMOKE++ Suite (www.opensmoke.polimi.it)\"" << std::endl;
		fOut << "fontsize=14;" << std::endl;
		fOut << "}" << std::endl;
		fOut.close();
	}

	template<typename map>
	void FluxAnalysisMap<map>::AddSpeciesToGraphFile(	const unsigned int index_j,
														std::vector<unsigned int>& local_indices,
														std::vector<double>& local_thickness,
														std::vector<double>& local_normal_fluxes,
														std::vector<double>& local_fluxes)
	{

		for(unsigned int j=0;j<local_indices.size();j++)
		{
			std::string attributes;

			// Tags
			if (normal_tags_ == true)
			{
				std::stringstream label;
				label.setf(std::ios::fixed, std::ios::floatfield);

				if (local_normal_fluxes[j]>10.)		label.precision(1);
				else if (local_normal_fluxes[j]>1.)	label.precision(2);
				else if (local_normal_fluxes[j]>0.1)	label.precision(3);
				else label.precision(3);
				label << local_normal_fluxes[j];
				attributes = "[label = \"" + label.str() + "%\", labelfontcolor = red];";
			}
			else
			{
				std::stringstream label;
				label.setf(std::ios::scientific);
				label.precision(2);
				label << local_fluxes[j];
				attributes = "[label = \"" + label.str() + "\", labelfontcolor = red];";
			}

			// Thickness
			std::stringstream thickness;
			if (logarithmic_thickness_ == false)	thickness << 0.5 + 15.*local_thickness[j];
			else									thickness << 0.5 + 5.*std::max(0.,3.+log10(local_thickness[j]));

			// Write


			if (destruction_analysis_ == true)
			{
				std::string name1 = thermodynamicsMapXML_.NamesOfSpecies()[index_j];
				std::string name2 = thermodynamicsMapXML_.NamesOfSpecies()[local_indices[j]];
				boost::replace_all(name1, "-", "_");
				boost::replace_all(name2, "-", "_");
				fOut << "edge [color=red, penwidth = " + thickness.str() << "];" << std::endl;
				fOut << name1 << "->" << name2 << " " << attributes << std::endl;
			}
			else
			{
				std::string name1 = thermodynamicsMapXML_.NamesOfSpecies()[local_indices[j]];
				std::string name2 = thermodynamicsMapXML_.NamesOfSpecies()[index_j];
				boost::replace_all(name1, "-", "_");
				boost::replace_all(name2, "-", "_");
				fOut << "edge [color=blue, penwidth = " + thickness.str() << "];" << std::endl;
				fOut << name1 << "->" << name2 << " " << attributes << std::endl;
			}
		}
	}

	template<typename map>
	void FluxAnalysisMap<map>::GloballyAnalyze(	const std::vector<unsigned int>& important_indices, const int current_depth )
	{
		if (current_depth <= max_depth_)
		{
			std::cout << "/* Entering  GloballyAnalyze in depth */" <<current_depth<< '\n';
			std::cout << "/* list_of_analyzed_species_ size */" <<list_of_analyzed_species_.size()<< '\n';
			for (unsigned int AUXI=0;AUXI<list_of_analyzed_species_.size();AUXI++)
			std::cout << "/* list_of_analyzed_species_ */" <<list_of_analyzed_species_[AUXI]<< '\n';

			std::cout << "/* important_indices size */" <<important_indices.size()<< '\n';
			for (unsigned int AUXI=0;AUXI<important_indices.size();AUXI++)
			std::cout << "/* important_indices */" <<important_indices[AUXI]<< '\n';

			for (unsigned int j=0;j<important_indices.size();j++)
			{
				const unsigned int index_jj = important_indices[j];
				std::cout << "/* important_indices j */" <<important_indices[j]<< '\n';

				if (std::find(list_of_analyzed_species_.begin(), list_of_analyzed_species_.end(), index_jj) == list_of_analyzed_species_.end())
				{
					AnalyzeNetFluxes(index_jj, global_important_indices_[index_jj], global_important_normal_fluxes_[index_jj], global_important_fluxes_[index_jj]);
					list_of_analyzed_species_.push_back(index_jj);
					GloballyAnalyze(global_important_indices_[index_jj], current_depth+1);
				}
			}
		}
	}

	template<typename map>
	void FluxAnalysisMap<map>::CalculateThickness()
	{
		if (normal_thickness_ == false)
		{
			double max_flux = 0.;
			for (unsigned int j=0;j<list_of_analyzed_species_.size();j++)
			{
				const unsigned int index_j = list_of_analyzed_species_[j];
				for (unsigned int k=0;k<global_important_indices_[index_j].size();k++)
					if ( global_important_fluxes_[index_j][k] > max_flux) max_flux = global_important_fluxes_[index_j][k];
			}

			for (unsigned int j=0;j<list_of_analyzed_species_.size();j++)
			{
				const unsigned int index_j = list_of_analyzed_species_[j];
				global_relative_thickness_[index_j].resize(global_important_indices_[index_j].size());
				for (unsigned int k=0;k<global_important_indices_[index_j].size();k++)
					global_relative_thickness_[index_j][k] = global_important_fluxes_[index_j][k]/max_flux;
			}
		}

		if (normal_thickness_ == true)
		{
			for (unsigned int j=0;j<list_of_analyzed_species_.size();j++)
			{
				const unsigned int index_j = list_of_analyzed_species_[j];
				global_relative_thickness_[index_j].resize(global_important_indices_[index_j].size());
				for (unsigned int k=0;k<global_important_indices_[index_j].size();k++)
					global_relative_thickness_[index_j][k] = global_important_normal_fluxes_[index_j][k]/100.;
			}
		}
	}

	template<typename map>
	void FluxAnalysisMap<map>::Plot(const std::string file_name)
	{
		OpenGraphFile(file_name);
		for (unsigned int j=0;j<list_of_analyzed_species_.size();j++)
		{
			const unsigned int index_j = list_of_analyzed_species_[j];
			AddSpeciesToGraphFile(index_j,	global_important_indices_[index_j], global_relative_thickness_[index_j],
											global_important_normal_fluxes_[index_j], global_important_fluxes_[index_j]);
		}
		CloseGraphFile();
	}

	template<typename map>
	void FluxAnalysisMap<map>::WriteStoichiometricMatrix( const std::string file_name)
	{
		std::ofstream fOut(file_name.c_str(), std::ios::out);

		fOut << NC << std::endl;
		fOut << NR << std::endl;

		for (unsigned int j=0;j<NC;j++)
			fOut << thermodynamicsMapXML_.NamesOfSpecies()[j] << std::endl;

		unsigned int number_of_reactant_entries = 0;
		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants(),k); it; ++it)
				number_of_reactant_entries++;

		fOut << number_of_reactant_entries << std::endl;
		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants(),k); it; ++it)
				fOut << k+1 << " " << it.row()+1 << " " << it.value() << std::endl;

		unsigned int number_of_product_entries = 0;
		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_products(),k); it; ++it)
				number_of_product_entries++;

		fOut << number_of_product_entries << std::endl;
		for (int k=0; k<kineticsMapXML_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_products(),k); it; ++it)
				fOut << k+1 << " " << it.row()+1 << " " << it.value() << std::endl;

		fOut.close();
	}
}
