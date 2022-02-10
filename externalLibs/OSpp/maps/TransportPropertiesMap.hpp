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

namespace OpenSMOKE
{
	template<typename map> 
	void TransportPropertiesMap<map>::ThermalConductivity(map& lambdamix, OpenSMOKEVectorDouble& moleFractions)
	{
		lambda();
		lambdaMix(lambdamix, moleFractions);
	}


	template<typename map> 
	void TransportPropertiesMap<map>::DynamicViscosity(map& etamix, OpenSMOKEVectorDouble& moleFractions)
	{
		eta();
		etaMix(etamix, moleFractions);
	}

	template<typename map> 
	void TransportPropertiesMap<map>::MassDiffusionCoefficients(OpenSMOKEVectorDouble& gammamix, OpenSMOKEVectorDouble& moleFractions, const bool bundling)
	{
		if (bundling == false)
		{
			gamma();
			gammaMix(gammamix, moleFractions);
		}
		else
		{
			bundling_gamma();
			bundling_gammaMix(gammamix, moleFractions);
		}
	}

	template<typename map> 
	void TransportPropertiesMap<map>::ThermalDiffusionRatios(OpenSMOKEVectorDouble& tetamix, OpenSMOKEVectorDouble& moleFractions)
	{
		teta();
		tetaMix(tetamix, moleFractions);
	}
}
