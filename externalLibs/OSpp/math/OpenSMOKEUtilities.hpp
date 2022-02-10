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
|   License                                                               |
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

#ifndef OpenSMOKE_OpenSMOKEUtilities_Cpp
#define OpenSMOKE_OpenSMOKEUtilities_Cpp

//#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEFunctions.h"
#include <boost/date_time.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace OpenSMOKE
{
	template <typename T>
	void Swap(T* x, T *y)
	{
		T temp = *x;
		*x = *y;
		*y = temp;
	}


	template <typename T>
	inline T Abs(T const& a)
	{
	   return a < 0 ? -a : a;
	}


	template <typename T>
	inline T const& Max(T const& a, T const& b)
	{
		return (a > b ? a : b);
	}


	template <typename T>
	inline T const& Min(T const& a, T const& b)
	{
		return (a < b ? a : b);
	}


	template <typename T>
	inline T MaxAbs(T const& a, T const& b)
	{
		return ( Abs(a) > Abs(b) ? Abs(a) : Abs(b));
	}

	template <typename T>
	inline T MinAbs(T const& a, T const& b)
	{
		return ( Abs(a) < Abs(b) ? Abs(a) : Abs(b) );
	}

	template <typename T>
	inline T Max(const int n, T *x, int *im)
	{
		if(n < 0)
			return x[0];

		T temp;
		temp = x[0];
		if(im != 0)
			*im = 0;

		for(int i=1;i<n;i++)
			if(temp<x[i])
			{
				temp = x[i];
				if(im != 0)
					*im = i;
			}

		return temp;
	}

	template <class T>
	inline T Max(const int n, T *x)
	{
		if(n < 0)
			return x[0];

		T temp;
		temp = x[0];
		for(int i=1;i<n;i++)
			if(temp < x[i])
				temp = x[i];

		return temp;
	}

	template <class T>
	inline T MaxAbs(const int n, T *x,int *im)
	{
		T temp = Abs(x[0]);
		if(n < 0)
			return temp;

		if(im != 0)
			*im = 0;

		for(int i=1;i<n;i++)
			if(temp < Abs(x[i]))
			{
				temp = Abs(x[i]);
				if(im != 0)
				*im = i;
			}

		return temp;
	}

	template <class T>
	inline T MaxAbs(const int n, T *x)
	{
		T temp = Abs(x[0]);
		if(n < 0)
			return temp;

		for(int i=1;i<n;i++)
			if(temp < Abs(x[i]))
				temp = Abs(x[i]);

		return temp;
	}

	template <class T>
	inline T Min(const int n, T *x, int *im)
	{
		if(n < 0)
			return x[0];

		T temp = x[0];
		if(im != 0)
			*im = 0;

		for(int i = 1;i < n;i++)
			if(temp > x[i])
			{
				temp = x[i];
				if(im != 0)
					*im = i;
			}

		return temp;
	}

	template <class T>
	inline T Min(const int n, T *x)
	{
		if(n < 0)
			return x[0];

		T temp = x[0];
		for(int i=1;i<n;i++)
			if(temp > x[i])
				temp = x[i];

		return temp;
	}

	template <class T>
	inline T MinAbs(const int n, T *x, int *im)
	{
		T temp = Abs(x[0]);
		if(n < 0)
			return temp;
		if(im != 0)
			*im = 0;
		for(int i = 1;i < n;i++)
			if(temp > Abs(x[i]))
			{
				temp = Abs(x[i]);
				if(im != 0)
				*im = i;
			}
		return temp;
	}

	template <class T>
	inline T MinAbs(const int n, T *x)
	{
		T temp = Abs(x[0]);

		if(n < 0)
			return temp;

		for(int i = 1;i < n;i++)
			if(temp > Abs(x[i]))
				temp = Abs(x[i]);

		return temp;
	}

	template <class T>
	void Sum(const int n, T* lval, T* rval, T* result)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lval++) + (*rval++);
			*result++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lval, const T rval, T* result)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lval++) + rval;
			*result++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, const T rval, T* lvalAndResult)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalAndResult) + rval;
			*lvalAndResult++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lvalAndResult, T* rval)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalAndResult) + (*rval++);
			*lvalAndResult++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lvalRvalAndResult)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalRvalAndResult) + (*lvalRvalAndResult);
			*lvalRvalAndResult++ = sum;
		}
	}

	template <class T>
	void Difference(const int n, T* lval, T* rval, T* result)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lval++) - (*rval++);
			*result++ = diff;
		}
	}

	template <class T>
	void Difference(const int n, T* lvalAndResult, T* rval)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lvalAndResult) - (*rval++);
			*lvalAndResult++ = diff;
		}
	}

	template <class T>
	void DifferenceBis(const int n, T* lval, T* rvalAndResult)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lval++) - (*rvalAndResult);
			*rvalAndResult++ = diff;
		}
	}

	template <class T>
	T Dot(const int n, T* lval, T* rval)
	{
		T result = 0;
		for(int i=0;i<n;i++)
			result += (*lval++) * (*rval++);
		return result;
	}

	template <class T>
	T UDot(const int n, T* lval, T* rval)
	{
		T result = 0;
		for(int i=0;i<n;i++)
			result += (*lval++) / (*rval++);
		return result;
	}

	template <class T>
	void Prod(const int n, T lval, T* rval, T* result)
	{
		for(int i=0;i<n;i++)
			*result++ = lval * (*rval++);
	}

	template <class T>
	void Prod(const int n, T lval, T* rvalAndResult)
	{
		for(int i=0;i<n;i++)
		{
			*rvalAndResult = (lval) * (*rvalAndResult);
			 rvalAndResult++;
		}
	}

	template <class T>
	void Div(const int n, T* lval, T rval, T* result)
	{
		if(rval == 0)
			ErrorMessage("Div(const int n, T* lval, T rval, T* result)", "Division with 0");

		for(int i=0;i<n;i++)
			*result++ = (*lval++)/double(rval);
	}

	template <class T>
	void Div(const int n, T* lvalAndResult, T rval)
	{
		if(rval == 0)
			ErrorMessage("Div(const int n, T* lvalAndResult, T rval)", "Division with 0");
		for(int i=0;i<n;i++)
		{
			*lvalAndResult = (*lvalAndResult)/double(rval);
			 lvalAndResult++;
		}
	}

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, T& value)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fInput.read(reinterpret_cast < char * > (&value), sizeof(T)))
				ErrorMessage(	"Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, T& value)",
								"I was unable to read an std::string coefficient from binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fInput >> value;
		}
	}

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, T& value)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fOutput.write( reinterpret_cast < char * > (&value), sizeof(T)))
				ErrorMessage(	"Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, T const& value)",
								"I was unable to write on a binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << value << std::endl;
		}
	}

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			for(int i=0;i<n;i++)
				if(!fInput.read(reinterpret_cast < char * > (&values[i]), sizeof(T)))
					ErrorMessage(	"Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)",
									"I was unable to read from binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			for(int i=0;i<n;i++)
				fInput >> values[i];
		}
	}

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			for(int i=0;i<n;i++)
				if(!fOutput.write( reinterpret_cast < char * > (&values[i]), sizeof(T)))
					ErrorMessage(	"Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, const int n, const T* values)",
									"I was unable to write on a binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			for(int i=0;i<n;i++)
				fOutput << values[i] << std::endl;
		}
	}

	template <typename T>
	inline unsigned char Compare_LE(const T a, const T b)
	{
		return (a <= b);
	}

	template <typename T>
	inline unsigned char Compare_LT(const T a, const T b)
	{
		return (a < b);
	}

	template <typename T>
	inline unsigned char Compare_GE(const T a, const T b)
	{
		return (a >= b);
	}

	template <typename T>
	inline unsigned char Compare_GT(const T a, const T b)
	{
		return (a > b);
	}

	template <typename T>
	void Sort(const int n, T *x, int *iS)
	{
		if(n <= 1)
			return;

		int j,k,ik,jk;

		for(int node=1;node<n;node++)
		{
			int i = node;
			j = ((i + 1)/2) - 1;
			while(i != 0 && Compare_LE(x[j],x[i]))
			{
				Swap(x + j,x + i);
				Swap(iS + j,iS + i);
				i = j;
				j = ((i + 1)/2) - 1;
			}
		}

		for(int i=n-1;i>=1;i--)
		{
			Swap(x + i,x);
			Swap(iS + i,iS);
			k = i - 1;
			ik = 0;
			jk = 1;

			if(k >= 2 && Compare_GT(x[2],x[1]))
				jk = 2;

			while(jk <= k && Compare_GT(x[jk],x[ik]))
			{
				Swap(iS + jk,iS + ik);
				Swap(x + jk,x + ik);
				ik = jk;
				jk = (2*(ik + 1)) - 1;
				if(jk + 1 <= k)
					if(Compare_GT(x[jk + 1],x[jk]))
						jk++;
			}
		}
	}

	template <typename T>
	std::vector<size_t> sort_and_track_indices(std::vector<T> const& values)
	{
		std::vector<size_t> indices(values.size());

		#if defined(_WIN32) || defined(_WIN64)

		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort(begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] > values[b]; });

		#else

		struct IdxCompare
		{
			const std::vector<T>& target;
			IdxCompare(const std::vector<T>& target) : target(target) {}
			bool operator()(size_t a, size_t b) const { return target[a] > target[b]; }
		};

		for (size_t i = 0; i < indices.size(); ++i)
			indices[i] = i;

		std::sort(indices.begin(), indices.end(), IdxCompare(values));

		#endif

		return indices;
	}

    template <typename T>
	std::vector<size_t> sort_and_track_indices_increasing(std::vector<T> const& values)
	{
		std::vector<size_t> indices(values.size());



		struct IdxCompare
		{
			const std::vector<T>& target;
			IdxCompare(const std::vector<T>& target) : target(target) {}
			bool operator()(size_t a, size_t b) const { return target[a] < target[b]; }
		};

		for (size_t i = 0; i < indices.size(); ++i)
			indices[i] = i;

		std::sort(indices.begin(), indices.end(), IdxCompare(values));



		return indices;
	}

	void EM(const std::string msg)
	{
		std::cout << msg << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	template <typename T>
	void Save(T value, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
				fOutput << value << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fOutput.write( reinterpret_cast < char * > (&value), sizeof(T)))
				EM("I was unable to write on binary file");
		}
	}

	template <typename T>
	void Load(T* value, std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
				fInput >> *value;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fInput.read(reinterpret_cast<char *>(value), sizeof(T)))
				EM("I was unable to read from binary file");
		}
	}

				// Reading vector size


	template<typename T>
	void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << vector.size() << std::endl;
			for(unsigned int i=0;i<vector.size();i++)
				fOutput << vector[i] << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			// Writing vector size
			if(!fOutput.write( reinterpret_cast < char * > (&vector.size()), sizeof(vector.size())))
				ErrorMessage("void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)", "I was unable to write on binary file");

			// Writing vector elements
			for(unsigned int i=0;i<vector.size();i++)
			if(!fOutput.write( reinterpret_cast < char * > (&vector[i]), sizeof(T)))
				ErrorMessage("void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)", "I was unable to write on binary file");
		}
	}

	template<typename T>
	void CheckIfFileIsOpen(const T& fOut, const std::string file_name)
	{
		if (!fOut.is_open())
		{
			std::cout << " ! Fatal error: I was unable to open this file: " << std::endl;
			std::cout << "   " << file_name << std::endl;
			std::cout << "   Press enter to exit..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}
	}

	void CreateDirectory(const boost::filesystem::path& path_output)
	{
		boost::filesystem::create_directories(path_output);

		if (!boost::filesystem::exists(path_output))
		{
			std::cout << "Error while attempting to create the " << path_output.string() << " directory." << std::endl;
			std::cout << "Please check the permissions you have on the folder." << std::endl;
			std::cout << "Press enter to exit..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}
	}

	std::string GetCurrentTime()
	{
		std::string current_time =  boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());
		current_time.erase(0,9);
		current_time.insert(2,":");
		current_time.insert(5,":");

		return current_time;
	}

	std::string GetRawCurrentTime()
	{
		std::string current_time =  boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());
		current_time.erase(0,9);
		return current_time;
	}


	std::string GetCurrentDate()
	{
			std::string current_date = boost::gregorian::to_iso_string(boost::posix_time::second_clock::local_time().date());
			current_date.insert(4,"-");
			current_date.insert(7,"-");

			return current_date;
	}

	#include "Crypto-SHA256.h"
	#include<boost/tokenizer.hpp>

	void OpenSMOKE_logo(const std::string application_name, const std::string author_name)
	{
		std::string current_time = __TIME__;
		std::string current_date = __DATE__;
		std::string author_complete = "Author: " + author_name;
		std::string compilation_time = "Compilation date: " + current_date + " at " + current_time;
		std::string version = "Version: "; version += __OPENSMOKE_VERSION__;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
		std::cout << "          ___                   ____  __  __  ___  _  _______                 " << std::endl;
		std::cout << "         / _ \\ _ __   ___ _ __ / ___||  \\/  |/ _ \\| |/ / ____| _     _        " << std::endl;
		std::cout << "        | | | | '_ \\ / _ \\ '_ \\\\___ \\| |\\/| | | | | ' /|  _| _| |_ _| |_      " << std::endl;
		std::cout << "        | |_| | |_) |  __/ | | |___) | |  | | |_| | . \\| |__|_   _|_   _|     " << std::endl;
		std::cout << "         \\___/| .__/ \\___|_| |_|____/|_|  |_|\\___/|_|\\_\\_____||_|   |_|       " << std::endl;
		std::cout << "              |_|                                                             " << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "           Department of Chemistry, Materials and Chemical Engineering        " << std::endl;
		std::cout << "                              Politecnico di Milano                           " << std::endl;
		std::cout << "                         http://www.opensmoke.polimi.it/                      " << std::endl;
		std::cout << "                      http://creckmodeling.chem.polimi.it/                    " << std::endl;
		std::cout << std::endl;
		for (unsigned int i=1;i<=39-(unsigned int)(application_name.size())/2;i++)	std::cout << " ";
		std::cout << application_name << std::endl;
		for (unsigned int i=1;i<=39-(unsigned int)(version.size())/2;i++)	std::cout << " ";
		std::cout << version << std::endl;
		for (unsigned int i=1;i<=39-(unsigned int)(author_complete.size())/2;i++)	std::cout << " ";
		std::cout << author_complete << std::endl;
		for (unsigned int i=1;i<=39-(unsigned int)(compilation_time.size())/2;i++)	std::cout << " ";
		std::cout << compilation_time << std::endl;
		std::cout << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                                  WARNING                                    " << std::endl;
		std::cout << "   This version of OpenSMOKE++ Suite can be used for educational purposes    " << std::endl;
		std::cout << "              only and cannot be distributed to third parties.               " << std::endl;
		std::cout << "       The software is and remains the sole property of Alberto Cuoci.       " << std::endl;
		std::cout << "      Whenever the OpenSMOKE++ Suite is used to produce any publication,     " << std::endl;
		std::cout << "       a detailed reference to the OpenSMOKE++ code should be reported       " << std::endl;
		std::cout << "                            (see User's Guide).                              " << std::endl;
		std::cout << "    Use for commercial purposes is not permitted. For any commercial issue   " << std::endl;
		std::cout << "         please contact Alberto Cuoci (email: alberto.cuoci@polimi.it)       " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                            LIMITED WARRANTY                                 " << std::endl;
		std::cout << "     This software is provided \"as is\" and without warranties as to        " << std::endl;
		std::cout << "  performance of merchantability or any other warranties whether expressed   " << std::endl;
		std::cout << "    or implied. Because of the various hardware and software environments    " << std::endl;
		std::cout << "   into which this library may be installed, no warranty of fitness for a    " << std::endl;
	        std::cout << "   particular purpose is offered. The user must assume the entire risk of    " << std::endl;
		std::cout << "                          using  the library.                                " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
	}

	std::string SplitStringIntoSeveralLines(std::string source, const std::size_t width, const std::string whitespace)
	{
		std::size_t  currIndex = width - 1;
		std::size_t  sizeToElim;
		while ( currIndex < source.length() )
		{
			currIndex = source.find_last_of(whitespace,currIndex + 1);
			if (currIndex == std::string::npos)
				break;
			currIndex = source.find_last_not_of(whitespace,currIndex);
			if (currIndex == std::string::npos)
				break;
			sizeToElim = source.find_first_not_of(whitespace,currIndex + 1) - currIndex - 1;
			source.replace( currIndex + 1, sizeToElim , "\n");
			currIndex += (width + 1);
		}
		return source;
	}

	void OpenInputFileXML(rapidxml::xml_document<>& doc, std::vector<char>& xml_copy, const boost::filesystem::path& file_name)
	{
		std::ifstream fInput;
		fInput.open(std::string(file_name.string()).c_str(), std::ios::in);
		const std::string string_file = std::string(std::istreambuf_iterator<char>(fInput), std::istreambuf_iterator<char>());
		fInput.close();

		xml_copy.assign(string_file.begin(), string_file.end());
		xml_copy.push_back('\0');

		doc.parse<rapidxml::parse_declaration_node | rapidxml::parse_no_data_nodes>(&xml_copy[0]);
	}

        void OpenInputFileASCII(std::ifstream &fASCII, const boost::filesystem::path input_file_ascii)
        {
            fASCII.open(input_file_ascii.c_str(), std::ios::in);

            if (!fASCII.is_open())
            {
                std::string msg = std::string("The file ") + std::string(input_file_ascii.string()).c_str() + std::string(" could not be open!");
                FatalErrorMessage(msg.c_str());
            }
        }

	void SetXMLFile(std::ostream& xml_string)
	{
		xml_string << std::setprecision(8);
		xml_string.setf(std::ios::scientific);

		xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
	}

        void OpenOutputFileXML(std::ofstream &fXML, const boost::filesystem::path output_file_xml)
        {
            if (!fXML.is_open())
                fXML.open(output_file_xml.c_str(), std::ios::out);

            if (!fXML.is_open())
            {
                std::string msg = std::string("The file ") + std::string(output_file_xml.string()).c_str() + std::string(" could not be open!");
                FatalErrorMessage(msg.c_str());
            }
        }

        void OpenOutputFileASCII(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii)
        {
            if (!fASCII.is_open())
                fASCII.open(output_file_ascii.c_str(), std::ios::out);

            if (!fASCII.is_open())
            {
                std::string msg = std::string("The file ") + std::string(output_file_ascii.string()).c_str() + std::string(" could not be open!");
                FatalErrorMessage(msg.c_str());
             }

            fASCII.setf(std::ios::scientific);
        }

        void OpenOutputFileASCII_Append(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii)
        {
            if (!fASCII.is_open())
                fASCII.open(output_file_ascii.c_str(), std::ios::out | std::ios::app);

            if (!fASCII.is_open())
            {
                std::string msg = std::string("The file ") + std::string(output_file_ascii.string()).c_str() + std::string(" could not be open!");
                FatalErrorMessage(msg.c_str());
             }

            fASCII.setf(std::ios::scientific);
        }

	void PrintTagOnASCIILabel(const unsigned int width, std::ostream &fOut, const std::string tag, unsigned int &counter)
	{
		std::stringstream number;
		number << counter++;
		std::string label = tag + "(" + number.str() + ")";
		fOut << std::setw(width) << std::left << label;
	}

	unsigned int CalculateSpeciesFieldWidth(const std::string species_name, const int max_number_species)
	{
		const double coefficient = 2.3*max_number_species;
		const std::size_t max_number_of_digits = boost::lexical_cast<std::size_t>(std::ceil(std::log10(coefficient)));
		const std::size_t safety_digits = 2;
		const std::size_t minimum_number = 16;

		// name + _ + W + ( max_digits ) + safety
		std::size_t number_to_compare = species_name.size() + 3 + max_number_of_digits + 1 + safety_digits;
		return boost::lexical_cast<unsigned int>(std::max(minimum_number, number_to_compare));
	}

	void CheckKineticsFolder(const boost::filesystem::path& path_output)
	{
		if ( boost::filesystem::is_directory(path_output) == true)
		{
			if ( boost::filesystem::exists(path_output) == false)
			{
				std::string message = "The folder " + path_output.string() + " you specified through the @KineticsFolder option does not exist!";
				FatalErrorMessage(message);
			}

			{
				const boost::filesystem::path kinetics_file = path_output / "kinetics.xml";
				if ( boost::filesystem::exists(kinetics_file) == false)
				{
					std::string message = "The " + kinetics_file.string() + " file does not exist!\n";
					message += "Please check the content of the folder you specified through the @KineticsFolder option";
					FatalErrorMessage(message);
				}
			}

			{
				const boost::filesystem::path reaction_names_file = path_output / "reaction_names.xml";
				if ( boost::filesystem::exists(reaction_names_file) == false)
				{
					std::string message = "The " + reaction_names_file.string() + " file does not exist!\n";
					message += "Please check the content of the folder you specified through the @KineticsFolder option";
					FatalErrorMessage(message);
				}
			}
		}
		else
		{
			std::string message = "The folder " + path_output.string() + " you specified through the @KineticsFolder option does not exist!";
			FatalErrorMessage(message);
		}
	}

	boost::filesystem::path GetExecutableFileName(char** argv)
	{
		#ifdef __linux
		{
			char buf[1024];
			memset(buf, 0, sizeof(buf));
			if (!readlink("/proc/self/exe", buf, sizeof(buf)-1))
			{
				perror("readlink");
				exit(OPENSMOKE_FATAL_ERROR_EXIT);
			}
			boost::filesystem::path executable_file = buf;
			return executable_file;
		}
		#elif __APPLE__
		{
			char path[1024];
			uint32_t size = sizeof(path);
			if (_NSGetExecutablePath(path, &size) == 0)
			{
				boost::filesystem::path executable_file = path;
				return executable_file;
			}
			else
			{
				perror("readlink");
				exit(OPENSMOKE_FATAL_ERROR_EXIT);
			}
		}
		#else
		{
			return boost::filesystem::system_complete(argv[0]);
		}
		#endif
	}

	template<typename T>
	double ErrorControl(const T& w, const T& e)
	{
		const std::size_t ne = w.size();

		double sum = 0.;
		for (unsigned int i = 0; i < ne; i++)
		{
			const double coeff = w(i)*e(i);
			sum += coeff*coeff;
		}

		return std::sqrt( sum/double(ne) );
	}
}

#endif	// OpenSMOKE_OpenSMOKEUtilities_Cpp
