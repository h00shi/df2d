/** @file error.cpp
 * cpp file for error.hpp.
 * Written by Shayan Hoshyari.
 * As a part of discrete fracture project.
 */

#include "error.hpp"

int Error::code = 1 ;
int Error::line = -2;
std::string Error::func = "";
std::string Error::file = "";
bool Error::ini = false;
std::ostringstream Error::mess;

void Error::sendtopetsc()
{
	PetscError(PETSC_COMM_SELF,
						 line,
						 func.c_str(),
						 file.c_str(),
						 code,
						 (ini ? PETSC_ERROR_INITIAL : PETSC_ERROR_REPEAT),
						 mess.str().c_str() );
}

void Error::print(){
	std::cout << "******************************\n"
						<< "* code: " << code << std::endl
						<< "* line: " << line << std::endl
						<< "* func: " << func << std::endl
						<< "* file: " << file << std::endl
						<< "* initial: " <<  ini  << std::endl
						<< "* mess: " << mess.str()  << std::endl
						<< "******************************\n" ;
}

