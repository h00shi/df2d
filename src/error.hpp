/** @file error.hpp
 * contains macros, and a class for error propagation.
 * Written by Shayan Hoshyari.
 * As a part of discrete fracture project.
 * @ingroup dr_module
 */

#if !defined(ERROR_H)
#define ERROR_H

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <petscsys.h>

/** @brief Data and functions for an error.
 *  @ingroup dr_module
 */
class Error
{
public:
	/** @brief Petsc error code . */
  static int code ;
	/** @brief The line at which the error has occured.
	 *
	 * due to the nature of try catch only works for the first time.
	 */
	static int line ;
	/** @brief The func in which the error has occured. */
	static std::string func;
	/** @brief The file at which error has occured. */
	static std::string file;
	/** @brief Show if we are at the beginning of the error propagation.  */
	static bool ini;
	/** @brief Description of why the error has occured.  */
	static std::ostringstream mess;
	/** @brief Send the error to Petsc so it can be printed.  */
	static void sendtopetsc();
	/** @brief Print the data on screen for debugging. */
	static void print();	    
};

/** @def ERRCHK
 * @brief Replaces Petsc's CHKERRQ(ierr).
 *
 * usage: Error::code=PetscFunc(...);ERRCHK();
 @ingroup dr_module
 */
#define ERRCHK() do{													\
		if(Error::code != 0){											\
			Error::line = __LINE__;										\
			Error::ini = false;											\
			throw *new Error;											\
		}																\
	}	while(0)

/** @def ERRSET
 * @brief  Replaces Petsc's SETERRQ(...).
 
 * usage:
 * Error::mess << (your error message);
 * ERRSET();
 @ingroup dr_module
 */
#define ERRSET() do{													\
		Error::line = __LINE__;											\
		Error::ini = true;												\
		throw *new Error;												\
	}while(0)

/** @def FuncBegin
 * @brief  Put at the beginning of every function.
 
 @ingroup dr_module
 */
#define FuncBegin() try{

/** @def FuncEnd
 * @brief  Put at the end of every function.
 
 @ingroup dr_module
 */
#define FuncEnd() 	}													\
		catch(Error &e){												\
			Error::file = __FILE__;										\
			Error::func = __func__;										\
			Error::sendtopetsc();										\
			Error::line = -1;											\
			Error::ini = false;											\
			throw e;													\
		}																\
		catch(std::exception &e){										\
			Error::file = __FILE__;										\
			Error::func = __func__;										\
			Error::ini = true;											\
			Error::mess.str("");										\
			Error::mess << e.what();									\
			Error::sendtopetsc();										\
			Error::ini = false;											\
			throw *new Error;											\
		}																\
		catch(...){														\
			Error::file = __FILE__;										\
			Error::func = __func__;										\
			Error::ini = true;											\
			Error::mess.str("");										\
			Error::mess << "UNKNOWN EXCEPTION OCCURED";					\
			Error::sendtopetsc();										\
			Error::ini = false;											\
			throw *new Error;											\
		}																\
		
#endif /*ERROR_H*/
