/** @file driver.hpp
	input-output-timemarcher functions.

	Written by Shayan Hoshyari as part of DF_2d.
	@ingroup dr_module
*/

#ifndef DRIVER_HPP
#define DRIVER_HPP

#include "mdata.hpp"
#include "mesh.hpp"


/** @namespace driver
	@brief contains functions for input, output and timemarcher.
	@ingroup dr_module
*/
namespace driver{

	/** @brief sets the current dir and inits petsc
		@ingroup dr_module
	 */
	void initialize(int *argc, char **argv[], MData &md);
	/** @brief read the region data 
		@ingroup dr_module
	*/
	void readregion (MData &md, Mesh &msh);
	/** @brief read the fixed data 
		@ingroup dr_module
	*/
	void readfixed (MData &md);
	/** @brief read mesh general 
		@ingroup dr_module
	*/
	void readmesh(MData &md, Mesh &msh);
	/** @brief read initial condition 
		@ingroup dr_module
	*/
	void readinitial(MData &md, Mesh &msh);
	/** @brief change the data in the field 
		@ingroup dr_module
	*/
	void setfield(MData &md, Mesh &msh);
	/** @brief assign memory to vectors 
		@ingroup dr_module
	*/
	void preparedata (MData &md, Mesh &msh);
	/** @brief march in time 
		@ingroup dr_module
	*/
	void marchintime(MData &md, Mesh &msh);
	/** @brief write the results if the time has come 
		@ingroup dr_module
		@param force if force is true the data will be written anyways
	*/
	void writeintime(MData &md, Mesh &msh, bool force);
	

	/** @brief write initial condition and restart file 
		@ingroup dr_module
	*/
	void writeinitial(MData &md, Mesh &msh,const std::string& adr,
					  const bool restart, const bool octave);
	/** @brief write visualization file general 
		@ingroup dr_module
	*/
	void writevisual(const char * address,MData &md,Mesh &msh);	
	/** @brief create a file with details for debugging 
		@ingroup dr_module
	*/
	void writetest (MData &md, Mesh &msh, const uint num);

}


#endif /*DRIVER_HPP*/
