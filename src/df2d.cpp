/** @file df2d.cpp
	.cpp file containing the main function for DF_2d program.
*/
	
#include "driver.hpp"

/** The main function, where everything the program does is called from . */
int main(int argc, char *argv[]){
	FuncBegin();

	MData md; //mesh
	Mesh msh; //external variables
	
	driver::initialize(&argc,&argv,md); //initialize all the data
	driver::readregion(md,msh);         //read the regions
	driver::readfixed(md);              //read other solver properties
	driver::readmesh(md,msh);           //read the mesh
	driver::readinitial(md, msh);       //read the initial conditions

    if(md.setfield){                    //we should change the initial conditions
		driver::setfield(md, msh);     
	}
	else{                               //we should run a simulation
		driver::preparedata(md,msh);
		do{
			driver::marchintime(md, msh);        //solve the system once in time
			driver::writeintime(md, msh,false);  //write the data if required
		}	while ( md.t < md.tEnd);
	}

	//anounce that everything is done
	std::cout << "Simulation finished successfully in " << md.cT << " seconds." << std::endl
			  << "The results can be found in " << md.dir << "result/(*.vtk and *.flow)" << std::endl ;
	
	//finalize petsc and other data
	md.finalize();
	PetscFinalize();
	return 0;
	FuncEnd();
}
