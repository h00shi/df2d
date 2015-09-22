/** @file mdata.hpp
 * Header file for main data structure.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 @ingroup edat_module
 */

#ifndef MDATA_HPP
#define MDATA_HPP

#include "jkfunc.hpp"
#include "region.hpp"
#include <list>
#include <vector>
#include <petscksp.h>
#include "gravity.hpp"


/** @ingroup edat_module
 * @brief This structure keeps the main data needed by our program.
 *
 * After the program starts, an instance of this class is created. All the data
 * (except mesh data) related to the program are stored in this class. All the
 * @ref dr_module functions pass this instance to each other to work with it.
 *
 * Alternatively these parameters could have been declared global but I
 * prefered this strategy.
 */
struct MData{
	
	/*******************************************************************
	 * Nodal Values
	 ******************************************************************/	
	/** @brief Saturation at each node - Discontinuous */
	std::vector<double> S;
	/** @brief Pointer to read Pvec - Continuous  */
	double const *P;
	/** @brief PetscVector used to store WettingPhasePotential - Continuous  */
	Vec Pvec;
	/** @brief Capillary Potential at each node - Discontinuous. */
	std::vector<double> Pc;
	/** @brief stores (gn-gw)*h for each node */
	std::vector<double> DgH;

	
	/** @brief Volume of each subcontrol volume * porosity  - Discontinuous */
	std::vector<double> VPhi;
	/** @brief Sum of V * phi * ds_master/ds (almost LHS of S equation) - Continuous */
	std::vector<double> SPhiV;
	/** @brief Wetting Phase Mobility - Discontinuous */
	std::vector<double> Lw;
	/** @brief Nonwetting Phase Mobility - Discontinuous */
	std::vector<double> Ln;

	
	/** @brief Sum of S fluxes (RHS of S equation) - Continuous */
	std::vector<double> Fs;
	/** @brief Change in saturation - Discontinuous  */
	std::vector<double> dS;
	/** @brief LHS of P equation - Petsc Mat	 */
	Mat A;
	/** @brief RHS of P equation- Petsc Vector  */
	Vec b;
	/** @brief KSP for p equation - Petsc Krylov SubsPace solver*/
	KSP ksp;
	

	/*******************************************************************
	 * Single Values
	 ******************************************************************/
	Gravity grav; /**< @brief data related to gravity field and densities */
	
	double dm,  /**< @brief dimensionless number mu_n/mu_w */
		dn,       /**< @brief dimensionless number L*mu_w/P/K * L/t */
		dp;       /**< @brief dimensionless number L*mu_w/P/K * Q */

	double dtm, /**< @brief minimum time step */
		dtM,      /**< @brief maximum time step */
		dsm,      /**< @brief minimum change in saturation */
		dsM,      /**< @brief maximum change in saturation */ 
		beta;     /**< @brief coefficient to change time step */

	double t,   /**< @brief current time */
		dt,       /**< @brief current time step */
		t0,       /**< @brief initial time */
		tEnd ,    /**< @brief end time */
		tWrite;   /**< @brief period of time used for writing files */

	int nIt,   /**< @brief number of times that s equation is solved from beginning */
		dnIt,    /**< @brief number of times that s equation is solved in this timestep*/
		dnItM;  /**< @brief maximum allowed number of times to solve S in one timestep*/

	double dcT, /**< @brief time passed for one time step */
		cT;       /**< @brief time passed since the program has started */

	double qIn, /**< @brief total injected fluid to reservoir*/
		qOut,     /**< @brief total extracted fluid from reservoir */
		qWin,     /**< @brief total water injected to reservoir */
		qWout;    /**< @brief total water extracted from reservoir */

	int nFile0,   /**< @brief first number of .vtk (octave for 1d) file */
		nFile,      /**< @brief current number of .vtk (octave for 1d) file */
		NFile;      /**< @brief final number of .vtk (octave for 1d) file */

	JFunc *J;        /**< @brief Cappilary Curve */

	std::string dir;  /**< @brief current directory */

	/** @brief indicates the type of mesh input file */
	enum MeshType{MeshTriangle, /**< @brief Triangle .ele + .node + .edge file */
				  MeshGmsh      /**< @brief gmsh .msh file */
	};
	/** @brief type of mesh we have to read  */
	MeshType meshtype;
	
	/** @brief indicates type of output file */
	enum VisualType {VisualVtk,     /**< vtk file */
					 VisualTecplot,  /**< Tecplot file */
					 VisualGmsh     /**< msh file */
	};
	/**@brief  type of output */
	VisualType visualtype;
	/**@brief  How to write discontinuous data into vtk files.
	 *
	 * Multiple values of saturation may exist at each node, So there are multiple options
	 * available for writing them into a visualization, e.g. vtk, file. If this variable is
	 * set to 2 the max(S) will be written as the nodes saturation. If set to 1 the min(S)
	 * will be written. Finally if this variable is 0, df2d will create as many nodes as
	 * DuplDatas in the output files and writes each saturation for one of them.
	 */
	int visualduplicate;

	/** @brief setfield mode or solver mode.
	 *
	 * If set to 1 df2d only changes the initial file and exits.
	 * If set to 0 df2d will run a simulation
	*/
	int setfield;
	

	/*******************************************************************
	 * Functions 
	 ******************************************************************/
	/** @brief Sets some pointers to NULL and some vars to 0.
	 */
	void initialize();
	/** @brief Deallocates memmory for the dynamic or petsc types.
	 */
	void finalize();
};

#endif /*MDATA*/
