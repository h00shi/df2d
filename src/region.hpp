/** @file region.hpp
 * Header file for Region data structure.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 @ingroup edat_module
 */

#ifndef REGION_HPP
#define REGION_HPP

#include "armadillo"
#include <string>

/** @addtogroup edat_module
 *	@{
*/

// The header file for the KFunc class is included in the .cpp file to prevent
// mutual including of header files.
class KFunc;

/***************************************************************************
 * Abstract class
 **************************************************************************/

/** @brief Abstract structure for storing a region.
 *
 * This is a region
 */
struct Region{
	/** @brief The number that the mesh generator has given to the member elements
		of this region.
	*/
	int ID;   
	/** @brief Priority to be master region.
		
		As I explained in my thesis, one of the porous media regions
		(not the boundary region of course) that are in contact, should be
		chosen as the master region. Which region is the master, depends on
		the J curve used. For example the Brooks-Corey J curves require the
		regions with the highest Pd to be the master region. The higher this
		idx number, the higher the priority of the region to be the master
		region.

		Because the priority depends on the J curve, this number is assigned
		by the JFunc class, through its member named RegionPointerComparer.
	*/
	int	idx;  
	/** @brief returns true if the region is a boundary codition region and
		false if the region is a porous media region.
	 */
	virtual bool isBoundary() const = 0;
	/** @brief returns a string that describes the region */
	virtual std::string name() const = 0;
	/** @brief constructor that sets idx to zero.*/
	Region();
	/** @brief virtual destructor that is present because the class is abstract */
	virtual ~Region() {}    
};

/***************************************************************************
 * Porous media regions
 **************************************************************************/

/** @brief Abstract class for a porous media region.
	
	May be Matrix or Fracture.
 */
struct RegionPorous : public Region{
	double phi;	/**< @brief Porosity. */
	double pd;  /**< @brief Cappilary Coefficient. */
	/** @brief Relative permeability.

		Because the relative permeability classes (KFunc) are polymorphic
		classes, this pointer can point to any of them.
	 */
	KFunc *kr;  
	bool isBoundary() const;
	/** @brief returns the dimensions of the region.

	    Matrix regions are 2d while fracture regions are 1d.
	 */
	virtual int dim() const = 0;
protected:
	/** @brief Sets kr to NULL . */
	RegionPorous();
	/** @brief Destroys kr. */
	virtual ~RegionPorous();
};

/** @brief region for matrix.

	@note
	To create one you have to
    \li set k using constructor\n
	\li set ID\n
	\li set phi\n
	\li set kr\n
    \li After all regions are created they should be indexed.
	
	To destroy one using the destructor is sufficient (takes care of kr).
 */
struct RegionPorousMat : public RegionPorous{
	arma::mat::fixed<2,2> k; /**< @brief permeability tensor */
	int dim() const ;
	std::string name() const;
	/**@brief Initializes k */
	RegionPorousMat(const std::string&);
};

/** @brief region for fracture.

	@note
	To create one you have to
	\li set k 
	\li set e 
	\li set ID
	\li set phi
	\li set kr
	After all was created the list should be sorted and indexed
	
	To destroy one using delete is sufficient (takes care of kr).
 */
struct RegionPorousFrac : public RegionPorous{	
	/**< @brief permeability scalar

	   While K is a matrix for the matrice is a scalar for fractures because they are
	   one-dimensional.
	 */
	double k; 
	double e; /**< @brief fracture thickness */
	int dim() const ;
	std::string name() const;				
};

/***************************************************************************
 * Boundary region
 **************************************************************************/

/** @brief struct for boundary region.

	@note
	To create one you have to
	\li set ID
	\li set idx
	\li set ptype
	\li set stype
	\li set value
	\li After all the regions are created they should be numbered.
	
	To destroy one, using delete is sufficient. 
 */
struct RegionBoundary : public Region{
	/** @brief length of val[] */
	const static uint VAL_MEM = 1;	
	/** @brief Saturation boundary condition enum. */
	enum SType{ SSConst,            /**< @brief constant saturation */
				SGPZero             /**< @brief zero cappilary surface gradient */
	};	
	/** @brief Pressure boundary condition enum */
	enum PType{ PPConst,        /**< @brief constant pressure */
				PQConst         /**< @brief constant total in(out)flow */
	};	
	SType stype; 	/**< @brief Saturation Boundary Condition */		
	PType ptype;    /**< @brief Pressure boundary condition   */
	/** @brief value of boundary patch.
	 *	
	 *  stores the inflow for cq and pressure for cp boundary conditions.
	 */
	double val[VAL_MEM];
	bool isBoundary() const;
	std::string name() const;
};

/**
 * @}
 */


#endif /*REGION_HPP*/
 
