/** @file element_bone.hpp
 * Header file for element classes - with only declarations.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 * @ingroup mesh_module
 */

#ifndef ELEMENT_BONE_HPP
#define ELEMENT_BONE_HPP

#include "node.hpp"
#include "formula.hpp"
#include "cell.hpp"

/*****************************************************************************
 * Element - first level
 ****************************************************************************/
/** @ingroup mesh_module
	@brief Abstract element class.
	All the element pointers in the code must be of this type.

	@note
	To create one: 
	@li constructor
	@li constructDuplData    - you may now tag node dd
	@li constructBVertices   - you may now constructGeoParams for bvertices
	@li constructGeoParams   
	@li fndUp                - you may now call rhs and lhs funcs
	
	To use:
	@li fndUp
	@li lhs,rhs
	@li matKD - matVol
	
	To destory:
	@li deleting is sufficient
 */
class Element{
public:
	/** @brief Creates a new element.

	 @param reg a pointer to the region to which the element should belong
	 @param nd  an array containing pointers to member nodes
	 @param celltype the celltype of the element (point, line, triangle, quad)
	 @returns allocates memory for a new element and returns the pointer to it.
	*/
	static Element* neww(RegionPorous *reg, Node *nd[], const CellType celltype);
	/** @brief returns element's cell type */
	virtual CellType cellType() const = 0;
	/** @brief returns number of corner points	 */
	virtual int nNode() const = 0;
	/** @brief returns number of faces (flux integration points) */
	virtual int nFace() const = 0;
	/** @brief evaluate upwind nodes for wetting phase.
		
		@param P the pointer to pressure data array.
	 */
	virtual void fndUpW(double const *P) = 0;
	/** @brief evaluate upwind nodes for non-wetting phase.
		
		@param P the pointer to pressure data array
		@param Pc the vector of cappilary pressure
	 */
	virtual void fndUpN(double const *P, const std::vector<double> &Pc) = 0;
	/** @brief finds the volume of sub-control volumes in the element.
		
		@returns a rowvec containing the volume of each sub-CV
	*/
	virtual const arma::rowvec& matVolume() = 0; 
	/** @brief find the D matrix.
		
		@returns the D matrix evaluated at the center point. D * P = Grad P.
		@note
		The retuned value is a reference value. It means that there is a static (much like
		global) arma::rowvec member of ElementBase class and calling this function will
		change that static member and return a reference value to it. So
		if your have previously called matJ, matB or matX and stored pointers to them, i.e.:
		@verbatim
		arma::mat *J;
		J = &element.matJ(xi,eta);
		@endverbatim
		your pointers will point to new values after calling this function!
	*/
	virtual const arma::mat& matKD() = 0; 
	/** @brief global indices of local nodes.
		
		@returns the indices of element corner nodes.
		@note
		The retuned value is a reference value. It means that there is a static (much like
		global) arma::rowvec member of ElementBase class and calling this function will
		change that static member and return a reference value to it. 
	 */
	virtual const arma::ivec& idxGlob() = 0;
	/** @brief local data of corner nodes for a dis-continuous (stored at dupldatas) global data.
		
		@param dat the global data to extract local data from
		@param i the index of the static vector to which the data is written
		@returns a vector containing local data!
	*/
	virtual const arma::vec& lDatCnDis (const std::vector<double> &dat, const uint i)= 0; 
	/** @brief local data of corner nodes for a continuous global data (stored at nodes) - std::vector.
		
		@param dat the global data to extract local data from
		@param i the index of the static vector to which the data is written on.
		@returns a vector containing local data.
		
		@note
		The retuned value is a reference value. It means that there is a static (much like
		global) arma::rowvec member of ElementBase class and calling this function will
		change that static member and return a reference value to it. So this is a right
		usage:
		@verbatim
		arma::vec *vec1,*vec2;
		vec1 = &elem.lDatCnDis(lw,1);
		vec2 = &elem.lDatCnDis(ln,2);
		@endverbatim
		While this is a wrong one and makes the value of vec1 invalid:
		@verbatim
		vec1 = &elem.lDatCnDis(lw,1);
		vec2 = &elem.lDatCnDis(ln,1);
		@endverbatim
		Note that if you do not use pointers everything will be ok since the data will
		be copied. However this will take additional run-time
	 */
	virtual const arma::vec& lDatCnCon (const std::vector<double> &dat, const uint i)= 0; 
	/** @brief local data for corner nodes for a continuous global data (stored at nodes) - double*.
		
		@param dat the global data to extract local data from
		@param i the index of the static vector to which the data is written
		@returns a vector containing local data
		
		@note
	    refer to the note in lDatCnDis.
	*/
	virtual const arma::vec& lDatCnCon (const double* dat, const uint i)= 0; 	
	/** @brief local lhs matrix for p equation.
		
		@param Lw vector containing the wetting phase mobilites
		@param Ln vector containing the non-wetting phase mobilities.

		@note
	    refer to the note in lDatCnDis.
	 */
	virtual const arma::mat& lhsP(const std::vector<double> &Lw, const std::vector<double> &Ln)= 0; 
	/** @brief local lhs matrix for p equation.
		
		@param Ln vector containing the non-wetting phase mobilities.
		@param Pc vector containing the cappilary pressure.
		@param i the index of the static vector to which the data is written
	 */
	virtual const arma::vec& rhsP(const std::vector<double> &Ln, const std::vector<double> &Pc, const uint i) =0;
	/** @brief local lhs matrix for p equation.
		
		@param Lw vector containing the wetting phase mobilites
		@param P pointer to pressure data
		@param i the index of the static vector to which the data is written

		@note
	    refer to the note in lDatCnDis.
	 */
	virtual const arma::vec& rhsS(const std::vector<double> &Lw, const double *P, const uint i) = 0;
	/** @brief creates dupldata list for nodes and link it to element */
	virtual void constructDuplData() = 0;
	/** @brief checks the orientation of the bvertices if the element has any.
		@returns the result of checking neighbourhood
	 */
	virtual void constructBVertices(int *res) = 0;
	/** @brief creates other internal data. like H matrices for polygons and ke_l for lines */
	virtual void constructGeoParams() = 0;
	/** @brief destructor for no compiler-warning.
	 */
	virtual ~Element() {}
	/** @brief returns the description of element as string */
	virtual std::string name(const std::vector<double> *S = NULL,
							 const double *P = NULL, const std::vector<double> *Pc = NULL,
							 const std::vector<double> *Lw = NULL, const std::vector<double> *Ln = NULL) = 0;
	/** @brief returns the regions id */
	virtual int regionID() const = 0;
	/** @brief returns the ith dupldata */
	virtual const DuplData* dupl(const int i) const = 0;
};
/** @typedef Element eleblank
	@brienf abstract element aliasing.
 */
typedef Element eleblank;

/*****************************************************************************
 * ElementBase - second level
 ****************************************************************************/

/** @brief Common base class for all element types.

	At first I wanted to create only one template class for all elements (frac,tri,quad).
	While the tri and quad elements could be perfectly programmed in one template class,
	the fracture element was a bit different. So first I created the ElementBase class,
	which is common between fractures and poly (tri, quad) elements and then inheritted
	ElementPoly class from it. The ElementPoly class is specialized for the fracture elements,
	but is common between quads and triangles.
	@ingroup mesh_module
 */
template<CellType C>
class ElementBase: public Element{
protected:
	typename Cell<C>::Reg *reg_;         /**< @brief The region that the element belongs to */
	Node *nd_[ Cell<C>::nPoint ];        /**< @brief pointer to corner nodes */
	DuplData *dd_[ Cell<C>::nPoint ];    /**< @brief pointer to corner dupldata */
	int upwetidx_[Cell<C>::nFace],          /**< @brief upwind index for wetting phase*/
		upnonidx_[ Cell<C>::nFace ];       /**< @brief upwind index for non-wetting phase */
	
	const static uint nSafe_ = 5;            /**< @brief number of internal local data vectors */
	static arma::mat::fixed<Cell<C>::nPoint,Cell<C>::nPoint> matLdCn_;   /**< @brief corner node local matrix */
	static arma::vec::fixed<Cell<C>::nPoint> vecLdCn_[nSafe_];   /**< @brief corner node local data double */ 
	static arma::vec::fixed<Cell<C>::nFace> vecLdFc_[nSafe_];    /**< @brief face local data double */
	static arma::ivec::fixed<Cell<C>::nPoint> ivecLdCn_;         /**< @brief corner node local data int */
	static arma::mat::fixed<2, Cell<C>::nPoint> KD_;              /**< @brief derivative matrix: KD*P = K\\nabla\\cdot P */
	static arma::rowvec::fixed<Cell<C>::nPoint> V_;              /**< @brief volume row vector */
	/** @brief initializes all data to zero and null and sets region and nodes	 */
	ElementBase(RegionPorous* reg, Node *nd[]);
	/** @brief string name for that part of element which is related to ElementBase  */
	std::string nameBase(const std::vector<double> *S = NULL,
						  const double *P = NULL, const std::vector<double> *Pc = NULL,
						  const std::vector<double> *Lw = NULL, const std::vector<double> *Ln = NULL) ;
public:
	CellType cellType() const;
	int nNode() const;
	int nFace() const;
	const arma::ivec& idxGlob();
	const arma::vec& lDatCnDis (const std::vector<double> &dat, const uint i);
	const arma::vec& lDatCnCon (const std::vector<double> &dat, const uint i);
	const arma::vec& lDatCnCon (const double* dat, const uint i);
	/** @brief local upwinded data of faces for a dis-continuous global data - std::vector.
		
		@param dat the global data to extract local data from
		@param idx the index of upwind nodes
		@param i the index of the static vector to which the data is written
		@returns a vec containing the local data
	 */
	const arma::vec& lDatUpDis (const std::vector<double> &dat, const int idx[], const uint i);
	void constructDuplData();
	int regionID() const ;
	const DuplData* dupl(const int i) const;
};

/*****************************************************************************
 * ElementPoly - third level
 ****************************************************************************/

/** @brief Class for n>=3 elements.
	@ingroup mesh_module
 */
template<CellType C>
class ElementPoly: public ElementBase<C>{
private:
	typedef ElementBase<C> f;                                 /**< @brief makes writing the code easier */
protected:
	/** @brief H matrix same as the thesis.
		
		This member is not static and is saved for each element.
	 */
	arma::mat::fixed< Cell<C>::nPoint , Cell<C>::nPoint > H_; /**< @brief H  matrix */
	/** @brief N matrix same as the \\psi matrix in the thesis.

		This static variable saves some memmory for the N matrix. Whenever
		an element's matN() function is called it will fill this static
		variable and return a reference to it. 
	 */
	static arma::rowvec::fixed< Cell<C>::nPoint > N_;        
	/** @brief B matrix same as the \\frac{d\\psi}{d \\xi} matrix in the thesis.

		This static variable saves some memmory for the B matrix. Whenever
		an element's matB() function is called it will fill this static
		variable and return a reference to it. 
	 */
	static arma::mat::fixed <Cell<C>::nPoint, 2> B_;         
	/** @brief X matrix same as the thesis.

		This static variable saves some memmory for the X matrix. Whenever
		an element's matX() function is called it will fill this static
		variable and return a reference to it. 
	 */
	static arma::mat::fixed <2, Cell<C>::nPoint> X_;         
	/** @brief J matrix same as the thesis.

		This static variable saves some memmory for the J matrix. Whenever
		an element's matJ() function is called it will fill this static
		variable and return a reference to it. 
	 */
	static arma::mat::fixed <2, 2> J_;                      
 
public:
	void fndUpW(double const *P);
	void fndUpN(double const *P, const std::vector<double> &Pc);
	const arma::mat& matKD();
	const arma::rowvec& matVolume();
	
	const arma::mat& lhsP(const std::vector<double> &Lw, const std::vector<double> &Ln);
	const arma::vec& rhsP(const std::vector<double> &Ln, const std::vector<double> &Pc, const uint i);
	const arma::vec& rhsS(const std::vector<double> &Lw, const double *P, const uint i);
	/** @brief Shape function rowvec.
	 * @note modifies no one
	 */
	 const arma::rowvec& matN(const double z, const double e);
	/** @brief Shape function derivative.
		@note modifies no one
	 */
	 const arma::mat& matB(const double z, const double e);
	/** @brief coordinate matrix.
		@note modifies no one
	 */
	 const arma::mat& matX();
	/** @brief jacobian matrix.
		@note modifies X and B and replaces them with the value at z,e
	 */
	 const arma::mat& matJ(const double z, const double e);
	/** @brief initializes to null and zero and sets nodes and reg.
	 */
	ElementPoly(RegionPorous* reg, Node *nd[]);
	void constructBVertices(int *res);
	void constructGeoParams();
	std::string name(const std::vector<double> *S = NULL,
					 const double *P = NULL, const std::vector<double> *Pc = NULL,
					 const std::vector<double> *Lw = NULL, const std::vector<double> *Ln = NULL);
};
/** @typedef ElementPoly<CellTri> eletri
	@brief triangle element aliasing.
*/
typedef ElementPoly<CellTri> eletri;
/** @typedef ElementPoly<CellQuad> elequad
	@brief quadrilateral element aliasing.
*/
typedef ElementPoly<CellQuad> elequad;

/** @brief ElementPoly class specialized for n=2 (fracture) elements.
	
	@ingroup mesh_module
 */
template<>
class ElementPoly<CellLine>: public ElementBase<CellLine>{
protected:
    double KE_L_;             /**< @brief k*e/l a paremeter used to create local matrices */
 
public:
	void fndUpW(double const *P);
	void fndUpN(double const *P, const std::vector<double> &Pc);
	const arma::rowvec& matVolume(); 
	const arma::mat& matKD();
	
	const arma::mat& lhsP(const std::vector<double> &Lw, const std::vector<double> &Ln);
	const arma::vec& rhsP(const std::vector<double> &Ln, const std::vector<double> &Pc, const uint i);
	const arma::vec& rhsS(const std::vector<double> &Lw, const double *P, const uint i);
	
	/** @brief initializes to null and zero and sets nodes and reg.
	 */
	ElementPoly(RegionPorous* reg, Node *nd[]);
	
	void constructBVertices(int *res) {}
	void constructGeoParams();

	std::string name(const std::vector<double> *S = NULL,
					 const double *P = NULL, const std::vector<double> *Pc = NULL,
					 const std::vector<double> *Lw = NULL, const std::vector<double> *Ln = NULL) ;
};
/** @typedef ElementPoly<CellLine> elefrac
	@brief fracture  element aliasing.
 */
typedef ElementPoly<CellLine> elefrac;


#endif /*ELEMENT_BONE_HPP*/

