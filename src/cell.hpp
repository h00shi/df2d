/** @file cell.hpp
 * Header file for cell classes.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 * @ingroup mesh_module
 */

#ifndef CELL_HPP
#define CELL_HPP

#include "region.hpp"
#include <armadillo>

/**
 * @addtogroup mesh_module
 * @{
 */

/** @enum CellType
	@brief indicates the geometrical type of a cell.
	Also, the (int)enum returns the number of points of the cell.
*/
enum CellType {CellPoint = 1,      /**< @brief a point cell */
			   CellLine = 2,       /**< @brief a line cell */
			   CellTri = 3,        /**< @brief a triangle cell */
			   CellQuad = 4        /**< @brief a quadrilateral cell */
};
/** @brief stores other properties of a Cell, to be used in template classes.
	the generic type is for n>=3.
 */
template <CellType C>
struct Cell{
	typedef RegionPorousMat Reg;                  /**< @brief region type */
	const static int nPoint = (int)C;             /**< @brief number of corner points */
	const static int nFace =  (int)C;             /**< @brief number of faces */
	const static double rawVol;                   /**< @brief volume in z-e space */
	const static arma::mat::fixed<2,nPoint> del;  /**< @brief del vectors */
	const static arma::mat::fixed<2,nFace> fIp;   /**< @brief flux integ points */
	const static arma::mat::fixed<2,nPoint> vIp;  /**< @brief vol integ points */
	const static arma::vec::fixed<2> center;      /**< @brief center point */
	const static arma::ivec::fixed<nPoint> idxPlus1; /**< @brief indices shifted right */
	const static arma::ivec::fixed<nPoint> idxMin1;  /**< @brief indices shifted left */
	/** @brief Derivative of shape functions matrix	*/
	static void B(const double z,const double e,arma::mat &bm);
	/** @brief shape function matrix */
	static void N(const double z,const double e,arma::rowvec &nm);
private:
	static double inirawvol();                  /**< @brief the value of rawvol */
	static arma::mat inidel();                  /**< @brief the value of del */
	static arma::mat inifip();                  /**< @brief the value of fip */
	static arma::mat inivip();                  /**< @brief the value of vip */
	static arma::vec inicenter();               /**< @brief the value of center */
	static arma::ivec iniidxplus1();          /**< @brief the value of idxplus1 */
	static arma::ivec iniidxmin1();           /**< @brief the value of idxmin1  */
};

/** @brief CellPoint properties.
 */
template <>
struct Cell<CellPoint>{};
/** @brief CellLine properties.
 */
template <>
struct Cell<CellLine>{
	typedef RegionPorousFrac Reg;            /**< region type */
	const static int nPoint = (int)CellLine ;/**< number of corner points */
	const static int nFace  = 1 ;            /**< number of faces */
};
/** @typedef Cell<CellTri> celtri
	@brief triangle cell short name. 
*/
typedef Cell<CellTri> celtri;
/** @typedef Cell<CellQuad> celquad
	@brief quadrilateral cell short name.
*/
typedef Cell<CellQuad> celquad;

/** @typedef Cell<CellLine> cellin
    @brief line cell short name.
*/
typedef Cell<CellLine> cellin;
/** @typedef Cell<CellPoint> celpoi
    @brief point cell short name.
*/
typedef Cell<CellPoint> celpoi;

/** @cond */

/*****************************************************************************
 * common static variables
 ****************************************************************************/
template<CellType C> const double  Cell<C>::rawVol = Cell<C>::inirawvol();
template<CellType C> const arma::mat::fixed<2,Cell<C>::nPoint> Cell<C>::del( Cell<C>::inidel().memptr() ) ; 
template<CellType C> const arma::mat::fixed<2,Cell<C>::nFace> Cell<C>::fIp ( Cell<C>::inifip().memptr() ) ;  
template<CellType C> const arma::mat::fixed<2,Cell<C>::nPoint> Cell<C>::vIp( Cell<C>::inivip().memptr() ) ;  
template<CellType C> const arma::vec::fixed<2> Cell<C>::center  ( Cell<C>::inicenter().memptr() );
template<CellType C> const arma::ivec::fixed<Cell<C>::nPoint> Cell<C>::idxPlus1(Cell<C>::iniidxplus1().memptr());  
template<CellType C> const arma::ivec::fixed<Cell<C>::nPoint> Cell<C>::idxMin1(Cell<C>::iniidxmin1().memptr()) ;

/*****************************************************************************
 * special functions
 * triangle - n=3
 ****************************************************************************/
template<> inline arma::mat celtri::inidel(){
	double pans[] = {1/3.,-1/6. ,-1/6.,1/3. ,-1/6.,-1/6.};
	arma::mat ans(pans, 2, 3);
	return ans;
}
template<> inline arma::mat celtri::inifip(){
	double pans[] = {1/6.,5/12. ,5/12.,1/6. ,5/12.,5/12.};
	arma::mat ans(pans, 2, 3);
	return ans;
}
template<> inline arma::mat celtri::inivip(){
	arma::mat ans("0,0;0,0;0,0");
	return ans;
}
template<> inline arma::vec celtri::inicenter(){
	arma::vec ans("0 0");
	return ans;
}
template<> inline double celtri::inirawvol(){
    return .5;
}
template<> inline arma::ivec celtri::iniidxplus1(){
	arma::ivec ans("1 2 0");
	return ans;
}
template<> inline arma::ivec celtri::iniidxmin1(){
	arma::ivec ans("2 0 1");
	return ans;
}
template<> inline void celtri::B(const double z,const double e,arma::mat &bm){
	bm << -1 << -1 << arma::endr
	   << 1 << 0 << arma::endr
	   << 0 << 1 << arma::endr;
}
template<> inline void celtri::N(const double z,const double e,arma::rowvec &nm){
	nm << 1-z-e << z << e << arma::endr;
}

/*****************************************************************************
 * special functions
 * quad - n=4
 ****************************************************************************/
template<> inline arma::mat celquad::inidel(){
	arma::mat ans(".5,0;0,.5;-.5,0;0,-.5");
	return ans.t();
}
template<> inline arma::mat celquad::inifip(){
	arma::mat ans(".25,.5;.5,.25;.75,.5;.5,.75");
	return ans.t();
}
template<> inline arma::mat celquad::inivip(){
	arma::mat ans(".25,.25;.75,.25;.75,.75;.25,.75");
	return ans.t();
}
template<> inline arma::vec celquad::inicenter(){
	arma::vec ans("0.5 0.5");
	return ans;
}
template<> inline double celquad::inirawvol(){
    return 1;
}
template<> inline arma::ivec celquad::iniidxplus1(){
	arma::ivec ans("1 2 3 0");
	return ans;
}
template<> inline arma::ivec celquad::iniidxmin1(){
	arma::ivec ans("3 0 1 2");
	return ans;
}
template<> inline void celquad::B(const double z,const double e,arma::mat &bm){
	bm << -(1-e) << -(1-z) << arma::endr
	   << 1-e << -z << arma::endr
	   << e << z << arma::endr
	   << -e << 1-z << arma::endr;
}
template<> inline void celquad::N(const double z,const double e,arma::rowvec &nm){
	nm << (1-z)*(1-e) << z*(1-e) << z*e << (1-z)*e <<  arma::endr;
}
/**< @endcond */

/** @}
 */

#endif /*CELL_HPP*/
