/* @file region.cpp
 * .cpp file for Region classes.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 */

#include "jkfunc.hpp"
#include "region.hpp" 
#include <sstream>

/***************************************************************************
 * Region
 **************************************************************************/
Region::Region():ID(0),idx(0){}
/*****************************************************************************
 * RegionPorous
 ****************************************************************************/
bool RegionPorous::isBoundary() const{
	return false;
}

RegionPorous::~RegionPorous(){
	if(kr) delete kr;
	kr = NULL;
}

RegionPorous::RegionPorous(){
	kr = NULL;
}

/*****************************************************************************
 * RegionPorousMat
 ****************************************************************************/

int RegionPorousMat::dim() const {
	return 2 ;
}

std::string RegionPorousMat::name() const{
	std::stringstream ss;
	ss << "MatrixRegion: |ID: " << ID << " |idx: " << idx
		 << " |phi: " << phi << " |pd: " << pd
		 << " |kr: " << kr->name()
		 << " |k: "
		 << k(0,0) << ","<< k(1,0)<< "," << k(0,1) << ","<< k(1,1) << " |";
	return ss.str();
}

RegionPorousMat::RegionPorousMat(const std::string &k_val):k(k_val){
}

/*****************************************************************************
 * RegionPorousFrac
 ****************************************************************************/

int RegionPorousFrac::dim() const {
	return 1 ;
}

std::string RegionPorousFrac::name() const{
	std::stringstream ss;
	ss << "FractureRegion: |ID: " << ID << " |idx: " << idx
		 << " |phi: " << phi << " |pd: " << pd
		 << " |kr: " << kr->name()
		 << " |k: " << k << " |e: " << e << "|";
	return ss.str();
}

/*****************************************************************************
 * RegionBoundary
 ****************************************************************************/

std::string RegionBoundary::name() const{
	std::stringstream ss;
	ss << "BoundaryRegion: |ID: " << ID << " |idx: " << idx
		 << " |SType: " << (stype==SSConst ? "constant_s" : "(grad_pc).n=0")
		 << " |PType: " << (ptype==PPConst ? "constant_p" : "constant_q")
		 << " |value: " << val[0] << "|";
	return ss.str();
}

bool RegionBoundary::isBoundary() const{
	return true;
}
