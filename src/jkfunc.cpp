/* @file jkfunc.cpp
 * .cpp file for jkfunc.
 */

#include "jkfunc.hpp"
#include "region.hpp"

/***************************************************************************
 * RegionPointerComparer
 **************************************************************************/

bool RegionPointerComparer::operator() (Region const* r1, Region const* r2) const{
	// if both are boundary we do not care about their priority, I just
	// sorted them by mesh generator id.	
	if (r1->isBoundary() && r2->isBoundary()) return (r1->ID < r2->ID);
	// if one of them is boundary the porous media should have a bigger
	// index.
	else if( r1->isBoundary() ) return false;
	else if( r2->isBoundary() ) return true;
	// if both regions are porous media it is up to the jfunc to decide
	// which one has the priority to be the master region
	else return master_->compareRegion((RegionPorous*)r1,(RegionPorous*)r2);			
}

RegionPointerComparer::RegionPointerComparer(JFunc const* master_val):master_(master_val){
}

/***************************************************************************
 * JFunc
 **************************************************************************/
JFunc::JFunc():cmp(this){
}

bool JFuncFirooz::compareRegion(RegionPorous const *r1, RegionPorous const *r2) const {
	return (r1->pd > r2->pd);
}
std::string JFuncFirooz::name() const{
	return "FiroozCappilaryCurve";
}

bool JFuncZero::compareRegion(RegionPorous const *r1, RegionPorous const *r2) const {
	return (r1->ID < r2->ID);
}
std::string JFuncZero::name() const{
	return "JFuncZero";
}

bool JFuncLinear::compareRegion(RegionPorous const *r1, RegionPorous const *r2) const {
	return (r1->pd > r2->pd);
}
std::string JFuncLinear::name() const{
	return "JFuncLinear";
}

JFuncVang::JFuncVang(const double m_,const double e_):m(m_),e(e_){
	j0 = j(e);
}
bool JFuncVang::compareRegion(RegionPorous const *r1, RegionPorous const *r2) const {
	return (r1->pd > r2->pd);
}
std::string JFuncVang::name() const{
	std::stringstream ss;
	ss << "JFuncVang: m = " << m
	   << " e = " << e
	   << " j0 = " << j0 ;
	return ss.str();
}

JFuncBrooks::JFuncBrooks(const double lambda_):lambda(lambda_){
}
bool JFuncBrooks::compareRegion(RegionPorous const *r1, RegionPorous const *r2) const {
	return (r1->pd < r2->pd);
}
std::string JFuncBrooks::name() const{
	std::stringstream ss;
	ss << "JFuncBrooks: lambda = " << lambda;
	return ss.str();
}



/***************************************************************************
 * KFunc
 **************************************************************************/
std::string KFuncFirooz::name() const{
	std::stringstream ss;
	ss << "FiroozRelativePerm vw: " << vw_
	   << " vn: " << vn_
	   << " kw0: " << kw0_
	   << " kn0: " << kn0_;
	return ss.str();
}

KFuncFirooz::KFuncFirooz(const double vw, const double vn,
						 const double kw0, const double kn0):vw_(vw),vn_(vn),kw0_(kw0),kn0_(kn0){
}

std::string KFuncVang::name() const{
	std::stringstream ss;
	ss << "VangRelativePerm m: " << m_
	   << " kw0: " << kw0_
	   << " kn0: " << kn0_;
	return ss.str();
}

KFuncVang::KFuncVang(const double m, const double kw0, const double kn0):m_(m),kw0_(kw0),kn0_(kn0){
}

std::string KFuncBrooks::name() const{
	std::stringstream ss;
	ss << "BrooksRelativePerm lambda: " << lambda_
	   << " kw0: " << kw0_
	   << " kn0: " << kn0_;
	return ss.str();
}

KFuncBrooks::KFuncBrooks(const double lambda, const double kw0, const double kn0):lambda_(lambda),kw0_(kw0),kn0_(kn0){
}
