/** @file geom.hpp
 * cpp file for polygons, used for changing field values in setfield.
 * Written by Shayan Hoshyari
 * As a part of DF_1d
 */

#include "geom.hpp"
#define EPSILON .00000000001

Rectangle::Rectangle(const double x0_,const double y0_,
					 const double x1_, const double y1_):
	x0(x0_-EPSILON),y0(y0_-EPSILON),x1(x1_+EPSILON),y1(y1_+EPSILON){}
	
std::string Rectangle::name() const{
	std::stringstream ss;
	ss << "Rectangle x0-y0-x1-y1: "
	   << x0 << " "<< y0 << " "<< x1 << " "<< y1;
	return ss.str();
}
	
bool Rectangle::isin(const double x , const double y) const{
	return ( ((x-x0) > 0) && ((x-x0) < (x1-x0)) &&
			 ( (y-y0)>0 ) && ((y-y0) < (y1-y0)) );
}

Circle::Circle(const double x0_,const double y0_,
			   const double r_):
	x0(x0_),y0(y0_),r(r_+EPSILON){}
	
std::string Circle::name() const{
	std::stringstream ss;
	ss << "Circle x0-y0-r: " << x0 << " "<< y0 << " " << r;
	return ss.str();
}
	
bool Circle::isin(const double x , const double y) const{
		return ( ( pow(x-x0,2)  + pow(y-y0,2) ) < pow(r,2) );
}


