/** @file geom.hpp
 * Header file for polygons, used for changing field values in setfield.
 * Written by Shayan Hoshyari
 * As a part of DF_1d
 * @ingroup dr_module
 */

#ifndef GEOM_HPP
#define GEOM_HPP

#include <cmath>
#include <string>
#include <sstream>

/** @brief Abstract class for a polygon.
	@ingroup dr_module

	This class and it's inherited classes store the data for a geometrical
	shape. Then they tell you if a point is inside or outside them. They
	are used for changing the initial conditions according to @setfield.	
 */
class Poly{
public:
	/** @brief virtual destructor for removing warnings */
	virtual ~Poly(){}
	/** @brief returns the description of the polygon */
	virtual std::string name() const = 0;
	/** @brief tells if a point is within the polygon */
	virtual bool isin(const double x , const double y) const = 0;
};

/** @brief Class for a rectangle.
	@ingroup dr_module
 */
class Rectangle : public Poly{
protected:
	double x0, /**< @brief lower left point x coordinate */
		y0,      /**< @brief lower left point y coordinate */
		x1,      /**< @brief upper right point x coordinate */
		y1;      /**< @brief upper right point y coordinate */
public:
	/** @brief sets the lower-left and upper-right points.
	* 
	* The orientation of the points:	
	*       X-------------------x1\n
	*       |||||||||||||||||||||\n
	*       |||||||||||||||||||||\n
	*       x0------------------X
	 */
	Rectangle(const double x0_,const double y0_,
			  const double x1_, const double y1_);
	std::string name() const;
	bool isin(const double x , const double y) const;
};


/** @brief Class for a circle.
	@ingroup dr_module
 */
class Circle : public Poly{
protected:
	double x0, /**< @brief center x coordinate */
		y0,      /**< @brief  center y coordinate */
		r;      /**< @brief radius */
public:
	/** @brief sets center co-ordinates and radius.
			@param x0_ x coordinate of center
			@param y0_ y coordinate of center
			@param r_ radius
	 */
	Circle(const double x0_,const double y0_, const double r_);
	std::string name() const;
	bool isin(const double x , const double y) const;
};



 
#endif /*GEOM_HPP*/
