/** @file gravity.hpp
 *  The structure defining the gravity field.
 *  Written by Shayan Hoshyari
 * @ingroup edat_module
 */

#ifndef GRAVITY_HPP
#define GRAVITY_HPP

/** @ingroup edat_module
 * @brief stores the data for gravity field.	
*/
struct Gravity{
	double x0,  /**< @brief x coordinate of a point with h = 0 */
		y0,     /**< @brief y coordinate of a point with h = 0 */
		gw,     /**< @brief dimless number rho_n*g*L/P    */
		gn,     /**< @brief dimless number rho_w*g*L/P    */
		xdir,     /**< @brief x direction of gravity vector */
		ydir;     /**< @brief y direction of gravity vector */
};

#endif /*GRAVITY_HPP*/
