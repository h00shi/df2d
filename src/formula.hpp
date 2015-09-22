/** @file formula.hpp
 * The functions that hold fundamental formulas are here;
 * @ingroup mesh_module
 */

#ifndef FORMULA_HPP
#define FORMULA_HPP

#include <armadillo>
#include "gravity.hpp"
#include "error.hpp"


/** @ingroup mesh_module
	@namespace fml
	@brief Holds the functions for fundamental formulas.
*/
namespace fml{
	/** @brief checks i < I.
		
		@param i the index to be checked
		@param I the maximum amount
		@note Actually the real purpose of this function was to be able to
		enable and disable bound checking using a compiler switch, which is
		not implemented up to now.
	 */
	void chkIdx(const int i , const int I);
	/** @brief creates local data from global data.
		@param n_glob the maximum size of glob
		@param glob the container of global data
		@param n_loc the number of local data
		@param idx the local indices
		@returns loc the local data
		@note this func does not allocate data, loc should already have data.
	 */
	template <typename G, typename L>
	inline void locFromGlob(const int n_glob, const G &glob,
							const int n_loc, const L *loc,
							arma::vec& ans){
		FuncBegin();
		for (int i = 0 ; i < n_loc ; i++){
			chkIdx(loc[i]->idx, n_glob);
			ans(i) = glob[loc[i]->idx];
		}
		FuncEnd();
	}
	/** @var arma::mat22 fml::Rot
		@brief 90 degrees rotation vector.
	*/
	extern arma::mat22 Rot;	
	/** @brief Finds the length of a line.
		@returns ans the length of a line from r1 to r2.
	*/
	double lineLength(const double x1,const double y1,
					  const double x2, const double y2);
	/** @brief Finds the normal vector of a line.
		@note point two must be after point one according to RRR.
		@returns ans the normal vector
	*/
	void lineNA(const double x1,const double y1,
				const double x2, const double y2,
				arma::vec2 &ans);
	/** @brief Finds the height of a point */
	double findHeight (const double x, const double y, const Gravity &g);
}

#endif /*FORMULA_HPP*/
