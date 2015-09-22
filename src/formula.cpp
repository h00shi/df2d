/** @file formula.cpp
 * The functions that hold fundamental formulas are here;
 */

#include "formula.hpp"
#include <cmath>

namespace fml{

	void chkIdx(const int i , const int I){
		FuncBegin();
		if ( (i >= I) && (I > 0) ) {
			Error::mess << "index bigger than or equal to bound, index: " << i << " bound: " << I;
			ERRSET();
		}
		FuncEnd();
	}

	/** @cond */
	arma::mat22 Rot ("0,-1;1,0");
	/** @endcond */

	double lineLength(const double x1,const double y1,
					  const double x2, const double y2){
		return sqrt(pow(x1-x2,2) + pow(y1-y2,2));
	}

	void lineNA(const double x1,const double y1,
				const double x2, const double y2,
				arma::vec2 &ans){
		ans.resize(2);
		ans(0) = x2 - x1;
		ans(1) = y2 - y1;
		ans = -Rot * ans;
	}

	double findHeight (const double x, const double y, const Gravity &g){
		double norm_direc;
		norm_direc = sqrt( g.xdir*g.xdir + g.ydir*g.ydir);
		return ( g.xdir * (x-g.x0) + g.ydir * (y - g.y0) ) / norm_direc ;
	}

}
