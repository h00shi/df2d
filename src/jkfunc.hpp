/** @file jkfunc.hpp
 * Classes concerning cappilary pressure and relative permeability curves.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 * @ingroup edat_module
 */


#ifndef JKFUNC_HPP
#define JKFUNC_HPP

#include <cmath>
#include <string>
#include <sstream>

// The file region.hpp is included in the .cpp file to prevent mutual including
// of header files. So these lines have to be added here.
class Region;
class RegionPorous;
class JFunc;

/*****************************************************************************
 *RegionPointerComparer
 *
 ****************************************************************************/
/** @ingroup edat_module
	@brief A functor to sort regions.

	As I mentioned in my thesis the regions should be sorted according to their
	priority to be master regions. As it is up to JFunc to decide how the regions
	have to be sorted, the JFunc has a functor member class. It sends the functor
	member class to STL sorting functions. This class is the functor.
 */
class RegionPointerComparer{
protected:
	JFunc const *master_; /**< @brief The JFunc this class belongs to */
public:
	/** @brief compare function. */
	bool operator() (Region const*, Region const*) const;

	/** @brief gives value to master_. */
	RegionPointerComparer(JFunc const*);
};

/*****************************************************************************
 *JFunc
 *
 ****************************************************************************/

/** @ingroup edat_module
	@brief Abstract class for the J curve.
	
	As we have multiple J curves, this is an abstract class that stands for all
	of them. Using polymorphism a pointer to this class, can point to the concrete
	curves.
	
	@note To create one using the constructor is sufficient.

 */
class JFunc{
public:
	/** @brief The comparer specific for this class.
	 *
	 * Any JFunc sends this functor to STL sorting functions. The functor's sole
	 * responsibility is to call the compareRegion function.
	 */
	RegionPointerComparer cmp; 
	/** @brief function that indicates how regions have to be compared.
		
		@param r1 first pointer
		@param r2 second pointer
		@return true if r1 < r2 and false otherwise
		@note reference region is the smallest of all for firooz model and the
		largest of all for brooks model.
	 */
	virtual bool compareRegion(RegionPorous const* r1, RegionPorous const* r2) const = 0;
	/** @brief returns value of J.
	 *
	 * @param s saturation
	 * @note that J * pd = pc
	 */
	virtual double J(const double s) const = 0;
	/** @brief returns ds_slave/ds_master = f(s,p_master,p_slave).
	 *
	 * @param s is s
	 * @param p1 is p_master
	 * @param p2 is p_slave
	 */
	virtual double ds(const double s, const double p1 ,const double p2) const = 0;
	/** @brief returns s_slave = f(s,p_master,p_slave).
	 *
	 *	@param s is s
	 *	@param p1 is p_master
	 *	@param p2 is p_slave
	*/
	virtual double sopp(const double s, const double p1, const double p2) const = 0;
	/** @brief returns description of the curve as readable string.
	 */
	virtual std::string name() const = 0;
	/** @brief destructor for removing warning.
	 */
	virtual ~JFunc() {};
protected:
	/** @brief Initializes the comparer.
	 */
	JFunc();	
};

/** @ingroup edat_module
 * @brief Firooz J curve.
 *	
 *	J(s) = -ln(s) and regions are sorted with biggest pd first.
 */
class JFuncFirooz : public JFunc{
public:
	/**@brief big pd must be first */
	bool compareRegion(RegionPorous const *r1, RegionPorous const *r2) const;
	inline double J(const double s) const{
		return -log( fmax(s,.001) );
	}
	inline double ds(const double  s, const double p1, const double p2) const{
		double r = p1/p2;
		return r * pow(s , r - 1);
	}
	inline double sopp(const double s, const double p1, const double p2) const{
		return pow(s , p1/p2);
	}
	std::string name() const;
};

/** @brief Zero cappilary pressure.
 *
 * J(s) = 0 .
 */
class JFuncZero : public JFunc{
public:
	/**@brief priority is unimportant */
	bool compareRegion(RegionPorous const *r1, RegionPorous const *r2) const;
	inline double J(const double s) const{
		return 0;
	}
	inline double ds(const double s, const double p1, const double p2) const{
		return 1;
	}
	inline double sopp(const double s, const double p1, const double p2) const{
		return s;
	}
	std::string name() const;
};

/** @brief linear cappilary curve.
 *
 * J(s) = 1-s
*/
class JFuncLinear : public JFunc{
public:
	/**@brief big pd must be first */
	bool compareRegion(RegionPorous const *r1, RegionPorous const *r2) const;
	inline double J(const double s) const{
		return (1-s);
	}
	inline double ds(const double s, const double p1, const double p2) const{
		return (s < 1 - p2/p1 ? 0 : p1/p2 );
	}
	inline double sopp(const double s, const double p1, const double p2) const{
		return (s < 1 - p2/p1 ? 0 : 1 - p1/p2*(1-s) );
	}
	std::string name() const;
};

/** @ingroup edat_module
 *  @brief VanGenuchten J curve.
 *
 * for this kind of curve we have: 
 * \li j = (s .^ (-1 / m) - 1) .^ (1 - m) 
 * \li sminus = ( r ^ (1/(1-m)) * (s .^ (-1/m) - 1) + 1 ) .^ (-m) 
 * \li jminus = (1+ (y .^ (1/(1-m)))) .^ (-m) 
 * \li dsminus = (r ^ (1/(1-m)) *(s .^ (-1/m)-1)+1).^(-1-m) * r ^ (1/(1-m)) .* s .^ (-1/m - 1)
 * 
 * @note This curve is not tested enough and has bugs. 
 * because these curves might produce NaN numbers near s=0 and s=1. 
 * currently I manually use a fmin(fmax()) function to bound s between 1+epsilon 1-epsilon 
 * but this is not a good solution. interpolation might be a better solution. 
 * so please use with caution.
 */
class JFuncVang : public JFunc{
protected:
	double m, /**< @brief m parameter of vangenuchten curve */
		j0,   /**< @brief J(epsilon)*/
		e;    /**< @brief truncation parameter epsilon */
	/** @brief finds J(s) without NaN checking */
	inline double j(const double s) const{
		return pow( pow(s,-1/m) - 1 , 1-m );
	}
	/** @brief finds J(RJinv(s)) but without NaN checking */
	inline double sminus(const double s, const double r) const{
		return pow( pow(r, (1/(1-m))) * ( pow(s, -1/m) - 1 ) + 1  , -m);
	}
	/** @brief finds Jinv(s) but without NaN checking */
	inline double jminus(const double y) const{
		return pow( 1 + pow(y,1/(1-m)) , -m);
	}
	/** @brief finds d J(RJinv(s))/ds but without NaN checking */
	inline double dsminus(const double s,const double r) const{
		return pow( pow(r, 1/(1-m)) * (pow(s,-1/m)-1) + 1 , (-1-m) ) * pow(r, 1/(1-m)) * pow(s, -1/m - 1);
	}
public:
	/** @brief gives value to m,e and j0 */
	JFuncVang(const double m_,const double e_ = .01);	
	/** @brief big pd must be first */
	bool compareRegion(RegionPorous const *r1, RegionPorous const *r2) const;
	inline double J(const double s) const{
		return ( s < e ? j0 : j( fmax(0.001,fmin(s,.999)) )  );
	}
	inline double ds(const double  s, const double p1, const double p2) const{
		double r, sl;
		r = p1/p2;
		sl = jminus(j0/r);
		return ( s < sl ? 0 : dsminus(fmax(0.001,fmin(s,.999)),r) ); 
	}
	inline double sopp(const double s, const double p1, const double p2) const{
		double r, sl;
		r = p1/p2;
		sl = jminus(j0/r);
		return ( s < sl ? 0 : sminus(fmax(0.001,fmin(s,.999)),r) ); 
	}
	std::string name() const;
};

/** @ingroup edat_module
	@brief Brooks-Corey J curve.
 */
class JFuncBrooks : public JFunc{
protected:
	double lambda; /**< @brief sole parameter of brooks-corey curves */
public:
	/**@brief  gives value to lambda */
	JFuncBrooks(const double lambda_);	
	/** @brief small pd must be first */
	bool compareRegion(RegionPorous const *r1, RegionPorous const *r2) const;
	inline double J(const double s) const{
		return pow( fmax(0.01,s)  ,-1/lambda);
	}
	inline double ds(const double  s, const double p1, const double p2) const{
		double r, sr;
		r = p1/p2;
		sr = pow(r,lambda);
		return ( s > sr ? 0 : pow(r,-lambda) ); 
	}
	inline double sopp(const double s, const double p1, const double p2) const{
		double r, sr;
		r = p1/p2;
		sr = pow(r,lambda);
		return ( s > sr ? 1 : pow(r,-lambda)*s ); 
	}
	std::string name() const;
};

/*****************************************************************************
 *KFunc
 *
 ****************************************************************************/

/** @ingroup edat_module
	@brief Abstract class for relative permeability curves.

	@note To construct one the constructor is sufficient.
 */
class KFunc{
public:
	/** @brief returns k_rw.
	 * @param s saturation
	 */
	virtual  double w(const double s) const = 0;
	/** @brief returns k_rnw
	 * @param s saturation
	 */
	virtual  double nw(const double s) const = 0;
	/** @brief returns the name of the model.
	 */
	virtual std::string name() const = 0;
	/** @brief destructor.
	 */
	virtual ~KFunc() {};
};

/**  @ingroup edat_module
 *  @brief Firooz relative permeability curves.
 *
 * In this curves the relative permeabilities are found from: 
 * \li krw(S) = kw0(S)^vw 
 * \li krn(S) = kn0(S)^vn
 */
class KFuncFirooz : public KFunc{	
protected:
	double vw_, /**< @brief exponent for wetting phase relative permeability function */
		vn_, /**< @brief exponent for non-wetting phase relative permeability function */
		kw0_, /**< @brief coefficient for wetting phase relative permeability function */
		kn0_; /**< @brief coefficient for non-wetting phase relative permeability function */
public:	
	inline double w(const double s) const{
		return kw0_*pow(s,vw_);
	}
	inline double nw(const double s) const{
		return kn0_*pow(1-s,vn_);
	}
	std::string name() const;
	/**  sets vw,vn,kw0 and kn0.
	 */
	KFuncFirooz(const double vw, const double vn,
				const double kw0, const double kn0);	
};

/**  @ingroup edat_module
 * @brief VanGenuchten - Parker relative permeability curves.
 *
 * \li In this curves the relative permeabilities are found from: 
 * \li krw(S) = kw0 sqrt(S)   ( 1-( 1-S^(1/m) )^m )^2 
 * \li krn(S) = kn0 sqrt(1-S) ( 1-S^(1/m) ) ^ 2m
 */
class KFuncVang : public KFunc{	
protected:
	double m_, /**< @brief exponent for relative permeability functions */
		kw0_, /**< @brief coefficient for wetting phase relative permeability function */
		kn0_; /**< @brief coefficient for non-wetting phase relative permeability function */
public:	
	inline double w(const double s) const{
		double ss = fmax(0.001,fmin(s,.999));
		return kw0_ * sqrt(ss) * pow( 1 - pow(1- pow(ss,1/m_), m_), 2);
	}
	inline double nw(const double s) const{
		double ss = fmax(0.001,fmin(s,.999));
		return kn0_ * sqrt(1-ss) * pow( 1 - pow(ss,1/m_), 2*m_ );
	}
	std::string name() const;
	/**  sets m, kw0 and kn0.
	 */
	KFuncVang(const double m,	const double kw0, const double kn0);	
};

/** @ingroup edat_module
 *  @brief BrooksCorey  relative permeability curves.
 * In this curves the relative permeabilities are found from: 
 * \li krw(S) = kw0 s^(3+2*lambda) 
 * \li krn(S) = kn0 (1-s)^2 * ( 1 - s^(1+2*lambda) )
 */
class KFuncBrooks : public KFunc{	
protected:
	double lambda_, /**< @brief exponent for relative permeability functions */
		kw0_, /**< @brief coefficient for wetting phase relative permeability function */
		kn0_; /**< @brief coefficient for non-wetting phase relative permeability function */
public:	
	inline double w(const double s) const{
		return kw0_ * pow(s,3+2*lambda_);
	}
	inline double nw(const double s) const{
	    return kn0_ * pow(1-s,2) * (1 - pow(s,1+2*lambda_));
	}
	std::string name() const;
	/**  @brief sets m, kw0 and kn0.
	 */
	KFuncBrooks(const double lambda,	const double kw0, const double kn0);	
};


#endif /*JKFUNC_HPP*/
 
