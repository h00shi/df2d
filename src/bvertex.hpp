/** @file bvertex.hpp
 * Header file for Boundary Vertices.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 * @ingroup mesh_module
 */

#ifndef BVERTEX_HPP
#define BVERTEX_HPP

#include "node.hpp"
#include "formula.hpp"
#include <vector>
#include <petscmat.h>
#include <armadillo>
#include <string>
#include "gravity.hpp"

/***************************************************************************
 * Constant Flux (base)
 **************************************************************************/
/** @ingroup mesh_module
    @brief Class for boundary vertex with constant in(out)flow.
	
	This will be the base class and other boundary vertex types will be
	inherited from it. The BVertex classes are responsible for imposing
	boundary conditions. For each pressure(potential) boundary condition
	there is a different BVertex class. However, The saturation boundary
	condition is only controlled using an if statement and an enum.

	@note
	To construct:
	\li use constructor. (through Mesh)
	\li add neighbours. (through Mesh) 
	\li check neighbours (through Mesh through Element) 
	\li constructBase (through Mesh)
	
	To use:
	\li call assemP before solving p equation (through solver)
	\li call findQAll after solving p equation and assembling s eqution(through solver)
	\li call assemS after calling findQAll (through solver) 
	\li use mGammaW and mGamma to update overall flux (through solver)
	
	To destroy:
	\li use delete over all list members
*/
class BVertexCQ{
protected:
	double length_;     /**< @brief associated boundary length of the node */
	arma::vec2 nA_;     /**< @brief nA = l1n1 + l2n2 */
	double mGamma_,     /**< @brief Total inflow, in thesis was named: u_{\Gamma i} */
		mGammaW_;         /**< @brief Wetting phase inflow in thesis was named: u_{\Gamma W i}*/
	Node *pre_,         /**< @brief The previous node connected to bvertex (RightHandRule) */
		*self_,           /**< @brief The node that is the bvertex */
		*next_;           /**< @brief the next node connected to bvertex (RightHandRule) */
	double kdgdzna;        /**< @brief stores the value K*(gn-gw)*Grad(z)*n*A */

   /** @brief calculates the wetting phase flux.
	   
	   @param F the RHS of S equation.
	   @param Lw the value of wetting phase mobilities.
	   @param Ln the value of non-wetting phase mobilites.
   */
	void findQW(const std::vector<double> &F,
				const std::vector<double> &Lw , const std::vector<double> &Ln);
	/** @brief Finds the length	 */
	void constructL();
	/** @brief Finds the normal vector and kdgdzna */
	void constructNA(const Gravity &grav);
	/** @brief Common function shared by all constructors
		
		@param nd the self node
		@param reg the boundary region the vertex belongs to.
	*/
	void constructBase(Node *nd, RegionBoundary *reg);
	/** @brief shared description	 */
	std::string nameB() const;
	/** @brief This constructor can not be called */
	BVertexCQ(){}

public:
	RegionBoundary *reg_; /**< @brief the boudnary condition region the bvertex belongs to */
	/** @brief gives you a new BVertexCQ based on the PType enum
	 *
	 * Depending on the boundary condition, i.e. pressure constant or inflow constant
	 * you would need and BVertexCP or BVertexCQ. This function takes a region, examines
	 * its content and allocates memory to  a new BVertex depending on the region type.
	 * Then it returns the pointer to this new memory.
	 **/
	static BVertexCQ* neww(Node *nd, RegionBoundary *reg, RegionBoundary::PType ptype);
	/** @brief access protected data.
		@returns mGamma
	 */
	const double& mGamma() const  {return mGamma_;}
	/** @brief access protected data.
		 @return mGammaW
	*/
	const double& mGammaW() const {return mGammaW_;}
	/** @brief calculates both the wetting flux and the non-wetting flux.
		
		@param F the RHS of S equation.
		@param P the pointer to P values.
		@param Lw the value of wetting phase mobilities.
		@param Ln the value of non-wetting phase mobilites.
	 */
	virtual void findQAll(const std::vector<double> &F, double const * const P,
						  const std::vector<double> &Lw , const std::vector<double> &Ln);
	/** @brief adds the contribution of the vertex to the S equation.
		
		@param F the RHS of S equation.
	*/
	void assemS(std::vector<double> &F) const;
	/** @brief adds the contribution of the vertex to the P equation.
		
		@param A LHS of P equation.
		@param b RHS of P equation.
	*/
	virtual void assemP(Mat A, Vec b);
	/** @brief adds a neighbour to the vertex.
		
		@param nd the neighbour to be added.
		@returns 1 if the node was set for pre, 2 for next and 0 if both were full.
		@note checks if the neighbour is already a bvertex and checks the boundary
		regions to be the same, generates error if they are not.
	*/
	int addNeigh(Node *nd);
	/** @brief checks the orientation of next and pre by using next.
		
		@param nd the node that is supposed to be next.
		@returns 1 if next was right, 2 if it got right by swapping and 0 if node not found.
	*/
	int checkNext(Node *nd);
	/** @brief checks the orientation of next and pre by using pre.
		
		@param nd the node that is supposed to be pre.
		@returns 1 if pre was right, 2 if it got right by swapping and 0 if node not found.
	*/
	int checkPre(Node *nd);
	/** @brief creates the geometric parameters related to vertex like length and normal vector.
		
		@param dp the dimensionless p number
		@param A a global assembles matrix (the content is not important)
		@param grav the gravity field.
		@note this function must be called after addNeigh, checkPre and checkNext. 
	 */
	virtual void constructGeoParams(const double dp, Mat A, const Gravity &grav);
	/** @brief returns the description of the vertex in text form.
	 */
	virtual std::string name() const;
	/** @brief creates a bvertex and informs the node it belongs to.
		
		@param nd the self node
		@param reg the boundary region the vertex belongs to.
	*/
	BVertexCQ(Node *nd, RegionBoundary *reg);
	/** @brief Does nothing just for no-warnings.
	 */
	virtual ~BVertexCQ() {}
};

/***************************************************************************
 * Constant Flux (base)
 **************************************************************************/
/** @ingroup mesh_module
    @brief constant pressure bvertex.
	
	@note
	To construct:
	\li use constructor. (through Mesh)
	\li add neighbours. (through Mesh) 
	\li check neighbours (through Mesh through Element) 
	\li constructBase (through Mesh)
	
	To use:
	\li call assemP before solving p equation (through solver)
	\li call findQAll after solving p equation and assembling s eqution(through solver)
	\li call assemS after calling findQAll (through solver) 
	\li use mGammaW and mGamma to update overall flux (through solver)
	
	To destroy:
	\li use delete over all list members
*/
class BVertexCP : public BVertexCQ{
private:

protected:
	double p_;           /**< @brief The constant flux that should be enforced on P equation */
	arma::rowvec lhs_;   /**< @brief The lhs of P equation before forcing constant p */
	arma::irowvec conn_; /**< @brief The index of vertex neighbours */
	double rhs_;         /**< @brief The rhs of P equation before forcing constant p */
	int nConn_;          /**< @brief Number of vertex neighbours */

	/** @brief This constructor should not be called */
	BVertexCP() {}
	
public:
	void findQAll(const std::vector<double> &F, double const * P,
				  const std::vector<double> &Lw , const std::vector<double> &Ln);
	void assemP(Mat A, Vec b) ;
	void constructGeoParams(const double dp, Mat A, const Gravity &grav);
	std::string name() const;
	/** @brief creates a bvertex and informs the node it belongs to.
		@param nd the self node
		@param reg the boundary region the vertex belongs to.
	*/
	BVertexCP(Node *nd, RegionBoundary *reg);
	
};
#endif /*BVERTEX_HPP*/

