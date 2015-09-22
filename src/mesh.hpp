/** @file mesh.hpp
	header file for Mesh class.
	This class is responsible for creating, managing and deleting
	regions, elements, nodes and bvertices.

	Written by Shayan Hoshyari as a part of DF_2d
	@ingroup mesh_module
*/

#ifndef MESH_HPP
#define MESH_HPP

#include "element.hpp"
#include "bvertex.hpp"
#include "jkfunc.hpp"
#include <list>
#include <vector>
#include <petscmat.h>

/** @brief all mesh data in one class.
	@ingroup mesh_module
 */
class Mesh{
protected:
	std::list<BVertexCQ*> lbvertex_ptr_; /**< @brief polymorphic list storing bvertices */
	std::list<eleblank*> lele_ptr_;      /**< @brief polymorphic list storing elements */
	std::list<Region*> lreg_ptr_;        /**< @brief polymorphic list storing regions  */
	std::vector<Node> vnode_;            /**< @brief normal list storing nodes         */
	int ndupldata_;                      /**< @brief number of dupldata */

	/** @brief find a region from ID.
		@param regid id of the region
		@returns pointer to the region
	*/
	Region* findRegion(const int regid);
public:
	/** @brief bvertex size */
	int nbvertex() const                 
		{ return lbvertex_ptr_.size(); }
	/** @brief ele size*/
	int nele() const                     
		{ return lele_ptr_.size(); }
	/** @brief region size*/
	int nreg() const                     
		 { return lreg_ptr_.size(); }
	/** @brief node size*/
	int nnode() const                    
		{ return vnode_.size(); }
	/** @brief total number of dupldata */
	int ndd() const                      
		{ return ndupldata_; }
	
	/** @brief first bvertex */
	std::list<BVertexCQ*>::iterator begbvertex()    
		{ return lbvertex_ptr_.begin(); }
	/** @brief first ele*/
	std::list<eleblank*>::iterator begele()          
		{ return lele_ptr_.begin(); }
	/** @brief first region*/
	std::list<Region*>::iterator begreg()            
		{return lreg_ptr_.begin(); }
	/** @brief first node*/
	std::vector<Node>::iterator begnode()            
		{ return vnode_.begin(); }

	/** @brief last bvertex */
	std::list<BVertexCQ*>::iterator endbvertex()   
		{ return lbvertex_ptr_.end(); }
	/** @brief last ele */
	std::list<eleblank*>::iterator endele()        
		{ return lele_ptr_.end(); }
	/** @brief last region */
	std::list<Region*>::iterator endreg()          
		{return lreg_ptr_.end(); }
	 /** @brief last node */
	std::vector<Node>::iterator endnode()         
		{ return vnode_.end(); }   

	/** @brief change the reserved size of node vector.
		
		@param sz the new reserved size
	 */
	void reserveNodes(const int sz);
	/** @brief add a node.
		
		@param x x co-ordinate
		@param y y co-ordinate
		@param idx index of the node
		@note indices should be added in order, and from zero
	*/
	void addNode(const double x, const double y, const int idx);
	/** @brief add an element or bvertex.
		
		@param regid ID of the region the element belongs to
		@param ndidx index of corner nodes
		@param celltype type of the cell
		@param linenumber number of the line the node was read from
		       used to generate error if anything goes wrong. is optional.
	*/
	void addElement(const int regid, const int ndidx[], const CellType celltype, const int linenumber= -1);
	/** @brief add a region.
		
		@param regionpointer pointer to the newly added region.
	*/
	void addRegion(Region *regionpointer);
	/** @brief pre-process everything.
		
	   @param cmp defines how to sort regions. use JFunc::cmp.
	   @param dp dimensionless p number
	   @param A a petsc matrix, the matrix should be null before call,
	   after the call the matrix will be fully defined.
	*/
	void constructGeoParams(const RegionPointerComparer& cmp,
							const double dp, Mat &A, const Gravity &grav);
	/** @brief sets ndupldata_ to zero */
	Mesh();
	/** @brief frees memory */
	~Mesh();

};


#endif /*MESH_HPP*/
