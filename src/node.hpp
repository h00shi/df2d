/** @file node.hpp
 * Header file for node and dupldata data structures.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 @ingroup mesh_module
 */

#ifndef NODE_HPP
#define NODE_HPP

#include "region.hpp"
#include <list>

// The bvertex is included via the cpp file
class BVertexCQ;

/**
 * @addtogroup mesh_module
 * @{
*/
 
/***************************************************************************
 * DuplData
 **************************************************************************/
/** @brief Handles multiple dof at a single node. Is called المثنی in the
	thesis.

	@note
	To create\n
	\li use constructBase - it is done by Node::addRegion
	\li set idx - it is done by the Mesh class
	\li To destroy destructor is enough
 */
struct DuplData{
	 /** @brief A pointer to the region, to which this DuplData belongs.
	  */
	RegionPorous *reg;
	 /** @brief Index of the DuplData.

		 The saturation and potential values are stored in the MData class.
		 This index allows each DuplData to access its relevant saturation,
		  mobilities, cappilary pressure, ... in the MData.
	  */
	int idx;           
	/** @brief Initializes the structure.
	 */
	void constructBase(RegionPorous *);
	/** @brief Sorting operator.

		As I mentioned in my thesis, at each node there can be multiple
		DuplDatas, but only one of them is the master. DuplData x will
		be less than y when y is a better candidate to be the master.
		The more the region's idx, the more the DuplData's priority.
	 */
	bool operator< (const DuplData& ) const;
};

/***************************************************************************
 * Node
 **************************************************************************/
/** @brief Data structure for a node.

   @note
   To construct use:
   \li constructBase (through Mesh class)
   \li addRegion (through Mesh through Element class)
   \li set bvertex (through BVertex)
   \li index dd (through mesh)
   \li To destroy destructor is enough
 */
struct Node{
	double x,          /**< @brief x coordinate */
		y;             /**< @brief y coordinate */
	int idx;           /**< @brief index of the node (used for access to cont edata)*/
	int n_dd;          /**< @brief number of associated DuplData s */
	/** @brief Stores the DuplDatas related to this node.

	    After the node is created, each element to whom the node belongs
		introduces itself to the node, so the node knows how many regions
		it belongs to and creates the required DuplDatas and stores them
		in this list.
	*/
	std::list<DuplData> dd;
	/** @brief The boundary vertex that the node may belong to

		If the node is located on a boundary, with permeable boundary conditions,
		this pointer will point to a BVertex object, which handles the boundary
		conditions. If not, its value will be NULL.
	*/
	BVertexCQ *bvertex; 
	/** @brief initialize the node.
	 */
	void constructBase(const int idx_, const double x_, const double y_);
	/** @brief Adds a DuplData to Node::dd list and returns its address.

		After all the elements, nodes and regions are created, each elements
		introduces itself to its member nodes through this function. This
		function does two main thing:
		\li Makes sure a DuplData is available at the node, which belongs to
		    the elements region (It will create it if the DuplData is not
			available.)
		\li Returns a pointer to the DuplData, so that by saving it the
		    elements knows the associated DuplDatas of its member nodes.
	 */
	DuplData* addRegion(RegionPorous *reg_);
};

/** @}
 */

#endif /*NODE_HPP*/
