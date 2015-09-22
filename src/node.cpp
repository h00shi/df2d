/* @file node.cpp
 * .cpp file for node.
 */

#include <iostream>
#include "node.hpp"
#include "error.hpp"
#include "bvertex.hpp"

/*****************************************************************************
 * DuplData
 ****************************************************************************/
 
void DuplData::constructBase(RegionPorous *reg_){
	reg = reg_;
	idx = 0;
}

bool DuplData::operator< (const DuplData& y) const{
	//the dupldatas are sorted according to their Region::idx
	//the dupldata with highest Region::idx is the master
	return (reg->idx < y.reg->idx);
}

/*****************************************************************************
 * Node
 ****************************************************************************/

void Node::constructBase(const int idx_, const double x_, const double y_){
	idx = idx_;
	x = x_;
	y = y_;
	n_dd = 0;
	bvertex = (BVertexCQ*) NULL;
}

DuplData* Node::addRegion(RegionPorous *reg_){
	FuncBegin();

	DuplData dum;
	dum.constructBase(reg_);	
	std::list<DuplData>::iterator it;

	//search for the region in the dupldata list
	it = std::lower_bound(dd.begin(), dd.end(), dum);
	//if a dupldata is not available create one in the right place
	if ( (it == dd.end()) || (it->reg->idx != dum.reg->idx) ){	
		it = dd.insert(it, dum);
		n_dd++;
	}
	//return the associated dupldata
	return &(*it);
	
	FuncEnd();
}
