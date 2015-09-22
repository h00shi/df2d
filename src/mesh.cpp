/** @file mesh.cpp
    defenition of Mesh class functions.

	Written by Shayan Hoshyari as a part of DF_2d
*/

#include "mesh.hpp"

using std::list;
using std::vector;


Region* Mesh::findRegion(const int regid){ 
	FuncBegin();
	for ( list<Region*>::iterator i = begreg() ; i != endreg() ; i++){
		if (regid == (*i)->ID) return (*i);
	}
	return (Region*) NULL;
	FuncEnd(); 
} 

void Mesh::reserveNodes(const int sz){ 
	FuncBegin(); 
	vnode_.reserve(sz);
	FuncEnd(); 
} 

void Mesh::addNode(const double x, const double y, const int idx){ 
	FuncBegin();
   	Node newnode;
	
	if ( nnode() == 0){
		if (idx != 0 ) {
			Error::mess << "first node should have index = 0";
			ERRSET();
		}
	}
	else if (idx-1 != vnode_.back().idx ){
		Error::mess << "nodes should be added in order, " << vnode_.back().idx << " -> " << idx;
		ERRSET();
	}
	vnode_.push_back(newnode);
	vnode_.back().constructBase(idx,x,y);
		
	FuncEnd(); 
} 

void Mesh::addElement(const int regid, const int ndidx[], const CellType celltype, const int linenumber){ 
	FuncBegin();
	Node *ndptr[10];
	Region *regptr;

	//find the region
	regptr = findRegion(regid);
	if (!regptr) return;

	//find the nodes
	for (int i = 0 ; i < (int)celltype ; i++)
		ndptr[i] = &vnode_.at(ndidx[i]);
	
	//add the element
	if ( regptr->isBoundary() ){
		RegionBoundary *regptrbnd = (RegionBoundary*) regptr;
	    switch (celltype){
			case CellPoint:
				lbvertex_ptr_.push_back( BVertexCQ::neww(ndptr[0], regptrbnd, regptrbnd->ptype) );
				break;
			case CellLine:
				if (!ndptr[0]->bvertex) 
					lbvertex_ptr_.push_back( BVertexCQ::neww( ndptr[0], regptrbnd , regptrbnd->ptype) );
				ndptr[0]->bvertex->addNeigh(ndptr[1]);
				if (!ndptr[1]->bvertex)
					lbvertex_ptr_.push_back( BVertexCQ::neww( ndptr[1], regptrbnd , regptrbnd->ptype) );
				ndptr[1]->bvertex->addNeigh(ndptr[0]);
				break;
			default:
				Error::mess << "Element belongs to region: " << regid << " which is boundary. "
							<< "However it is neither a line nor a point. It is a " << celltype
							<< " at line "<< linenumber << " of mesh file. " ;
				ERRSET();
				break;
			}
	}
	else{
		RegionPorous *regptrpour = (RegionPorous*) regptr;
		lele_ptr_.push_back( eleblank::neww(regptrpour, ndptr, celltype) );
	}
									 
	FuncEnd(); 
} 

void Mesh::addRegion(Region *regionpointer){ 
	FuncBegin(); 
	lreg_ptr_.push_back(regionpointer);
	FuncEnd(); 
} 

void Mesh::constructGeoParams(const RegionPointerComparer& cmp,
							  const double dp, Mat &A, const Gravity &grav){ 
	FuncBegin();
	int j,jbup;
	arma::vec::fixed<20> matloc;

	// sort regions
	lreg_ptr_.sort(cmp);
	j = 0;
	for ( list<Region*>::iterator i = begreg() ; i != endreg() ; i++){
		(*i)->idx = j;
		j++;
	}

	//create dupldata
	for ( list<eleblank*>::iterator i = begele() ; i != endele() ; i++){
		(*i)->constructGeoParams();
		(*i)->constructDuplData();
		(*i)->constructBVertices(NULL);
	}
	j = 0;
	jbup = 0;
	for (vector<Node>::iterator i = begnode() ; i < endnode() ; i++){
		jbup += i->n_dd;
		for (list<DuplData>::iterator k = i->dd.begin() ; k != i->dd.end() ; k++){
			k->idx = j;
			j++;
		}
	}
	ndupldata_ = j;
	if ( jbup != j ){
		Error::mess << "jbup != j" ;
		ERRSET();
	}

	//create the matrix
	Error::code=MatCreateSeqAIJ(PETSC_COMM_SELF,nnode(),nnode(),0,NULL,&A);ERRCHK();
	Error::code=MatZeroEntries(A);ERRCHK();
	Error::code=MatSetOption(A, MAT_ROW_ORIENTED, PETSC_FALSE);ERRCHK();
	Error::code=MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);ERRCHK();
	j = 1;
	for (list<eleblank*>::iterator i = begele() ; i != endele() ; i++){
		matloc.ones();
		matloc *= j;
		Error::code=MatSetValues(A ,(*i)->nNode(),(*i)->idxGlob().memptr()
								 ,(*i)->nNode(),(*i)->idxGlob().memptr()
								 ,matloc.memptr(),ADD_VALUES);ERRCHK();
		j++;
	}
	Error::code=MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);ERRCHK();
	Error::code=MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);ERRCHK();
	Error::code=MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS , PETSC_FALSE);ERRCHK();
	Error::code=MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR , PETSC_TRUE);ERRCHK();
	Error::code=MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);ERRCHK();

	//create the bvertices
	for ( list<BVertexCQ*>::iterator i = begbvertex() ; i != endbvertex() ; i++)
		(*i)->constructGeoParams(dp, A, grav);
	
	FuncEnd(); 
} 

Mesh::Mesh():ndupldata_(0){}

Mesh::~Mesh(){ 
	FuncBegin(); 
	for ( list<BVertexCQ*>::iterator i = begbvertex() ; i != endbvertex() ; i++)
		if (*i) delete (*i);
	for ( list<eleblank*>::iterator i = begele() ; i != endele() ; i++)
		if (*i) delete (*i);
	for ( list<Region*>::iterator i = begreg() ; i != endreg() ; i++)
		if (*i) delete (*i);
	FuncEnd(); 
} 
