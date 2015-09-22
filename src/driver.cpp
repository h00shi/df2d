/** @file driver.cpp
	defenitions for driver.hpp

	Shayan Hoshyari - DF_2d
	@ingroup dr_module
*/

#include "driver.hpp"
#include "asciifile.hpp"
#include "geom.hpp"
#include "visit_writer.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <ctime>

using std::string;
using std::cout;
using std::stringstream;
using std::fstream;
using std::endl;
using std::vector;
using std::list;
using std::setprecision;
using std::setw;
using std::left;
using std::setw;
using std::fixed;

/************************************************************************
 * addresses and switches
 ************************************************************************/

static const char swidir[] = "-d";
static const char swisetf[] = "-s";
static const char help[] = "DF_2d written by Shayan Hoshyari \n Please consult the doxygen documentation.";

static const string adrsol = "solver.config";
static const string adrpetsc = "petsc.config";
static const string adrmesh = "mesh";
static const string adrini = "initial";
static string adrrestart = "restart/restart";
static string adrresult = "result/result";

/************************************************************************
 * converting cell types
 ************************************************************************/

/** @brief convert gmsh int to CellType */
static CellType  gmsh2cell(const int c){
	FuncBegin();
	switch(c){
	case 1:
		return CellLine;
	case 2:
		return CellTri;
	case 3:
		return CellQuad;
	case 15:
		return CellPoint;
	default:
		ERRSET();
		return CellPoint;
	}
	FuncEnd();
}
/** @brief convert CellType to gmsh int */
static int cell2gmsh(const CellType c){
	FuncBegin();
	switch(c){
	case CellPoint:
		return 15;
	case CellLine:
		return 1;
	case CellTri:
		return 2;
	case CellQuad:
		return 3;
	default:
		ERRSET();
		return 0;
	}
	FuncEnd();
}
/** @brief convert CellType to visit int */
static int cell2visit(const CellType c){
	switch(c){
	case CellLine:
		return VISIT_LINE;
	case CellTri:
		return VISIT_TRIANGLE;
	case CellQuad:
		return VISIT_QUAD;
	default:
		ERRSET();
		return 0;
	}
}
/* @brief convert CellType to tecplot */
static int cell2tecplot(const CellType c){
	ERRSET();
}

/************************************************************************
 * reading stuff
 ************************************************************************/

/** @brief reads a kr from one line of a file */
static KFunc* readregion_kr(AsciiFile &cfile){
	FuncBegin();
	string type;
	double kw0,kn0;
	
    cfile(type, "kr_type");
	if ( type.compare("firooz") == 0 ){
		double vw,vn;
		cfile (vw, "vw_value");
		cfile (vn, "vn_value");
		cfile (kw0, "krw0_value");
		cfile (kn0, "krn0_value");
		return new KFuncFirooz(vw,vn,kw0,kn0);
	}
	else if ( type.compare("vang") == 0 ){ 
		double m;
		cfile (m, "m_value");
		cfile (kw0, "krw0_value");
		cfile (kn0, "krn0_value");
		return new KFuncVang(m,kw0,kn0);
	}
	else if ( type.compare("brooks") == 0 ){ 
		double lambda;
		cfile (lambda, "lambda_value");
		cfile (kw0, "krw0_value");
		cfile (kn0, "krn0_value");
		return new KFuncBrooks(lambda,kw0,kn0);
	}
	else {
		Error::mess << type << " not supported. "
					<< cfile.fn << ", line: " <<cfile.ln;
		ERRSET();
	}

	FuncEnd();
}

/** @brief  read a triangle mesh */
static void readmesh_triangle(MData& md,Mesh &msh){
	FuncBegin();

	AsciiFile fl;
	int size,idx0,ndidx[3],offset,bndmarker;
	double x,y;

	/*
	  read the node file

	  the format is
	  1: <# of vertices> <dim=2> <# atributes> <# boundary marker>
	  other: <vertex #> <x> <y> [attributes] [boundary marker]	  
	*/
	fl.open(md.dir+adrmesh+".node");
	fl(); fl(size,"<# of vertices>");
	msh.reserveNodes(size);
	fl(); fl(idx0,"<vertex #>"); fl(x,"x"); fl(y,"y");
	offset = idx0;
	msh.addNode(x, y, idx0-offset);
	for (int i = 1 ; i < size ; i++){
	    fl(); fl(idx0,"<vertex #>"); fl(x,"x"); fl(y,"y");
		msh.addNode(x, y, idx0-offset);
	}
	fl.close();

	/*
	  read the element file - triangles

	  the format is
	  1: <# of triangles> <nodes per triangle> <# of attributes>
	  other: <triangle #> <node> <node> <node> ... [attributes] 
	*/
	fl.open(md.dir+adrmesh+".ele");
	fl(); fl(size,"<# of triangles>");
	for (int i = 0 ; i < size ; i++){
	    fl(); fl(idx0,"<triangle #>");
		fl(ndidx[0],"<node1>"); fl(ndidx[1],"<node2>"); fl(ndidx[2],"<node3>");
		fl(bndmarker,"[attribute0]");
		ndidx[0]-=offset; ndidx[1]-=offset; ndidx[2]-=offset;
		msh.addElement(bndmarker,ndidx, CellTri, fl.ln);
	}
	fl.close();

	/*
	  read the poly file - fractures and boundary

	  the format is
	  1: <# of triangles> <nodes per triangle> <# of attributes>
	  other: <triangle #> <node> <node> <node> ... [attributes]
	  ---
	   1: <# of segments> <# of boundary markers (0 or 1)>
	  other: <segment #> <endpoint> <endpoint> [boundary marker] 
	*/
	fl.open(md.dir+adrmesh+".poly");
	fl(); fl(size,"<# of vertices>");
	for (int i = 0 ; i < size ; i++) fl();
	fl(); fl(size,"<# of segments>");
	for (int i = 0 ; i < size ; i++){
		fl(); fl(idx0,"<segment #>");
		fl(ndidx[0],"<start point>"); fl(ndidx[1],"<end point>"); 
		fl(bndmarker,"[bndmarker]");
		ndidx[0]-=offset; ndidx[1]-=offset; 
		msh.addElement(bndmarker,ndidx, CellLine, fl.ln);
	}
	fl.close();
	
	FuncEnd();
}

/** @brief  read a gmsh mesh */
static void readmesh_gmsh(MData& md,Mesh &msh){
	FuncBegin();

	AsciiFile fl;
	int size,idx0,ndidx[5],offset,bndmarker,sizetmp,bndmarkertmp, eletype;
	CellType celltype;
	double x,y;

	fl.open(md.dir+adrmesh+".msh");
	/*
	  read the nodes

	  the format is
	  * $Nodes
	  * no. of nodes
	  * idx x y z
	*/
	fl.efind("$Nodes");
	fl(); fl(size,"no. of nodes");
	msh.reserveNodes(size);
	fl(); fl(idx0,"idx"); fl(x,"x"); fl(y,"y"); 
	offset = idx0;
	msh.addNode(x, y, idx0-offset);
	for (int i = 1 ; i < size ; i++){
		fl(); fl(idx0,"idx"); fl(x,"x"); fl(y,"y"); 
		msh.addNode(x, y, idx0-offset);
	}

	/*
	  read the elements

	  the format is
	  * $Elements
	  * # of elements
	  * idx <type> <# of tags>  [tags-phys first] [nodes]
	*/
	fl.goto_beg();
	fl(); fl.efind("$Elements");
	fl(); fl(size,"<# of elements>");
	for (int i = 0 ; i < size ; i++){
		fl();
		fl(idx0,"idx");
		fl(eletype, "eletype"); celltype=gmsh2cell(eletype);
		fl(sizetmp, "tag size");
		fl(bndmarker, "physical tag");
		for (int j = 1 ; j < sizetmp ; j++)
			fl(bndmarkertmp,"unused tag");
		for(int j = 0 ; j < (int)celltype ; j++){
			fl(ndidx[j],"node id");
			ndidx[j]-=offset;
		}
		msh.addElement(bndmarker, ndidx,celltype, fl.ln);
	}
	fl.close();

	FuncEnd();
}

/************************************************************************
 * writing stuff
 ************************************************************************/
/** @brief writes a .vtk file.
 *
 * although this function allocates a lot of memmory it will not affect the
 * execution time much, as it will be called only a few times.
 */
static void writevisual_vtk	(const char * address,MData& md,Mesh &msh){
	FuncBegin();

	int usebinary = 1; 
	int npts;
	float *pts;
	int *celltypes;
	int *conn;
	int nvars = 2;
	int vardim[] = { 1,1, 1,1,1, 3,3,3};
	int centering[] = {0,1, 1,0,1, 0,0,0};
	int i,j;

	const char name1[] = "porous_region_no";
	const char name2[] = "bnd_region_no";
	const char name3[] = "Sw_node";
	const char name4[] = "Sw_cell";
	const char name5[] = "Pw";
	const char name6[] = "qw";
	const char name7[] = "qn";
	const char name8[] = "qtotal";
	const char *varnames[] = {name1,name2,name3,name4,name5,name6,name7,name8};

	float *vars[8];

	const arma::ivec *idxele;
	const arma::mat *kdarma;
	const arma::vec *parma,*pcarma,*lnarma,*lwarma;
	arma::vec2 vt , vn, vw;

	//give number of points
	npts = (md.visualduplicate ==0 ? msh.ndd() : msh.nnode());

	//allocate memory to mesh_node data
	if (md.visualduplicate == 0 ){
		vars[1] = new float [ msh.ndd() ];  //boundary region no
		pts = new float [ 3 * msh.ndd() ];  //points
		for (vector<Node>::iterator it = msh.begnode() ; it != msh.endnode() ; it++){
			for (list<DuplData>::iterator l = it->dd.begin() ; l != it->dd.end() ; l++){				
				pts[3*l->idx] = (float)it->x;
				pts[3*l->idx+1] = (float)it->y;
				pts[3*l->idx+2] = 0;
				vars[1][l->idx] = ( it->bvertex ?
									(float) it->bvertex->reg_->ID :
									-1 );
			}
		}
	}
	else{
		vars[1] = new float [ msh.nnode() ];  //boundary region no
		pts = new float [ 3 * msh.nnode() ];  //points
		i = 0;
		for (vector<Node>::iterator it = msh.begnode() ; it != msh.endnode() ; it++){
			pts[3*i] = (float)it->x;
			pts[3*i+1] = (float)it->y;
			pts[3*i+2] = 0;
			vars[1][i] = ( it->bvertex ?
						   (float) it->bvertex->reg_->ID :
						   -1 );
			i++;
		}
	}

	//allocate memory to mesh_element data
	conn = new int [ msh.nele() * 4 ] ;    //connectivity
	celltypes = new int [ msh.nele() ];    //cell type (vtk)
	vars[0] = new float [ msh.nele() ];    //region_id
	//fill mesh data
	i = 0; j = 0;
	for (list<eleblank*>::iterator it = msh.begele() ; it != msh.endele() ; it++){
		celltypes[i] = cell2visit ( (*it)->cellType() );
		vars[0][i] = (float) ( (*it)->regionID() ) ;
	    if (md.visualduplicate == 0 ){
			for (int k=0; k < (*it)->nNode() ; k++) {
				conn[j] = (*it)->dupl(k)->idx;
				j++;
			}			
		}
		else{
			idxele = &(*it)->idxGlob();
			for (int k=0; k < (*it)->nNode() ; k++) {
				conn[j] = (*idxele)(k);
				j++;
			}
		}
		i++;
	}

	
	//Saturation
	if (md.S.size() > 0){
		nvars+=2;
		//nodewise S
		if (md.visualduplicate == 0){
			vars[2] = new float[ msh.ndd() ] ;
			for (vector<Node>::iterator it = msh.begnode() ; it != msh.endnode() ; it++){
				for (list<DuplData>::iterator l = it->dd.begin() ; l != it->dd.end() ; l++)
					vars[2][l->idx] = (float)md.S.at(l->idx);
			}
		}
		else{
			vars[2] = new float[ msh.nnode() ] ;
			i = 0;
			for (vector<Node>::iterator it = msh.begnode() ; it != msh.endnode() ; it++){
				if (md.visualduplicate == 1)
					vars[2][i] = (float)fmin( md.S.at( it->dd.front().idx ) ,
											  md.S.at( it->dd.back().idx ) );
				else
					vars[2][i] = (float)fmax( md.S.at( it->dd.front().idx ) ,
											  md.S.at( it->dd.back().idx ) );				
				i++;
			}
		}
		//cellwise S
		vars[3] = new float[ msh.nele() ] ;
		i = 0;
		for (list<eleblank*>::iterator it = msh.begele() ; it != msh.endele() ; it++){
		    vars[3][i] = arma::mean( (*it)->lDatCnDis(md.S, 0) );
			i++;
		}
	}

	//Pressure
	if (md.P != NULL){
		nvars++;
	    if (md.visualduplicate == 0){
			vars[4] = new float[ msh.ndd() ] ;
			for (vector<Node>::iterator it = msh.begnode() ; it != msh.endnode() ; it++)
				for (list<DuplData>::iterator l = it->dd.begin() ; l != it->dd.end() ; l++)
					vars[4][l->idx] = (float)md.P[it->idx];
					
		}
		else{
			vars[4] = new float[ msh.nnode() ] ;
			for (i = 0 ; i < msh.nnode() ; i++ ) vars[4][i] = (float)md.P[i];
		}
	}
	
	
	//cellwise phase velocity
	if ((md.Lw.size() > 0) && (md.Ln.size() > 0)){
		nvars+=3;
		//allocate memmory for water, oil and total velocity
		vars[5] = new float[ msh.nele() * 3] ;
		vars[6] = new float[ msh.nele() * 3] ;
		vars[7] = new float[ msh.nele() * 3] ;
		//find v
		i = 0; 
		for (list<eleblank*>::iterator it = msh.begele() ; it != msh.endele() ; it++){
			parma = &(*it)->lDatCnCon(md.P, 0);
			pcarma = &(*it)->lDatCnDis(md.Pc, 1);
			lwarma = &(*it)->lDatCnDis(md.Lw, 2);
			lnarma = &(*it)->lDatCnDis(md.Ln, 3);
			kdarma = &(*it)->matKD();
			vw = -arma::mean(*lwarma) * (*kdarma) * (*parma) / md.dp ;
			vn = -arma::mean(*lnarma) * (*kdarma) * (  *parma  +  *pcarma  ) / md.dp ;
			vt = vw + vn ;
		    vars[5][3*i] = (float)vw(0); vars[5][3*i+1] = (float)vw(1); vars[5][3*i+2] = 0;
			vars[6][3*i] = (float)vn(0); vars[6][3*i+1] = (float)vn(1); vars[6][3*i+2] = 0;
			vars[7][3*i] = (float)vt(0); vars[7][3*i+1] = (float)vt(1); vars[7][3*i+2] = 0; 
			i++;
		}		
	}
	

	//write the mesh using VisIt writer
	write_unstructured_mesh(address, usebinary,
						    npts, pts,
							msh.nele(), celltypes, conn,
							nvars, vardim, centering, varnames, vars);

	//clear the data
	delete[] pts;
	delete[] celltypes;
	delete[] conn;
	for (i = 0 ; i < nvars ; i++ ) delete[] vars[i];

	FuncEnd();
}

/** @brief  writes a gmsh visualization file */
static void writevisual_gmsh(const char * address,MData& md,Mesh &msh){
	cell2gmsh(CellTri);
	ERRSET();
}

/** @brief  writes a tecplot visualization file */
static void writevisual_tecplot (const char * address,MData& md,Mesh &msh){
	cell2tecplot(CellTri);
	ERRSET();
}
/************************************************************************
 * computational stuff
 ************************************************************************/

/** @brief find a nodes slave saturations */
static void cmpnode_slave_s (MData &md, Node &node){
	FuncBegin();
	list<DuplData>::iterator i = node.dd.begin(); i++;
	for ( ; i != node.dd.end() ; i++){
		md.S.at( i->idx ) =
			md.J->sopp(md.S.at(node.dd.front().idx), node.dd.front().reg->pd,i->reg->pd);
	}
	FuncEnd();
}

/** @brief find a nodes slave delta saturation.
	
	not used (unstable probably because of round of errors ).
	instead use cmpnode_slave_s.
*/
static void cmpnode_slave_deltas (MData &md, Node &node){
	FuncBegin();
	list<DuplData>::iterator i = node.dd.begin(); i++;
	for ( ; i != node.dd.end() ; i++){
		md.dS.at( i->idx ) = md.dS.at( node.dd.front().idx ) *
			md.J->ds(md.S.at(node.dd.front().idx), node.dd.front().reg->pd,i->reg->pd);
	}
	FuncEnd();
}

/** @brief  find a nodes cappilary potential and mobilities*/
static void cmpnode_cappil_mobil (MData &md, Node &node){
	FuncBegin();
	for (list<DuplData>::iterator i = node.dd.begin(); i != node.dd.end() ; i++){
		md.Lw.at( i->idx ) = i->reg->kr->w( md.S.at (i->idx) );
	    md.Ln.at( i->idx ) = i->reg->kr->nw( md.S.at (i->idx) ) / md.dm;
		md.Pc.at(i->idx) =
			i->reg->pd*md.J->J(md.S.at(i->idx)) +  //cappilary term
			md.DgH.at(node.idx) ;                  //gravity
	}
	FuncEnd();
}

/** @brief find a nodes SIGMA( ds_i/ds_master * phi_i * v_i ) */
static void cmpnode_sphiv (MData &md, Node &node){
	FuncBegin();
	list<DuplData>::iterator i = node.dd.begin(); i++;
	md.SPhiV.at(node.idx) = md.VPhi.at(node.dd.front().idx);
	for ( ; i != node.dd.end() ; i++){
		md.SPhiV.at(node.idx) += md.VPhi.at(i->idx)
			* md.J->ds(md.S.at(node.dd.front().idx), node.dd.front().reg->pd,i->reg->pd);
	}
	FuncEnd();
}



/***************************************************************************
 * driver namespace
 **************************************************************************/

namespace driver{
	
	void initialize(int *argc,char **argv[], MData &md){
		FuncBegin();
		int i;
		PetscBool setf;
		string adrfullpetsc;
		
		//initialize md
		md.initialize();
		//get current directory
		for (i = 0 ; i < *argc ; i++){
			if ( strcmp( (*argv)[i] , swidir) == 0 ) break;
		}
		if ( (i == *argc) || (i+1 == *argc) ) {
			md.dir = "";
		}
		else{
			md.dir = (*argv)[i+1];
			for (i = md.dir.length() -1 ; i >= 0 ; i--){
				if (md.dir.at(i) == '/') md.dir.erase(i);
				else break;
			}			
			md.dir += '/';
		}
		//initialize petsc and read running mode
		adrfullpetsc = md.dir + adrpetsc;
		Error::code=PetscInitialize(argc,argv,adrfullpetsc.c_str(),help);ERRCHK();
		Error::code=PetscOptionsHasName(NULL,swisetf,&setf);ERRCHK();
		md.setfield = (setf == PETSC_TRUE ? 1 : 0);
		//report
		cout << "\nWorking directory found: " << md.dir << endl
			 << "Main data initialized successfuly" << endl
			 << "Petsc initialized successfuly" << endl
			 << "Running mode: " << (md.setfield ? "SetField" : "Solver" ) << endl;
	
		FuncEnd();
	}

	void readregion (MData &md, Mesh &msh){
		FuncBegin();
		
		AsciiFile cfile;
		string regtype;

		RegionBoundary *bndptr;
	    double bndval;
		string bndpstr, bndsstr;
		
		RegionPorousMat *matptr;
		string matk;
		
		RegionPorousFrac *fracptr;
		double frace, frack;

		double prpd, prphi;
		int prid;
		KFunc *prkr;

		cfile.open(md.dir+adrsol);
		cfile.efind("$regiondata");
		
		while ( cfile.find("$region") ){
			
			cfile(); cfile("type"); cfile(regtype,"type_value");
			//frac
			if (regtype.compare("frac") == 0){
				//read
				cfile(); cfile("id"); cfile(prid, "id_value");
				cfile(); cfile("phi"); cfile(prphi, "phi_value"); 
				cfile(); cfile("pd"); cfile(prpd, "pd_value"); 
				cfile(); cfile("kr"); prkr = readregion_kr(cfile); 
				cfile(); cfile("k"); cfile(frack, "k_value"); 
				cfile(); cfile("e"); cfile(frace, "e_value");
				//create
				fracptr = new RegionPorousFrac;
				fracptr->ID = prid;
				fracptr->phi = prphi;
				fracptr->pd = prpd;
				fracptr->kr = prkr;
				fracptr->k = frack;
				fracptr->e = frace;
				msh.addRegion(fracptr);
			}
			//mat
			else if (regtype.compare( "mat" ) == 0){
				//read
				cfile(); cfile("id"); cfile(prid, "id_value");
				cfile(); cfile("phi"); cfile(prphi, "phi_value");
				cfile(); cfile("pd"); cfile(prpd, "pd_value"); 
				cfile(); cfile("kr"); prkr = readregion_kr(cfile); 
				cfile(); cfile("k"); cfile(matk, "k_value"); 
				//create
				matptr = new RegionPorousMat(matk);
				matptr->ID = prid;
				matptr->phi = prphi;
				matptr->pd = prpd;
				matptr->kr = prkr;
				msh.addRegion(matptr);
			}
			//bnd
			else if (regtype.compare( "bnd" ) == 0 ){
				//read
				cfile(); cfile("id"); cfile(prid, "id_value");
				cfile(); cfile("stype"); cfile(bndsstr, "stype_value");
				cfile(); cfile("ptype"); cfile(bndpstr, "ptype_value");
				cfile(); cfile("value"); cfile(bndval, "value_value");
				//create
				bndptr = new RegionBoundary;
				bndptr->ID = prid;
				bndptr->val[0] = bndval;
				//ptype
				if ( bndpstr.compare( "pconst") == 0 )
					bndptr->ptype = RegionBoundary::PPConst;
				else if ( bndpstr.compare( "qconst") == 0 )
					bndptr->ptype = RegionBoundary::PQConst;
				else {
					Error::mess <<"ptype " << bndpstr  << " not valid. "
								<< cfile.fn << ", line: " <<cfile.ln;
					ERRSET();
				}
				//stype
				if ( bndsstr.compare( "sconst")  == 0)
					bndptr->stype = RegionBoundary::SSConst;
				else if ( bndsstr.compare( "gpczero") == 0)
					bndptr->stype = RegionBoundary::SGPZero;
				else {
					Error::mess <<"stype " << bndsstr  << " not valid. "
								<< cfile.fn << ", line: " <<cfile.ln;
					ERRSET();
				}
				//add
				msh.addRegion(bndptr);
			}
			
			else{
				Error::mess << "region type: " << regtype << " not found in "
							<< cfile.fn << ", line: " <<cfile.ln;
				ERRSET();
			}
		}

		//report
		cout << "\nRegions read successfuly: " << endl;
		for (list<Region*>::iterator i = msh.begreg() ; i != msh.endreg(); i++)
			cout << (*i)->name() << endl;
		
		FuncEnd();
	}
	
	void readfixed (MData &md){
		FuncBegin();

	    AsciiFile fl;
		string tstr;

		fl.open(md.dir + adrsol);
		fl.efind("$fixeddata");

		//dim-less numbers and gravity
		fl(); fl("dm"); fl(md.dm,"dm_value");
		fl(); fl("dn"); fl(md.dn,"dn_value");
		fl(); fl("dp"); fl(md.dp,"dp_value");
		fl(); fl("dgw"); fl(md.grav.gw, "dimless wetting phase gravity");
		fl(); fl("dgn"); fl(md.grav.gn, "dimless non-wetting phase gravity");
		fl(); fl("gorigin"); fl(md.grav.x0, "zero-g x"); fl(md.grav.y0, "zero-g y");
		fl(); fl("gradz");  fl(md.grav.xdir, "g-x-direc"); fl(md.grav.ydir, "g-y-direc");
 
		//time
		fl(); fl("starttime"); fl(md.t,"t_value");
		fl(); fl("stoptime"); fl(md.tEnd,"tend_value");
		fl(); fl("timestep"); fl(md.dt,"dt_value");
		
		//time-stepping
		fl(); fl("maxtimestep"); fl(md.dtM,"dtM_value");
		fl(); fl("mintimestep"); fl(md.dtm,"dtm_value");
		fl(); fl("maxdeltas"); fl(md.dsM,"dsM_value");
		fl(); fl("mindeltas"); fl(md.dsm,"ds_value");
		fl(); fl("beta"); fl(md.beta,"beta_value");
		fl(); fl("maxtimeiteration"); fl(md.dnItM,"dnItM_value");

		//file writing options
		fl(); fl("nfilebegin"); fl(md.nFile,"nFile_value");
		fl(); fl("nfileend"); fl(md.NFile,"NFile_value");

		//jmodel
		fl(); fl("jmodel"); fl(tstr,"jmodel_value");
		if ( tstr.compare("firooz") == 0 )	md.J = new JFuncFirooz;
		else if ( tstr.compare("zero") == 0 ) md.J = new JFuncZero;
		else if ( tstr.compare("linear") == 0 ) md.J = new JFuncLinear;
		else if ( tstr.compare("vang") == 0 ){
			double m,e;
			fl(m,"mvalue");
			fl(e,"evalue");
			md.J = new JFuncVang(m,e);
		}
		else if ( tstr.compare("brooks") == 0 ){
			double lambda;
			fl(lambda,"lambda_value");
			md.J = new JFuncBrooks(lambda);
		}
		else{
			Error::mess << "jmodel_" << tstr << " not supported. "
						<< fl.fn << " line " << fl.ln ;
			ERRSET();		
		}

		//meshtype
		fl(); fl("meshtype"); fl(tstr,"meshtype_value");
		if ( tstr.compare("gmsh") == 0 )	md.meshtype = MData::MeshGmsh;
		else if ( tstr.compare("triangle") == 0 ) md.meshtype = MData::MeshTriangle; 
		else {
			Error::mess << "mesh_" << tstr << "not supported. "
						<< fl.fn << " line " << fl.ln ;
			ERRSET();
		}

		//visualtype
		fl(); fl("visualtype"); fl(tstr,"visualtype_value");
		if ( tstr.compare("vtk") == 0 )	md.visualtype = MData::VisualVtk;
		else if ( tstr.compare("tecplot") == 0 ){
			Error::mess << "tecplot not supported yet. "
						<< fl.fn << " line " << fl.ln ;
			ERRSET();
		}
		else if ( tstr.compare("gmsh") == 0 ) {
			Error::mess << "gmsh not supported yet. "
						<< fl.fn << " line " << fl.ln ;
			ERRSET();
		}
		else {
			Error::mess << "mesh_" << tstr << "not supported. "
						<< fl.fn << " line " << fl.ln ;
			ERRSET();
		}
		fl(); fl("visualduplicate"); fl(md.visualduplicate, "visualduplicate_value");

		//set values which were not read
		md.t0 = md.t;
		md.nFile0 = md.nFile;
		md.tWrite = (md.tEnd - md.t0) / (md.NFile - 1);

		//report
		//print stuff
		cout << "\nFixed data was read successfuly " << endl
			 << "DM: " << md.dm << endl
			 << "DN: " << md.dn << endl
			 << "DP: " << md.dp << endl
			 << "DGW: " << md.grav.gw << endl
			 << "DGN: " << md.grav.gn << endl
			 << "GOrigin: " << md.grav.x0 << "\t" << md.grav.y0 << endl
			 << "GDirection: " << md.grav.xdir << "\t" << md.grav.ydir << endl			
			 << "StartTime: " << md.t << endl
			 << "StopTime: " << md.tEnd << endl
			 << "TimeStep: " << md.dt << endl
			 << "MaxTimeStep: " << md.dtM << endl
			 << "MinTimeStep: " << md.dtm << endl
			 << "MaxDeltaS: " << md.dsM << endl
			 << "MinDeltaS: " << md.dsm << endl
			 << "Beta: " << md.beta << endl
			 << "MaxTimeIteration: " << md.dnItM << endl
			 << "NFileBegin: " << md.nFile << endl
			 << "NFile: " << md.NFile << endl
			 << "JModel: " << md.J->name() << endl
			 << "MeshType: " << md.meshtype << endl
			 << "VisualType: " << md.visualtype << endl
			 << "TWrite: " << md.tWrite << endl;
				
		FuncEnd();
	}
	
	void readmesh(MData &md, Mesh &msh){
		FuncBegin();

		//read the file
		switch(md.meshtype){
		case MData::MeshTriangle:
			readmesh_triangle(md,msh);
			break;
		case MData::MeshGmsh:
			readmesh_gmsh(md,msh);
			break;
		default:
			Error::mess << "invalid mesh" ;
			ERRSET();
		}
		//construct the mesh
		msh.constructGeoParams(md.J->cmp, md.dp, md.A, md.grav);
		//report
		cout << "\nMesh file(s) was read successfuly." << endl;
		
		FuncEnd();
	}
	
	void readinitial(MData &md, Mesh &msh){
		FuncBegin();
		
		AsciiFile fl;
	    char mod;
		double s0;
		bool uni;

		//allocate memory
		md.S.resize(msh.ndd(),0);
		//read file
		fl.open(md.dir+adrini);
		fl.efind("%initialcondition");
		fl(); fl("%mod"); fl(mod,"mod_value");	
		switch(mod){
		case 'u':
			fl(s0,"uniform s0 value");
			uni = true;
			break;
		case'n':
			uni = false;
			break;		
		default:
			Error::mess << "In file " << fl.fn << "line " << fl.ln
						<< " mod \"" << mod << "\" is invalid";
			ERRSET();
		}	
		for (vector<Node>::iterator i = msh.begnode() ; i != msh.endnode() ; i++){
			if (uni){
				md.S.at(i->dd.front().idx) = s0;
			}
			else{
				fl(); fl( md.S.at(i->dd.front().idx) ,"non-uni s value");
			}
		}
		
		//print
		std::cout << "Initial condition read successfuly from " << fl.fn << std::endl;
		FuncEnd();
	}
	
	void setfield(MData &md, Mesh &msh){
		FuncBegin();
		
		const int n_max = 10;
		Poly *ply[n_max];
		double val[n_max],x0,y0,x1,y1,r;
		int n,m;
		string cmd;
		AsciiFile afile;
		
		//open the set field file
		afile.open(md.dir + adrsol);
		afile.efind("$setfield");
		afile(); afile(n,"Number of commands");
		if ( n > n_max ){
			Error::mess << "maximum number of commands in " << afile.fn << " is " << n_max;
			ERRSET();
		}	
		cout << "Setting field according to " << afile.fn << endl;
	
		//read polygons
		for (int i = 0 ; i < n ; i++){
			afile(); afile(cmd,"Polygon type");
			if (cmd.compare("circle") == 0){
				afile(x0,"x0");	afile(y0,"y0"); afile(r,"r");
				ply[i] = new Circle(x0,y0,r);
			}
			else if (cmd.compare("rectangle") == 0){
				afile(x0,"x0");	afile(y0,"y0");afile(x1,"x1");afile(y1,"y1");
				ply[i] = new Rectangle(x0,y0,x1,y1);
			}
			else{
				Error::mess << afile.fn << " line " << afile.ln << " poly " << cmd << " is invalid";
				ERRSET();
			}
			afile(val[i],"S value");
			cout << ply[i] ->name() << " with value = " << val[i] << endl;
		}
		afile.close();
		
		//change S according to given commands
		for(vector<Node>::iterator i = msh.begnode() ; i < msh.endnode() ; i++){
			m = n - 1;
			for (; m >= 0 ;m--)	if (ply[m]->isin(i->x, i->y)) break;
			if (m != -1){
				md.S.at(i->dd.front().idx) = val[m] ;
			}
		    cmpnode_slave_s(md, *i);
		}
		
		//destroy the poly's
		for (int i = 0 ; i < n ; i++) delete ply[i];
		
		//write the initialcondition
		writeinitial(md,msh,md.dir+adrini,false,false );
		//write a vtk file
		writevisual((md.dir+"newfield.vtk").c_str(),md,msh);
		
		//print
		cout << "Field was set and saved successfully in " << md.dir+adrini << endl;	
		
		FuncEnd();
		
	}

	void preparedata (MData &md, Mesh &msh){
		FuncBegin();

		const arma::rowvec *vol;
	    md.Pc.resize(msh.ndd(),0);
		md.SPhiV.resize(msh.nnode(),0);
		md.VPhi.resize(msh.ndd(),0);
		md.Lw.resize(msh.ndd(),0);
		md.Ln.resize(msh.ndd(),0);
		md.Fs.resize(msh.nnode(),0);
		md.dS.resize(msh.ndd(),0);
		md.DgH.resize(msh.nnode(),0); //gravity
		//Vectors
		Error::code=VecCreateSeq(PETSC_COMM_SELF, msh.nnode(), &md.Pvec);ERRCHK();	
		Error::code=VecDuplicate(md.Pvec, &md.b);ERRCHK();
		Error::code=VecGetArrayRead(md.Pvec, &md.P);ERRCHK();
		Error::code=VecSet(md.Pvec, 0);ERRCHK();                //gravity - starts with zero, might cause problems.
		//KSP
		Error::code=KSPCreate(PETSC_COMM_SELF, &md.ksp);ERRCHK();
		Error::code=KSPSetOperators(md.ksp, md.A, md.A);ERRCHK();
		Error::code=KSPSetFromOptions(md.ksp);ERRCHK();
		//SphiV
		for (list<eleblank*>::iterator i = msh.begele() ; i != msh.endele() ; i++ ){
			//calc SphiV
			for (int j = 0 ; j < (*i)->nNode() ; j++ ){
				vol = &(*i)->matVolume();
				md.VPhi.at( (*i)->dupl(j)->idx ) += (*i)->dupl(j)->reg->phi * (*vol)(j);
			}
		}
		//initial data
		for(vector<Node>::iterator i = msh.begnode(); i!= msh.endnode() ; i++){
			md.DgH.at(i->idx) = (md.grav.gn - md.grav.gw) * fml::findHeight(i->x, i->y, md.grav); //gravity
			cmpnode_slave_s(md, *i);
			cmpnode_cappil_mobil(md, *i);
			cmpnode_sphiv(md, *i);
		}
		//upwind stuff
		for (list<eleblank*>::iterator i = msh.begele() ; i != msh.endele() ; i++ ){
			(*i)->fndUpW(md.P);
			(*i)->fndUpN(md.P,md.Pc);
		}
		//report
		std::cout << "External data initialized successfuly ..." << std::endl;
		FuncEnd();
	}
	
	void marchintime(MData &md, Mesh &msh){
		FuncBegin();

		//set initial values
		bool flag = false;
		double ds = 0 ;
		int it ; double res;
		const arma::vec *rhs;
	    const arma::mat *lhs;
		const arma::ivec *idx;
		md.dcT = clock();
		md.dnIt = 0;

		//assemble P equation
		Error::code=VecSet(md.b, 0);ERRCHK();
		Error::code=MatZeroEntries(md.A);ERRCHK();
		for (list<eleblank*>::iterator i = msh.begele() ; i != msh.endele() ; i++){
			//update non wetting upwind
			(*i)->fndUpN(md.P, md.Pc);
			//get mats
			lhs = &(*i)->lhsP(md.Lw,md.Ln);
			rhs = &(*i)->rhsP(md.Ln,md.Pc, 0);
			idx = &(*i)->idxGlob();	
			Error::code=MatSetValues(md.A,
									 (*i)->nNode(), idx->memptr(),
									 (*i)->nNode(), idx->memptr(),
									 lhs->memptr(), ADD_VALUES);ERRCHK();
			Error::code=VecSetValues(md.b,
									 (*i)->nNode(), idx->memptr(),
									 rhs->memptr(),ADD_VALUES);ERRCHK();
		}
		Error::code=MatAssemblyBegin(md.A, MAT_FINAL_ASSEMBLY);ERRCHK();
		Error::code=MatAssemblyEnd(md.A, MAT_FINAL_ASSEMBLY);ERRCHK();
		Error::code=VecAssemblyBegin(md.b);ERRCHK();
		Error::code=VecAssemblyEnd(md.b);ERRCHK();

		//force boundary condition
		for (list<BVertexCQ*>::iterator i = msh.begbvertex() ; i != msh.endbvertex() ; i++)
			(*i)->assemP(md.A, md.b);
		
		//Solve the p equation
		Error::code=KSPSolve(md.ksp, md.b, md.Pvec);ERRCHK();
		
		//assemble S equation
		//make the flux zero
		for (int i = 0 ; i < msh.nnode() ; i++) md.Fs.at(i) = 0;
		for (list<eleblank*>::iterator i = msh.begele() ; i != msh.endele() ; i++){
			//update wetting upwind node
			(*i)->fndUpW(md.P);
			//get mat and idx
			rhs = &(*i)->rhsS(md.Lw, md.P, 0);
			idx = &(*i)->idxGlob();
            //assemble
			for (int j = 0 ; j < (*i)->nNode();  j++)
				md.Fs.at( (*idx)(j) ) += (*rhs)(j) ;
		}
		
		//force boundary condition
		for (list<BVertexCQ*>::iterator i = msh.begbvertex() ; i != msh.endbvertex() ; i++){
			(*i)->findQAll(md.Fs, md.P, md.Lw, md.Ln);
			(*i)->assemS(md.Fs);
		}
		
		//solve for ds
		do {
			md.dnIt++;
			ds = 0;
			for (vector<Node>::iterator j = msh.begnode() ; j < msh.endnode() ; j++){
				md.dS.at(j->dd.front().idx) =	md.Fs.at(j->idx) * md.dt /md.dn / md.SPhiV.at( j->idx ) ;
				ds = fmax( ds , fabs(md.dS.at(j->dd.front().idx)) );
			}			
			if (ds > md.dsM){
				md.dt /= md.beta;
				if ( (md.dt < md.dtm) || ( md.dnIt > md.dnItM ) ){
				Error::mess << "either min_dt or max_s_iter error." << endl
							<< "dt: " << md.dt << "\tdtm: " << md.dtm
							<< " dn_it: " << md.dnIt << "\tdn_max: " << md.dnItM;
				ERRSET();
				}
			}
			else {
				if (ds < md.dsm) flag = true;
				break;
			}
		}while(true);
		
		//update everything
		for (vector<Node>::iterator i = msh.begnode() ; i < msh.endnode() ; i++){
			md.S.at(i->dd.front().idx) += md.dS.at(i->dd.front().idx);
			cmpnode_slave_s(md, *i);
		    cmpnode_sphiv(md, *i);
			cmpnode_cappil_mobil(md, *i);
		}

	//update the fluxes
	for (list<BVertexCQ*>::iterator i = msh.begbvertex() ; i != msh.endbvertex() ; i++){
		md.qIn   += fmax( 0 , (*i)->mGamma()  ) * md.dt / md.dp;
		md.qOut  -= fmin( 0 , (*i)->mGamma()  ) * md.dt / md.dp;
		md.qWin  += fmax( 0 , (*i)->mGammaW() ) * md.dt / md.dp;
		md.qWout -= fmin( 0 , (*i)->mGammaW() ) * md.dt / md.dp;
	}
	
	//update time step
	md.t += md.dt;
    md.nIt += md.dnIt;
	md.dcT = (clock()-md.dcT) / CLOCKS_PER_SEC;
	md.cT += md.dcT ;
	Error::code = KSPGetResidualNorm(md.ksp,&res);ERRCHK();
	Error::code = KSPGetIterationNumber(md.ksp,&it);ERRCHK();
	if (flag ) md.dt = fmin ( md.dtM, md.dt * md.beta );

	//report
	cout << left ;
	cout << setw(10) << "t_COMP: " << setw(15) << md.t 
						<< setw(10)<< "t_clock: " << setw(15)<< md.cT 
						<< setw(10)<< "n_it: " << setw(15)<< md.nIt <<endl
						<< setw(10)<< "dt_comp: " << setw(15)<< md.dt
						<< setw(10)<< "dt_clock: " << setw(15)<< md.dcT
						<< setw(10)<< "dn_it: " << setw(15)<< md.dnIt <<endl
						<< setw(10)<< "ds_max: "<< setw(15) << ds
						<< setw(10)<< "ksp it: "<< setw(15) << it 
						<< setw(10)<< "ksp res: " << setw(15)<< res << "\n\n";

		FuncEnd();
	}
	
	void writeintime(MData &md, Mesh &msh, bool force){
    	FuncBegin();

		if ( force || ( (md.t - md.t0) > (md.nFile - md.nFile0) * md.tWrite ) ){
			fstream fl;
			stringstream ss;
			double vw;
			fl << left;
			
			//open result.flow file
			ss << md.dir+adrresult << ".flow";
			fl.open(ss.str().c_str(),fstream::app | fstream::out);
			if (!fl.is_open()){
				Error::mess << ss.str() << " could not be openned.";
				ERRSET();
			}
			
			//write header
			if(fl.tellp() == 0)
				fl << setw(5) << "# n"
				   << setw(15) << "t"
				   << setw(15) << "Q_in"
				   << setw(15) << "Q_out"
				   << setw(15) << "Q_w_in"
				   << setw(15) << "Q_w_out"
				   << setw(15) << "V_w" << endl ;
			if (md.nFile == 0) fl << endl;
			
			//find V_w
			vw = 0;
			for(vector<Node>::iterator i = msh.begnode() ; i < msh.endnode() ; i++)
				for(std::list<DuplData>::iterator j = i->dd.begin() ; j != i->dd.end() ; j++)
					vw += md.VPhi.at(j->idx) * md.S.at(j->idx);
			
			//write data
			fl << setw(5) << md.nFile
			   << setw(15) << md.t
			   << setw(15) << md.qIn
			   << setw(15) << md.qOut
			   << setw(15) << md.qWin
			   << setw(15) << md.qWout
			   << setw(15) << vw << endl;
			
			//close the file
			fl.close();
			
			//write restart file
			ss.str("");
			ss << md.dir+adrrestart << "." << md.nFile;
			writeinitial(md, msh, ss.str(), true, false);
			//write visual file
			ss.str("");
			ss << md.dir+adrresult << "." << md.nFile;
			writevisual(ss.str().c_str(), md, msh);
			
			//increase file number
			md.nFile++;
		}
		
		FuncEnd();
	}
	
	void writeinitial(MData &md, Mesh &msh,const string& adr,
					  const bool restart, const bool octave){
		FuncBegin();
		
		fstream ofile;
		ofile.open(adr.c_str(), fstream::out);
		if (!ofile.is_open()) {
			Error::mess << "could not open " << adr ;
			ERRSET();
		}
		ofile << setprecision(12);
		
		if (restart) ofile << "#Restart file @ time = " << md.t << endl;
		ofile << "%initialcondition" << endl;
		ofile << "%mod n" << endl;
		
		for(vector<Node>::iterator i = msh.begnode() ; i < msh.endnode() ; i++){
			ofile << md.S.at(i->dd.front().idx) << " ";
			if (octave){
				ofile << md.P[i->idx] << " "
					  << md.Pc.at(i->dd.front().idx) << " "
					  << i->x << " ";
			}
			ofile << endl;
		}
 
		ofile.close();
		
		FuncEnd();
	}
	
	void writevisual(const char * address,MData &md,Mesh &msh){
		FuncBegin();
		switch ( md.visualtype ){
		case MData::VisualVtk:
			writevisual_vtk(address, md, msh);
			break;
		case MData::VisualGmsh:
			writevisual_gmsh(address, md, msh);
			break;
		case MData::VisualTecplot:
			writevisual_tecplot(address, md, msh);
			break;
		}

		FuncEnd();
	}	

	void writetest (MData &md, Mesh &msh, const uint num){
		FuncBegin();

		fstream fl;
		stringstream fladdress;
		const vector<double> *ptrpc,*ptrlw,*ptrln,*ptrs;
		const double *ptrp;
		PetscViewer view;

		//set the ptrs
		ptrp = (md.P ? md.P : NULL);
		ptrpc = (md.Pc.size() > 0 ? &md.Pc : NULL);
		ptrlw = (md.Lw.size() > 0 ? &md.Lw : NULL);
		ptrln = (md.Ln.size() > 0 ? &md.Ln : NULL);
		ptrs = (md.S.size() > 0 ? &md.S : NULL);
		
		//open the file
		fladdress << md.dir + "mydata." << num <<".log";
		fl.open(fladdress.str().c_str(),fstream::out);

		//write the nodes
		fl << "XXXXXX_________Nodes__________XXXXXX" << endl;
		for (vector<Node>::iterator i = msh.begnode() ; i < msh.endnode() ; i++){
			fl << i->idx << " "  << i->x << " "  << i->y << " "  << i->n_dd << " ";
			for (list<DuplData>::iterator j = i->dd.begin(); j != i->dd.end() ; j++)
				fl << " | " <<j->idx << " " << j->reg->ID << " ";
			fl << endl;
		}
		fl << endl;

		//write the boundary
		fl << "XXXXXX_________BOUNDARY___VERTICES__________XXXXXX\n" ;
		for (list<BVertexCQ*>::iterator i= msh.begbvertex() ; i != msh.endbvertex() ; i++)
			fl << (*i)->name() << endl;
		fl << endl;

		//write the elements
	    fl << "XXXXXX________________ElEMENTS______________XXXXXX\n" ;
		for (list<eleblank*>::iterator i= msh.begele() ; i != msh.endele() ; i++){
			fl <<  (*i)->name(ptrs, ptrp,ptrpc, ptrlw, ptrln) ;
			if (ptrlw && ptrln) fl << "lhsP:\n" << (*i)->lhsP(*ptrlw,*ptrln);
			if (ptrln && ptrpc) fl << "rhsP:\n" << (*i)->rhsP(*ptrln, *ptrpc, 0);
			if (ptrlw && ptrp) fl <<  "rhsS:\n" << (*i)->rhsS(*ptrlw, ptrp, 0);
			fl << endl;
		}
		fl << endl << fixed << left;
		
		//print the external data with size n
		fl << setw(15) << "P"
		   << setw(15) <<"SPhiV"
		   << setw(15) << "Fs" <<  endl;
		for (int i = 0 ; i < (int)md.Fs.size() ; i++)
			fl << setw(15) << md.P[i] 
			   << setw(15) << md.SPhiV.at(i)
			   << setw(15) << md.Fs.at(i) <<  endl;
		fl << endl;
		
		//write data with length nd.nn
		fl << setw(15) << "S"			
		   << setw(15) <<"Pc"
		   << setw(15) << "Lw"
		   << setw(15) << "Ln"
		   << setw(15) << "VPhi"
		   << setw(15) << "dS"
		   << endl;
		if(md.dS.size() > 0)
			for (int i = 0 ; i < msh.ndd() ; i++)
				fl << setw(15) << md.S.at(i)
				   << setw(15) << md.Pc.at(i) 
				   << setw(15) << md.Lw.at(i)
				   << setw(15) << md.Ln.at(i) 
				   << setw(15) << md.VPhi.at(i) 
				   << setw(15) << md.dS.at(i) << endl;
		fl << endl;

		//close the file
		fl.close();

		//write lhs and rhs
		fladdress.str("");
		fladdress << md.dir + "petsc." << num <<".log";			;
		Error::code=PetscViewerASCIIOpen(PETSC_COMM_SELF,fladdress.str().c_str(),&view);ERRCHK();
		if (md.A) Error::code=MatView(md.A, view);ERRCHK();
		if (md.b) Error::code=VecView(md.b, view);ERRCHK();		
		Error::code=PetscViewerDestroy(&view);ERRCHK();
		
		FuncEnd();
	}
}
