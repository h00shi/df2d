/** @file bvertex.cpp.
	cpp file for bvertex.hpp
*/

#include "bvertex.hpp"
#include "formula.hpp"
#include "error.hpp"
#include <sstream>

/***************************************************************************
 * Constant Flux (base)
 **************************************************************************/
BVertexCQ* BVertexCQ::neww(Node *nd, RegionBoundary *reg,
						  RegionBoundary::PType ptype){
	FuncBegin();
	
	switch(ptype){
	case RegionBoundary::PPConst:
		return new BVertexCP(nd,reg);
	case RegionBoundary::PQConst:
		return new BVertexCQ(nd,reg);
	default:
		Error::mess << "ptype not recognized: " << ptype;
		ERRSET();
		return (BVertexCQ*) NULL;
	}
	
	FuncEnd();
}

void BVertexCQ::findQW(const std::vector<double> &F,
					   const std::vector<double> &Lw, const std::vector<double> &Ln){
	FuncBegin();
	
	if (reg_->stype == RegionBoundary::SSConst) mGammaW_ = -F.at(self_->idx);
	else if(reg_->stype == RegionBoundary::SGPZero)
		mGammaW_ = ( mGamma_  - Ln.at(self_->dd.front().idx) * kdgdzna) *
			Lw.at(self_->dd.front().idx) /
			( Lw.at(self_->dd.front().idx) + Ln.at(self_->dd.front().idx) );

	FuncEnd();
}

void BVertexCQ::constructL(){
	FuncBegin();

	length_ = 0;
	if (pre_) length_ += fml::lineLength(pre_->x, pre_->y, self_->x, self_->y) / 2.;
	if (next_) length_ += fml::lineLength(self_->x, self_->y, next_->x, next_->y) / 2.;
	if( (!pre_) && (!next_) ) length_ = 1;

	FuncEnd();
}

void BVertexCQ::constructNA(const Gravity &grav){
	FuncBegin();
	arma::vec2 tmp;
	RegionPorousMat *matptr;
	
	nA_.zeros();
	if (pre_) {
		fml::lineNA(pre_->x, pre_->y, self_->x, self_->y, tmp);
		nA_ += tmp / 2.;
	}
	if (next_){
		fml::lineNA(self_->x, self_->y, next_->x, next_->y, tmp);
		nA_ += tmp / 2.;
	}
	if( (!pre_) && (!next_) ){ //normal vector not supported in 1d or source terms
		nA_.zeros(); 
		std::cout << " BoundaryWarning: node " << self_->idx
			 << " is source-sink term , setting nA to zero\n";
	}

	//compute gravity
	switch( self_->dd.front().reg->dim() ){
	case 1:
		std::cout << " BoundaryWarning: node " << self_->idx
			 << " is boundary but has a fracture as master region"
			 << ", setting  K*(gn-gw)*Grad(z)*n*A = 0\n";
		kdgdzna = 0;
		break;
	case 2:
		matptr = (RegionPorousMat*)self_->dd.front().reg;
		tmp = matptr->k * (grav.gn - grav.gw) /	sqrt( pow(grav.xdir,2) + pow(grav.ydir,2)) * nA_;
		kdgdzna = tmp(0)*grav.xdir + tmp(1)*grav.ydir;
		break;
	}

	FuncEnd();
}

void BVertexCQ::constructBase(Node *nd, RegionBoundary *reg){
	FuncBegin();

	mGamma_ = mGammaW_ = 0;
	reg_ = reg;
	pre_ = next_ = (Node*) NULL;
	self_ = nd;
	length_ = 0;
	nA_.zeros();

	if (self_->bvertex) {
		Error::mess << "node " << nd->idx << " was considered as boundary twice."
					<< " once to " << self_->bvertex->reg_->ID
					<< " second time to " << reg_->ID ;
		ERRSET();
	}
	else self_->bvertex = this;

	FuncEnd();
}

std::string BVertexCQ::nameB() const{
	FuncBegin();

	std::stringstream ss;

	ss << "l: " << length_
		 << " nA: " << nA_(0) << ", " << nA_(1)
		 << " mGamma: " << mGamma_ << " mGammaW: " << mGammaW_
		 << " pre: " << (pre_ ? pre_->idx : -1)
		 << " self: " << (self_ ? self_->idx : -1)
		 << " next: " << (next_ ? next_->idx : -1)
		 << " reg: " << reg_->ID ;
	return ss.str();

	FuncEnd();
}

void BVertexCQ::findQAll(const std::vector<double> &F, double const * const P,
						 const std::vector<double> &Lw , const std::vector<double> &Ln){
	FuncBegin();
	findQW(F,Lw,Ln);
	FuncEnd();
}

void BVertexCQ::assemS(std::vector<double> &F) const{
	FuncBegin();
	F.at(self_->idx) += mGammaW_;
	FuncEnd();
}

void BVertexCQ::assemP(Mat A, Vec b) {
	FuncBegin();
	double rhs = -mGamma_;
	Error::code=VecSetValue(b , self_->idx, rhs, ADD_VALUES);ERRCHK();
	FuncEnd();
}

int BVertexCQ::addNeigh(Node *nd){
	FuncBegin();

	if (nd->bvertex){
		if (nd->bvertex->reg_ != reg_){
			Error::mess << "node " << nd->idx << " wanted to be added to node "
						<< self_->idx << " as a boundary neighbour. However region mismatch!"
						<< " ID1: " << nd->bvertex->reg_->ID << " ID2: " << reg_->ID;
			ERRSET();
		}
	}
	if(!pre_) {
		pre_ = nd ;
		return 1;
	}
	else if(!next_){
		next_ = nd;
		return 2;
	}
	else return 0;

	FuncEnd();
}

int BVertexCQ::checkNext(Node *nd){
	FuncBegin();

	Node *tmp;
	if (next_ == nd) return 1;
	else if (pre_ == nd){
		tmp = pre_;
		pre_ = next_;
		next_ = tmp;
		return 2;
	}
	else return 0;
															 
	FuncEnd();
}

int BVertexCQ::checkPre(Node *nd){
	FuncBegin();
	
	Node *tmp;
	
	if (pre_ == nd) return 1;
	else if (next_ == nd){
		tmp = pre_;
		pre_ = next_;
		next_ = tmp;
		return 2;
	}
	else return 0;

	FuncEnd();
}
	
void BVertexCQ::constructGeoParams(const double dp, Mat A, const Gravity &grav){
	constructL();
	constructNA(grav);
	mGamma_=length_ * reg_->val[0] * dp;
}

std::string BVertexCQ::name() const{
	FuncBegin();
	return "BVertexCQ " + nameB();
	FuncEnd();
}

BVertexCQ::BVertexCQ(Node *nd, RegionBoundary *reg){
	FuncBegin();
	if (reg->ptype != RegionBoundary::PQConst){
		Error::mess << "Region must be PQConst";
		ERRSET();
	}
	constructBase(nd, reg);
	FuncEnd();
}

/***************************************************************************
 * Constant Pressure 
 **************************************************************************/

void BVertexCP::findQAll(const std::vector<double> &F, double const * const P,
						 const std::vector<double> &Lw , const std::vector<double> &Ln){
	FuncBegin();

	mGamma_ = rhs_;
	for (int i = 0 ; i < nConn_ ; i++) mGamma_ -= lhs_(i) * P[ conn_(i) ] ;
	findQW(F,Lw,Ln);

	FuncEnd();
}

void BVertexCP::assemP(Mat A, Vec b) {
	FuncBegin();

	double const *vals;

	Error::code=MatGetRow(A,self_->idx,NULL,NULL,&vals);ERRCHK();
	for (int i = 0 ; i < nConn_ ; i++) lhs_(i) = vals[i] ;
	Error::code=MatRestoreRow(A,self_->idx,NULL,NULL,&vals);ERRCHK();
	Error::code=VecGetValues(b, 1, &self_->idx, &rhs_);ERRCHK();

	Error::code=MatZeroRows(A, 1, &self_->idx, 1, NULL, NULL);ERRCHK();
	Error::code=VecSetValue(b , self_->idx, p_, INSERT_VALUES);ERRCHK();
	
	FuncEnd();
}

void BVertexCP::constructGeoParams(const double dp, Mat A,const Gravity &grav){
	FuncBegin();

	int const *conn;
	int nconn;

	constructL();
	constructNA(grav);
	
	Error::code=MatGetRow(A,self_->idx,&nconn,&conn,NULL);ERRCHK();
	nConn_ = nconn;
	lhs_.resize(nConn_);
	lhs_.zeros();
	conn_.resize(nConn_);
	for (int i = 0 ; i < nConn_ ; i++) conn_(i) = conn[i];
	Error::code=MatRestoreRow(A,self_->idx,&nconn,&conn,NULL);ERRCHK();

	FuncEnd();
}

std::string BVertexCP::name() const{
	FuncBegin();

	std::stringstream ss;
	ss << "BVertexCP " << nameB() << "\n"
		 << "p: " << p_ << " nConn: " << nConn_ << " rhs: " << rhs_ << "\n"
		 << "conn: " << conn_ << "lhs: " << lhs_ ;
	return ss.str();

	FuncEnd();
}

BVertexCP::BVertexCP(Node *nd, RegionBoundary *reg){
	FuncBegin();

	if (reg->ptype != RegionBoundary::PPConst){
		Error::mess << "Region must be PPConst";
		ERRSET();
	}
	p_ = reg->val[0];
	rhs_ = 0;
	nConn_ = 0;
	constructBase(nd, reg);

	FuncEnd();
}
