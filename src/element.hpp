/** @file element.hpp
 * Header file for element classes - with only defenitions.
 * Written by Shayan Hoshyari
 * As a part of DF_2d
 * @ingroup mesh_module
 */

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "element_bone.hpp"
#include "error.hpp"
#include <sstream>

/****************************************************************************
 * Element - first level
 ***************************************************************************/
inline Element* Element::neww(RegionPorous *reg, Node *nd[], const CellType celltype){
	FuncBegin();
	
	switch (celltype){
	case CellQuad:
		return new elequad(reg, nd ) ;
		break;
	case CellTri:
		return new eletri(reg, nd ) ;
		break;
	case CellLine:
		return new elefrac(reg, nd ) ;
		break;
	default:
		Error::mess << "cell type invalid: " << celltype;
		ERRSET();
		break;
	}
	
	FuncEnd();
}
/*****************************************************************************
 * ElementBase - second level
 ****************************************************************************/

/** @cond */
template<CellType C>
arma::mat::fixed<Cell<C>::nPoint,Cell<C>::nPoint> ElementBase<C>::matLdCn_;
template<CellType C>
arma::vec::fixed<Cell<C>::nPoint> ElementBase<C>::vecLdCn_[nSafe_];
template<CellType C>
arma::vec::fixed<Cell<C>::nFace> ElementBase<C>::vecLdFc_[nSafe_];    
template<CellType C>
arma::ivec::fixed<Cell<C>::nPoint> ElementBase<C>::ivecLdCn_;
template<CellType C>
arma::mat::fixed<2, Cell<C>::nPoint> ElementBase<C>::KD_;              
template<CellType C>
arma::rowvec::fixed<Cell<C>::nPoint> ElementBase<C>::V_;             
/** @endcond */

template<CellType C>
inline CellType ElementBase<C>::cellType() const{
	FuncBegin();
	return C;
	FuncEnd();
}

template<CellType C>
inline int ElementBase<C>::nNode() const{
	FuncBegin();
	return Cell<C>::nPoint;
	FuncEnd();
}

template<CellType C>
inline int ElementBase<C>::nFace() const{
	FuncBegin();
	return Cell<C>::nFace;
	FuncEnd();
}

template<CellType C>
inline const arma::ivec& ElementBase<C>::idxGlob(){
	FuncBegin();
	for (int i = 0 ; i < nNode() ; i++) ivecLdCn_(i) = nd_[i]->idx;
	return ivecLdCn_;
	FuncEnd();
}

template<CellType C>
inline const arma::vec& ElementBase<C>::lDatCnDis (const std::vector<double> &dat, const uint i){
	FuncBegin();
	fml::chkIdx(i, nSafe_);
	fml::locFromGlob(dat.size(), dat, nNode(), dd_, vecLdCn_[i]);
	return vecLdCn_[i];
	FuncEnd();
}

template<CellType C>
inline const arma::vec& ElementBase<C>::lDatCnCon (const std::vector<double> &dat, const uint i){
	FuncBegin();
	fml::chkIdx(i, nSafe_);
	fml::locFromGlob(dat.size(), dat, nNode(), nd_, vecLdCn_[i]);
	return vecLdCn_[i];
	FuncEnd();
}

template<CellType C>
inline const arma::vec& ElementBase<C>::lDatCnCon (const double* dat, const uint i){
	FuncBegin();
	fml::chkIdx(i, nSafe_);
	fml::locFromGlob(-1, dat, nNode(), nd_, vecLdCn_[i]); //no bound checking for arrays
	return vecLdCn_[i];
	FuncEnd();
}

template<CellType C>
inline const arma::vec& ElementBase<C>::lDatUpDis (const std::vector<double> &dat, const int idx[],const uint i){
	FuncBegin();
	fml::chkIdx(i, nSafe_);
	for (int j = 0 ; j < nFace() ; j++)
		vecLdFc_[i](j) = dat.at( dd_[ idx[j] ]->idx );
	return vecLdFc_[i];
	FuncEnd();
}

template<CellType C>
inline void ElementBase<C>::constructDuplData(){
	FuncBegin();
	for(int i = 0 ; i < nNode() ; i++) dd_[i] = nd_[i]->addRegion(reg_);
	FuncEnd();
}

template<CellType C>
inline ElementBase<C>::ElementBase(RegionPorous* reg, Node *nd[]){
	FuncBegin();
	reg_ = (typename Cell<C>::Reg *) reg;
	for(int i = 0 ; i < nNode() ; i++){
		nd_[i] = nd[i];
		dd_[i] = (DuplData*) NULL;
	}	
	for(int i = 0 ; i < nFace() ; i++){
		upwetidx_[i] = i;
		upnonidx_[i] = i;
	}
	FuncEnd();
}

template<CellType C>
inline std::string ElementBase<C>::nameBase(const std::vector<double> *S,
											const double *P, const std::vector<double> *Pc,
											const std::vector<double> *Lw, const std::vector<double> *Ln){
	FuncBegin();
	std::stringstream ss;
	ss << "Region: " << reg_->idx << std::endl
	   << "Nodes: " << idxGlob().t();
	ss << "DuplData: ";
	for(int i = 0 ; i < nNode() ; i++) if(dd_[i]) ss << dd_[i]->idx << " " << dd_[i]->reg->idx << " | ";
	ss << std::endl << "UpNodeWetting: ";
	for(int i = 0 ; i < nFace() ; i++) ss << nd_[ upwetidx_[i] ]->idx << " ";
	ss << std::endl;
	ss << std::endl << "UpNodeNonWetting: ";
	for(int i = 0 ; i < nFace() ; i++) ss << nd_[ upnonidx_[i] ]->idx << " ";
	ss << std::endl;
	if(P)  ss << "p_corner:" << lDatCnCon(P, nSafe_-1).t();
	if(Pc) ss << "pc_corner:" << lDatCnDis(*Pc, nSafe_-1).t();
	if(S)  ss << "s_corner:"  << lDatCnDis(*S , nSafe_-1).t();
	if(Lw){
		ss << "lw_corner"   << lDatCnDis(*Lw , nSafe_-1).t();
		ss << "lw_up:"  << lDatUpDis(*Lw, upwetidx_, nSafe_-1).t();
	}
	if(Ln) {
		ss << "ln_corner:"  << lDatCnDis(*Ln , nSafe_-1).t();
		ss << "ln_up:"  << lDatUpDis(*Ln , upnonidx_, nSafe_-1).t();
	}

	return ss.str();
	FuncEnd();
}

template<CellType C>
int ElementBase<C>::regionID() const {
	FuncBegin();
	return reg_->ID ;
	FuncEnd();
}
template<CellType C>
const DuplData* ElementBase<C>::dupl(const int i) const{
	FuncBegin();
	fml::chkIdx(i, nNode());
	return dd_[i];
	FuncEnd();
}

/*****************************************************************************
 * ElementPoly n>=3 - third level
 ****************************************************************************/
/** @cond */
template<CellType C>
arma::rowvec::fixed< Cell<C>::nPoint > ElementPoly<C>::N_;
template<CellType C>
arma::mat::fixed <Cell<C>::nPoint, 2> ElementPoly<C>::B_;
template<CellType C>
arma::mat::fixed <2, Cell<C>::nPoint> ElementPoly<C>::X_;
template<CellType C>
arma::mat::fixed <2, 2> ElementPoly<C>::J_;                      
/** @endcond */

template < CellType C>
inline const arma::mat& ElementPoly<C>::matKD(){
	FuncBegin();

	matJ( Cell<C>::center(0), Cell<C>::center(1) ); // B is found implicitly
	f::KD_ = f::reg_->k * arma::trans( B_ * arma::inv(J_) );
	
	return f::KD_;	
	FuncEnd();
}


template < CellType C>
inline const arma::rowvec& ElementPoly<C>::matVolume(){
	FuncBegin();
	
	for (int i = 0 ; i < f::nNode() ; i++){
		f::V_(i) =	Cell<C>::rawVol * det ( matJ( Cell<C>::vIp(0,i), Cell<C>::vIp(1,i) ) ) / f::nNode();
	}
	
	return f::V_;	
	FuncEnd();
}

template < CellType C>
inline void ElementPoly<C>::fndUpW(double const *P){
	FuncBegin();
	
	const arma::vec &ploc = f::lDatCnCon(P,f::nSafe_-1);    
	for (int i = 0 ; i < f::nNode() ; i++){
		f::upwetidx_[i] = ( arma::as_scalar(H_.row(i) * (ploc)) < 0 ?
							i                                       :
							Cell<C>::idxMin1(i)                       );
	}
	
	FuncEnd();
}

template < CellType C>
inline void ElementPoly<C>::fndUpN(double const *P, const std::vector<double> &Pc){
	FuncBegin();
	
	arma::vec::fixed< Cell<C>::nPoint > ploc = f::lDatCnCon(P,f::nSafe_-1) + f::lDatCnDis(Pc,f::nSafe_-2);    
	for (int i = 0 ; i < f::nNode() ; i++){
		f::upnonidx_[i] =( arma::as_scalar(H_.row(i) * (ploc)) < 0 ?
						   i                                       :
						   Cell<C>::idxMin1(i)                       );
	}
	
	FuncEnd();
}

template < CellType C>
inline const arma::mat& ElementPoly<C>::lhsP(const std::vector<double> &Lw,
									  const std::vector<double> &Ln){
	FuncBegin();

	arma::vec::fixed<Cell<C>::nFace> lloc =
		f::lDatUpDis(Lw,f::upwetidx_,f::nSafe_ - 1) + f::lDatUpDis(Ln,f::upnonidx_,f::nSafe_ - 2);

	for (int i = 0 ; i < f::nNode() ; i++)
		f::matLdCn_.row(i) =
			lloc(i) * H_.row(i) -
			lloc( Cell<C>::idxPlus1(i) ) * H_.row( Cell<C>::idxPlus1(i) );

	return f::matLdCn_;
	FuncEnd();
}

template < CellType C>
inline const arma::vec& ElementPoly<C>::rhsP(const std::vector<double> &Ln,
											 const std::vector<double> &Pc,
											 const uint i){
	FuncBegin();

	fml::chkIdx(i, f::nSafe_- 1);
    const arma::vec &lnloc = f::lDatUpDis(Ln, f::upnonidx_ , f::nSafe_ - 1); 
	const arma::vec &pcloc = f::lDatCnDis(Pc , f::nSafe_ - 1);
	
	for (int j = 0 ; j < f::nNode() ; j++){
		f::vecLdCn_[i](j) = -arma::as_scalar ( ( lnloc(j) * H_.row(j) -
												 lnloc(Cell<C>::idxPlus1(j)) * H_.row(Cell<C>::idxPlus1(j) ) ) *
											   pcloc );	
	}

	return f::vecLdCn_[i];
	FuncEnd();
}

template < CellType C>
inline const arma::vec& ElementPoly<C>::rhsS(const std::vector<double> &Lw,
									  const double *P,
									  const uint i){
	FuncBegin();

	fml::chkIdx(i, f::nSafe_ - 1);
	const arma::vec &lwloc = f::lDatUpDis(Lw, f::upwetidx_, f::nSafe_ - 1);
	const arma::vec &ploc = f::lDatCnCon(P , f::nSafe_ - 1);
	for (int j = 0 ; j < f::nNode() ; j++)
		f::vecLdCn_[i](j) = arma::as_scalar (
			( lwloc( j )                   * H_.row( j ) -
			  lwloc( Cell<C>::idxPlus1(j) )* H_.row( Cell<C>::idxPlus1(j) ) ) *
			ploc );

	return f::vecLdCn_[i];
	FuncEnd();
}

template < CellType C>
inline const arma::rowvec& ElementPoly<C>::matN(const double z, const double e){
	FuncBegin();
	
	Cell<C>::N(z,e,N_);
	return N_;
		
	FuncEnd();
}

template < CellType C>
inline const arma::mat& ElementPoly<C>::matB(const double z, const double e){
	FuncBegin();

	Cell<C>::B(z,e,B_);
	return B_;
	
	FuncEnd();
}

template < CellType C>
inline const arma::mat& ElementPoly<C>::matX(){
	FuncBegin();

	for (int i = 0 ; i < f::nNode() ; i++) {
		X_(0,i) = f::nd_[i]->x;
		X_(1,i) = f::nd_[i]->y;
	}

	return X_;
	FuncEnd();
}

template < CellType C>
inline const arma::mat& ElementPoly<C>::matJ(const double z, const double e){
	FuncBegin();

	J_ = matX() * matB(z,e);
	
	return J_;
	FuncEnd();
}

template < CellType C>
inline ElementPoly<C>::ElementPoly(RegionPorous* reg, Node *nd[]):ElementBase<C>(reg,nd){}

template < CellType C>
inline void ElementPoly<C>::constructBVertices(int *res){
	FuncBegin();

	int res_bup[Cell<C>::nPoint*2]; //just for enabling to pass null as res
	if (!res) res = res_bup;        // use an if instead of this idiotic method

	for (int i = 0 ; i < f::nNode(); i++){
		if (f::nd_[i]->bvertex){
			res[2*i] = f::nd_[i]->bvertex->checkPre   ( f::nd_[ Cell<C>::idxMin1(i) ] );
			res[2*i+1] = f::nd_[i]->bvertex->checkNext( f::nd_[ Cell<C>::idxPlus1(i) ] );
		}
	}

	FuncEnd();
}

template < CellType C>
inline void ElementPoly<C>::constructGeoParams(){
	FuncBegin();

	for (int i = 0 ; i < f::nNode() ; i++){
		matJ( Cell<C>::fIp(0,i) , Cell<C>::fIp(1,i) ); //mat b is found automatically
		H_.row(i) = arma::trans( B_ * arma::inv(J_) * f::reg_->k.t() * fml::Rot * J_ * Cell<C>::del.col(i) );
	}
	
	FuncEnd();
}

template < CellType C>
inline std::string ElementPoly<C>::name(const std::vector<double> *S,
										const double *P, const std::vector<double> *Pc,
										const std::vector<double> *Lw, const std::vector<double> *Ln){
	FuncBegin();

	std::stringstream ss;

	ss << "ElementPoly<" << (int)C << ">:" << std::endl
	   << f::nameBase(S,P,Pc,Lw,Ln) << "H:\n" << H_ << "V:\n" << matVolume() << "KD:\n" << matKD();
	
	return ss.str();
	FuncEnd();
}


/*****************************************************************************
 * ElementPoly n=2 fracture - third level
 ****************************************************************************/
inline const arma::mat& elefrac::matKD(){
	FuncBegin();

	double dx, dy, l;
	l = fml::lineLength(nd_[0]->x, nd_[0]->y, nd_[1]->x, nd_[1]->y);
	dx = nd_[1]->x - nd_[0]->x;
	dy = nd_[1]->y - nd_[0]->y;
	KD_(0,0) = - ( KD_(0,1) = reg_->k * dx / l / l );
	KD_(1,0) = - ( KD_(1,1) = reg_->k * dy / l / l );
	
	return KD_;
	FuncEnd();
}

inline const arma::rowvec& elefrac::matVolume(){
	FuncBegin();

	V_(0) = V_(1) = .5 * fml::lineLength(nd_[0]->x, nd_[0]->y, nd_[1]->x, nd_[1]->y) * reg_->e;
	
	return V_;
	FuncEnd();
}


inline void elefrac::fndUpW(double const *P){
	FuncBegin();
	upwetidx_[0] = ( P[nd_[1]->idx] > P[nd_[0]->idx] ? 1 : 0 );
	FuncEnd();
}

inline void elefrac::fndUpN(double const *P, const std::vector<double> &Pc){
	FuncBegin();
	upnonidx_[0] = ( P[nd_[1]->idx] + Pc.at(dd_[1]->idx) > P[nd_[0]->idx] + Pc.at(dd_[0]->idx)
					 ? 1 : 0 );
	FuncEnd();
}

inline const arma::mat& elefrac::lhsP(const std::vector<double> &Lw,
							   const std::vector<double> &Ln){
	FuncBegin();

	arma::vec::fixed<cellin::nFace> lloc = lDatUpDis(Lw, upwetidx_,nSafe_-1) + lDatUpDis(Ln, upnonidx_, nSafe_- 2);
	matLdCn_(0,0) = matLdCn_(1,1) =
		- ( matLdCn_(0,1) = matLdCn_(1,0) = lloc(0) * KE_L_ );
	
	
	return matLdCn_;
	FuncEnd();
}

inline const arma::vec& elefrac::rhsP(const std::vector<double> &Ln,
							   const std::vector<double> &Pc,
							   const uint i){
	FuncBegin();

	fml::chkIdx (i , nSafe_ - 1 );
	const arma::vec &pcloc = lDatCnDis(Pc, nSafe_-1); 
	const arma::vec &lnloc = lDatUpDis(Ln, upnonidx_, nSafe_-1); 
	vecLdCn_[i](0) = -( vecLdCn_[i](1) = lnloc(0) * KE_L_ * ( pcloc(1) - pcloc(0) ) );

	return vecLdCn_[i];
	FuncEnd();
}

inline const arma::vec& elefrac::rhsS(const std::vector<double> &Lw,
							   const double *P,
							   const uint i){
	FuncBegin();

	fml::chkIdx (i , nSafe_ - 1 );
	const arma::vec &ploc = lDatCnCon(P, nSafe_-1); 
	const arma::vec &lwloc = lDatUpDis(Lw, upwetidx_, nSafe_-1); 
	vecLdCn_[i](0) = -( vecLdCn_[i](1) = lwloc(0) * KE_L_ * ( ploc(0) - ploc(1) ) );
	
	return vecLdCn_[i];
	FuncEnd();
}

inline elefrac::ElementPoly(RegionPorous* reg, Node *nd[]):ElementBase<CellLine>(reg,nd), KE_L_(0) {}

inline void elefrac::constructGeoParams(){
	FuncBegin();
	KE_L_ =  reg_->k * reg_->e / fml::lineLength(nd_[0]->x, nd_[0]->y, nd_[1]->x, nd_[1]->y) ;
	FuncEnd();
}

inline std::string elefrac::name(const std::vector<double> *S,
								 const double *P, const std::vector<double> *Pc,
								 const std::vector<double> *Lw, const std::vector<double> *Ln){
	FuncBegin();

	std::stringstream ss;
	ss << "ElementPoly<" << (int)CellLine << "> ie. elefrac: " << std::endl
	   << nameBase(S,P,Pc,Lw,Ln)  << "KE_L: " << KE_L_ << std::endl
	   << "V:\n" << matVolume() << "KD:\n" << matKD();
	return ss.str();
	FuncEnd();
}


#endif /*ELEMENT_HPP*/

