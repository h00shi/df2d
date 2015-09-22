/* @file mdata.cpp
 * .cpp file for mdata.
 */

#include "mdata.hpp"
#include "error.hpp"

void MData::initialize(){
	FuncBegin();
	
	ksp = (KSP) NULL;
    Pvec = (Vec) NULL;
	P = (double*) NULL;
	A = (Mat) NULL;
	b = (Vec) NULL;
	J = (JFunc*) NULL;
	qIn = qOut = qWin = qWout = 0;
	nIt = 0;
	cT = 0;

	FuncEnd();
}

void MData::finalize(){
	FuncBegin();

	if (ksp) KSPDestroy(&ksp);
	if (Pvec){
		if (P) VecRestoreArrayRead(Pvec, &P);
		VecDestroy(&Pvec);
	}
	if (A) MatDestroy(&A);
	if (b) VecDestroy(&b);
	if (J)	delete J ;

  FuncEnd();
}
