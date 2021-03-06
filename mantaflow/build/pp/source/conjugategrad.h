




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/conjugategrad.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Conjugate gradient solver
 *
 ******************************************************************************/

#ifndef _CONJUGATEGRADIENT_H
#define _CONJUGATEGRADIENT_H

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"

namespace Manta { 

static const bool CG_DEBUG = false;

//! Basic CG interface 
class GridCgInterface {
	public:
		enum PreconditionType { PC_None=0, PC_ICP, PC_mICP };
        
        GridCgInterface() : mUseResNorm(true) {};
		virtual ~GridCgInterface() {};

		// solving functions
		virtual bool iterate() = 0;
		virtual void solve(int maxIter) = 0;

		// precond
		virtual void setPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) = 0;

		// access
		virtual Real getSigma() const = 0;
		virtual Real getIterations() const = 0;
		virtual Real getResNorm() const = 0;
		virtual void setAccuracy(Real set) = 0;
		virtual Real getAccuracy() const = 0;

		void setUseResNorm(bool set) { mUseResNorm = set; }

	protected:

		// use norm of residual, or max value for threshold?
		bool mUseResNorm; 
};


//! Run single iteration of the cg solver
/*! the template argument determines the type of matrix multiplication,
    typically a ApplyMatrix kernel, another one is needed e.g. for the
    mesh-based wave equation solver */
template<class APPLYMAT>
class GridCg : public GridCgInterface {
	public:
        //! constructor
		GridCg(Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& residual, Grid<Real>& search, FlagGrid& flags, Grid<Real>& tmp, 
				Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk);
        ~GridCg() {}
        
        void doInit();
        bool iterate();
        void solve(int maxIter);
        //! init pointers, and copy values from "normal" matrix
        void setPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak);
        
        // Accessors        
        Real getSigma() const { return mSigma; }
        Real getIterations() const { return mIterations; }

        Real getResNorm() const { return mResNorm; }

        void setAccuracy(Real set) { mAccuracy=set; }
        Real getAccuracy() const { return mAccuracy; }

	protected:
		bool mInited;
		int mIterations;
		// grids
		Grid<Real>& mDst;
		Grid<Real>& mRhs;
		Grid<Real>& mResidual;
		Grid<Real>& mSearch;
		FlagGrid& mFlags;
		Grid<Real>& mTmp;

		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;

		PreconditionType mPcMethod;
		//! preconditioning grids
		Grid<Real> *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk;

		//! sigma / residual
		Real mSigma;
		//! accuracy of solver (max. residuum)
		Real mAccuracy;
		//! norm of the residual
		Real mResNorm;
}; // GridCg


//! Kernel: Apply symmetric stored Matrix



struct ApplyMatrix : public KernelBase {  ApplyMatrix (FlagGrid& _flags, Grid<Real >& _dst, Grid<Real >& _src, Grid<Real >& _A0, Grid<Real >& _Ai, Grid<Real >& _Aj, Grid<Real >& _Ak) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_dst(_dst), m_src(_src), m_A0(_A0), m_Ai(_Ai), m_Aj(_Aj), m_Ak(_Ak) { run(); } ApplyMatrix (const ApplyMatrix& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_dst(o.m_dst), m_src(o.m_src), m_A0(o.m_A0), m_Ai(o.m_Ai), m_Aj(o.m_Aj), m_Ak(o.m_Ak) {} inline void op(int idx,FlagGrid& flags, Grid<Real >& dst, Grid<Real >& src, Grid<Real >& A0, Grid<Real >& Ai, Grid<Real >& Aj, Grid<Real >& Ak)  {
    if (!flags.isFluid(idx)) {
        dst[idx] = src[idx];
        return;
    }    
    dst[idx] =  src[idx] * A0[idx]
                + src[idx-X] * Ai[idx-X]
                + src[idx+X] * Ai[idx]
                + src[idx-Y] * Aj[idx-Y]
                + src[idx+Y] * Aj[idx]
                + src[idx-Z] * Ak[idx-Z] 
                + src[idx+Z] * Ak[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_flags, m_dst, m_src, m_A0, m_Ai, m_Aj, m_Ak); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_dst; Grid<Real >& m_src; Grid<Real >& m_A0; Grid<Real >& m_Ai; Grid<Real >& m_Aj; Grid<Real >& m_Ak;  }; 

//! Kernel: Apply symmetric stored Matrix. 2D version



struct ApplyMatrix2D : public KernelBase {  ApplyMatrix2D (FlagGrid& _flags, Grid<Real >& _dst, Grid<Real >& _src, Grid<Real >& _A0, Grid<Real >& _Ai, Grid<Real >& _Aj, Grid<Real >& _Ak) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_dst(_dst), m_src(_src), m_A0(_A0), m_Ai(_Ai), m_Aj(_Aj), m_Ak(_Ak) { run(); } ApplyMatrix2D (const ApplyMatrix2D& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_dst(o.m_dst), m_src(o.m_src), m_A0(o.m_A0), m_Ai(o.m_Ai), m_Aj(o.m_Aj), m_Ak(o.m_Ak) {} inline void op(int idx,FlagGrid& flags, Grid<Real >& dst, Grid<Real >& src, Grid<Real >& A0, Grid<Real >& Ai, Grid<Real >& Aj, Grid<Real >& Ak)  {
    unusedParameter(Ak); // only there for parameter compatibility with ApplyMatrix
    
    if (!flags.isFluid(idx)) {
        dst[idx] = src[idx];
        return;
    }    
    dst[idx] =  src[idx] * A0[idx]
                + src[idx-X] * Ai[idx-X]
                + src[idx+X] * Ai[idx]
                + src[idx-Y] * Aj[idx-Y]
                + src[idx+Y] * Aj[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_flags, m_dst, m_src, m_A0, m_Ai, m_Aj, m_Ak); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_dst; Grid<Real >& m_src; Grid<Real >& m_A0; Grid<Real >& m_Ai; Grid<Real >& m_Aj; Grid<Real >& m_Ak;  }; 

//! Kernel: Construct the matrix for the poisson equation

struct MakeLaplaceMatrix : public KernelBase {  MakeLaplaceMatrix (FlagGrid& _flags, Grid<Real >& _A0, Grid<Real >& _Ai, Grid<Real >& _Aj, Grid<Real >& _Ak) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_A0(_A0), m_Ai(_Ai), m_Aj(_Aj), m_Ak(_Ak) { run(); } MakeLaplaceMatrix (const MakeLaplaceMatrix& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_A0(o.m_A0), m_Ai(o.m_Ai), m_Aj(o.m_Aj), m_Ak(o.m_Ak) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<Real >& A0, Grid<Real >& Ai, Grid<Real >& Aj, Grid<Real >& Ak)  {
    if (!flags.isFluid(i,j,k))
        return;
    
    // center
    if (!flags.isObstacle(i-1,j,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i+1,j,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j-1,k)) A0(i,j,k) += 1.;
    if (!flags.isObstacle(i,j+1,k)) A0(i,j,k) += 1.;
    if (flags.is3D() && !flags.isObstacle(i,j,k-1)) A0(i,j,k) += 1.;
    if (flags.is3D() && !flags.isObstacle(i,j,k+1)) A0(i,j,k) += 1.;
    
    if (flags.isFluid(i+1,j,k)) Ai(i,j,k) = -1.;
    if (flags.isFluid(i,j+1,k)) Aj(i,j,k) = -1.;
    if (flags.is3D() && flags.isFluid(i,j,k+1)) Ak(i,j,k) = -1.;    
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_A0, m_Ai, m_Aj, m_Ak); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_A0; Grid<Real >& m_Ai; Grid<Real >& m_Aj; Grid<Real >& m_Ak;  }; 




} // namespace

#endif 

