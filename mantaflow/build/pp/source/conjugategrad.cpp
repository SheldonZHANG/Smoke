




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/conjugategrad.cpp"
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

#include "conjugategrad.h"
#include "commonkernels.h"

using namespace std;
namespace Manta {

const int CG_DEBUGLEVEL = 4;
    
//*****************************************************************************
//  Precondition helpers

//! Preconditioning a la Wavelet Turbulence (needs 4 add. grids)
void InitPreconditionIncompCholesky(FlagGrid& flags,
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak,
				Grid<Real>& orgA0, Grid<Real>& orgAi, Grid<Real>& orgAj, Grid<Real>& orgAk) 
{
	// compute IC according to Golub and Van Loan
	A0 = orgA0;
    Ai = orgAi;
    Aj = orgAj;
    Ak = orgAk;
    
    FOR_IJK(A0) {
        if (flags.isFluid(i,j,k)) {
            const int idx = A0.index(i,j,k);
            A0[idx] = sqrt(A0[idx]);
            
            // correct left and top stencil in other entries
            // for i = k+1:n
            //  if (A(i,k) != 0)
            //    A(i,k) = A(i,k) / A(k,k)
            Real invDiagonal = 1.0f / A0[idx];
            Ai[idx] *= invDiagonal;
            Aj[idx] *= invDiagonal;
            Ak[idx] *= invDiagonal;
            
            // correct the right and bottom stencil in other entries
            // for j = k+1:n 
            //   for i = j:n
            //      if (A(i,j) != 0)
            //        A(i,j) = A(i,j) - A(i,k) * A(j,k)
            A0(i+1,j,k) -= square(Ai[idx]);
            A0(i,j+1,k) -= square(Aj[idx]);
            A0(i,j,k+1) -= square(Ak[idx]);
        }
    }
    
    // invert A0 for faster computation later
    InvertCheckFluid (flags, A0);
};

//! Preconditioning using modified IC ala Bridson (needs 1 add. grid)
void InitPreconditionModifiedIncompCholesky2(FlagGrid& flags,
				Grid<Real>&Aprecond, 
				Grid<Real>&A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) 
{
	// compute IC according to Golub and Van Loan
	Aprecond.clear();
    
    FOR_IJK(flags) {
        if (!flags.isFluid(i,j,k)) continue;

        const Real tau = 0.97;
        const Real sigma = 0.25;
            
        // compute modified incomplete cholesky
        Real e = 0.;
        e = A0(i,j,k) 
            - square(Ai(i-1,j,k) * Aprecond(i-1,j,k) )
            - square(Aj(i,j-1,k) * Aprecond(i,j-1,k) )
            - square(Ak(i,j,k-1) * Aprecond(i,j,k-1) ) ;
        e -= tau * (
                Ai(i-1,j,k) * ( Aj(i-1,j,k) + Ak(i-1,j,k) )* square( Aprecond(i-1,j,k) ) +
                Aj(i,j-1,k) * ( Ai(i,j-1,k) + Ak(i,j-1,k) )* square( Aprecond(i,j-1,k) ) +
                Ak(i,j,k-1) * ( Ai(i,j,k-1) + Aj(i,j,k-1) )* square( Aprecond(i,j,k-1) ) +
                0. );

        // stability cutoff
        if(e < sigma * A0(i,j,k))
            e = A0(i,j,k);

        Aprecond(i,j,k) = 1. / sqrt( e );
    }
};

//! Apply WT-style ICP
void ApplyPreconditionIncompCholesky(Grid<Real>& dst, Grid<Real>& Var1, FlagGrid& flags,
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak,
				Grid<Real>& orgA0, Grid<Real>& orgAi, Grid<Real>& orgAj, Grid<Real>& orgAk)
{
	
    // forward substitution        
    FOR_IJK(dst) {
        if (!flags.isFluid(i,j,k)) continue;
        dst(i,j,k) = A0(i,j,k) * (Var1(i,j,k)
                 - dst(i-1,j,k) * Ai(i-1,j,k)
                 - dst(i,j-1,k) * Aj(i,j-1,k)
                 - dst(i,j,k-1) * Ak(i,j,k-1));    
    }
    
    // backward substitution
    FOR_IJK_REVERSE(dst) {
        const int idx = A0.index(i,j,k);
        if (!flags.isFluid(idx)) continue;
        dst[idx] = A0[idx] * ( dst[idx] 
               - dst(i+1,j,k) * Ai[idx]
               - dst(i,j+1,k) * Aj[idx]
               - dst(i,j,k+1) * Ak[idx]);
    }
}

//! Apply Bridson-style mICP
void ApplyPreconditionModifiedIncompCholesky2(Grid<Real>& dst, Grid<Real>& Var1, FlagGrid& flags,
				Grid<Real>& Aprecond, 
				Grid<Real>& A0, Grid<Real>& Ai, Grid<Real>& Aj, Grid<Real>& Ak) 
{
    // forward substitution        
    FOR_IJK(dst) {
        if (!flags.isFluid(i,j,k)) continue;
        const Real p = Aprecond(i,j,k);
        dst(i,j,k) = p * (Var1(i,j,k)
                 - dst(i-1,j,k) * Ai(i-1,j,k) * p
                 - dst(i,j-1,k) * Aj(i,j-1,k) * p
                 - dst(i,j,k-1) * Ak(i,j,k-1) * p);
    }
    
    // backward substitution
    FOR_IJK_REVERSE(dst) {            
        const int idx = A0.index(i,j,k);
        if (!flags.isFluid(idx)) continue;
        const Real p = Aprecond[idx];
        dst[idx] = p * ( dst[idx] 
               - dst(i+1,j,k) * Ai[idx] * p
               - dst(i,j+1,k) * Aj[idx] * p
               - dst(i,j,k+1) * Ak[idx] * p);
    }
}


//*****************************************************************************
// Kernels    

//! Kernel: Compute the dot product between two Real grids
/*! Uses double precision internally */

struct GridDotProduct : public KernelBase {  GridDotProduct (const Grid<Real >& _a, const Grid<Real >& _b) : KernelBase(&_a, 0), parent((_a).getParent()), m_a(_a), m_b(_b), result(0.0) { run(); } GridDotProduct (const GridDotProduct& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_a(o.m_a), m_b(o.m_b), result(0.0) {} inline void op(int idx,const Grid<Real >& a, const Grid<Real >& b, double& result)  {
    result += (a[idx] * b[idx]);    
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_a, m_b, result); } operator double  () { return result; } const double& getRet() const { return result; } const Grid<Real >& getArg0() { return m_a; } typedef Grid<Real > type0; const Grid<Real >& getArg1() { return m_b; } typedef Grid<Real > type1; FluidSolver* parent; const Grid<Real >& m_a; const Grid<Real >& m_b; double result;  }; ;

//! Kernel: compute residual (init) and add to sigma


struct InitSigma : public KernelBase {  InitSigma (FlagGrid& _flags, Grid<Real >& _dst, Grid<Real >& _rhs, Grid<Real >& _temp) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_dst(_dst), m_rhs(_rhs), m_temp(_temp), sigma(0) { run(); } InitSigma (const InitSigma& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_dst(o.m_dst), m_rhs(o.m_rhs), m_temp(o.m_temp), sigma(0) {} inline void op(int idx,FlagGrid& flags, Grid<Real >& dst, Grid<Real >& rhs, Grid<Real >& temp, double& sigma)  {    
    const double res = rhs[idx] - temp[idx]; 
    dst[idx] = (Real)res;

    // only compute residual in fluid region
    if(flags.isFluid(idx)) 
        sigma += res*res;
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_flags, m_dst, m_rhs, m_temp, sigma); } operator double  () { return sigma; } const double& getRet() const { return sigma; } FlagGrid& getArg0() { return m_flags; } typedef FlagGrid type0; Grid<Real >& getArg1() { return m_dst; } typedef Grid<Real > type1; Grid<Real >& getArg2() { return m_rhs; } typedef Grid<Real > type2; FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_dst; Grid<Real >& m_rhs; Grid<Real >& m_temp; double sigma;  }; ;

//! Kernel: update search vector

struct UpdateSearchVec : public KernelBase {  UpdateSearchVec (Grid<Real >& _dst, Grid<Real >& _src, Real _factor) : KernelBase(&_dst, 0), parent((_dst).getParent()), m_dst(_dst), m_src(_src), m_factor(_factor) { run(); } UpdateSearchVec (const UpdateSearchVec& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_dst(o.m_dst), m_src(o.m_src), m_factor(o.m_factor) {} inline void op(int idx,Grid<Real >& dst, Grid<Real >& src, Real factor)  {
    dst[idx] = src[idx] + factor * dst[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_dst, m_src, m_factor); } FluidSolver* parent; Grid<Real >& m_dst; Grid<Real >& m_src; Real m_factor;  }; 

//*****************************************************************************
//  CG class

template<class APPLYMAT>
GridCg<APPLYMAT>::GridCg(Grid<Real>& dst, Grid<Real>& rhs, Grid<Real>& residual, Grid<Real>& search, FlagGrid& flags, Grid<Real>& tmp, 
               Grid<Real>* pA0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk) :
    GridCgInterface(), mInited(false), mIterations(0), mDst(dst), mRhs(rhs), mResidual(residual),
    mSearch(search), mFlags(flags), mTmp(tmp), mpA0(pA0), mpAi(pAi), mpAj(pAj), mpAk(pAk),
    mPcMethod(PC_None), mpPCA0(pA0), mpPCAi(pAi), mpPCAj(pAj), mpPCAk(pAk), mSigma(0.), mAccuracy(VECTOR_EPSILON), mResNorm(1e20) 
{
    dst.clear();
    residual.clear();
    search.clear();
    tmp.clear();            
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::doInit() {
    mInited = true;

    mResidual = mRhs; // p=0, residual = b
    
    if (mPcMethod == PC_ICP) {
        assertMsg(mDst.is3D(), "ICP only supports 3D grids so far");
        InitPreconditionIncompCholesky(mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
        ApplyPreconditionIncompCholesky(mTmp, mResidual, mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
    } else if (mPcMethod == PC_mICP) {
        assertMsg(mDst.is3D(), "mICP only supports 3D grids so far");
        InitPreconditionModifiedIncompCholesky2(mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
        ApplyPreconditionModifiedIncompCholesky2(mTmp, mResidual, mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
    } else {
        mTmp = mResidual;
    }
    
    mSearch = mTmp;
    
    mSigma = GridDotProduct(mTmp, mResidual);    
}

template<class APPLYMAT>
bool GridCg<APPLYMAT>::iterate() {
    if(!mInited) doInit();

    mIterations++;

    // create matrix application operator passed as template argument,
    // this could reinterpret the mpA pointers (not so clean right now)
    // tmp = applyMat(search)
    
    APPLYMAT (mFlags, mTmp, mSearch, *mpA0, *mpAi, *mpAj, *mpAk);
    
    // alpha = sigma/dot(tmp, search)
    Real dp = GridDotProduct(mTmp, mSearch);
    Real alpha = 0.;
    if(fabs(dp)>0.) alpha = mSigma / (Real)dp;
    
    mDst.scaledAdd(mSearch, alpha);    // dst += search * alpha
    mResidual.scaledAdd(mTmp, -alpha); // residual += tmp * -alpha
    
    if (mPcMethod == PC_ICP)
        ApplyPreconditionIncompCholesky(mTmp, mResidual, mFlags, *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk, *mpA0, *mpAi, *mpAj, *mpAk);
    else if (mPcMethod == PC_mICP)
        ApplyPreconditionModifiedIncompCholesky2(mTmp, mResidual, mFlags, *mpPCA0, *mpA0, *mpAi, *mpAj, *mpAk);
    else
        mTmp = mResidual;
        
    // compute norm of the residual?
    if(this->mUseResNorm) { 
        mResNorm = GridSumSqr(mResidual).sum;
    } else {
        // TODO test...
        mResNorm = mResidual.getMaxAbsValue();        
    }
    if(mIterations % 10 == 9) debMsg("GridCg::Iteration i="<<mIterations<<", resNorm="<<mResNorm<<" accuracy="<<mAccuracy, 1);

    // abort here to safe some work...
    if(mResNorm<mAccuracy) return false;

    Real sigmaNew = GridDotProduct(mTmp, mResidual);
    Real beta = sigmaNew / mSigma;
    
    // search =  tmp + beta * search
    UpdateSearchVec (mSearch, mTmp, beta);

    debMsg("PB-Cg::iter2 i="<<mIterations<<" sigma="<<mSigma<<" alpha="<<alpha<<" beta="<<beta<<" ", CG_DEBUGLEVEL);
    mSigma = sigmaNew;
    
    /* debMsg("PB-CG-Norms::p"<<sqrt( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<" search"<<sqrt( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue(), CG_DEBUGLEVEL ) 
            <<" res"<<sqrt( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) <<" tmp"<<sqrt( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() ), CG_DEBUGLEVEL ); // debug
    if(CG_DEBUG) debMsg("PB-CG-Norms","p"<<( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<" search"<<( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue(), CG_DEBUGLEVEL ) 
            <<" res"<<( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) <<" tmp"<<( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() ), CG_DEBUGLEVEL ); // debug, no sqrt!*/
    return true;
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::solve(int maxIter) {
    for (int iter=0; iter<maxIter; iter++) {
        if (!iterate()) iter=maxIter;
    } 
    return;
}

template<class APPLYMAT>
void GridCg<APPLYMAT>::setPreconditioner(PreconditionType method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) {
    mPcMethod = method;
    mpPCA0 = A0;
    mpPCAi = Ai;
    mpPCAj = Aj;
    mpPCAk = Ak;
}

// explicit instantiation
template class GridCg<ApplyMatrix>;
template class GridCg<ApplyMatrix2D>;

}; // DDF


