




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/pressure.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction:
 * - solve_pressure
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "kernel.h"
#include "conjugategrad.h"

using namespace std;
namespace Manta {

//! Kernel: Construct the right-hand side of the poisson equation



struct MakeRhs : public KernelBase {  MakeRhs (FlagGrid& _flags, Grid<Real >& _rhs, MACGrid& _vel, Grid<Real >* _perCellCorr) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_rhs(_rhs), m_vel(_vel), m_perCellCorr(_perCellCorr), cnt(0), sum(0) { run(); } MakeRhs (const MakeRhs& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_rhs(o.m_rhs), m_vel(o.m_vel), m_perCellCorr(o.m_perCellCorr), cnt(0), sum(0) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<Real >& rhs, MACGrid& vel, Grid<Real >* perCellCorr, int& cnt, double& sum)  {
    if (!flags.isFluid(i,j,k)) {
        rhs(i,j,k) = 0;
        return;
    }
       
    // compute divergence 
    // assumes vel at obstacle interfaces is set to zero
    /* Real set = 0;
    if (!flags.isObstacle(i-1,j,k)) set += vel(i,j,k).x;
    if (!flags.isObstacle(i+1,j,k)) set -= vel(i+1,j,k).x;
    if (!flags.isObstacle(i,j-1,k)) set += vel(i,j,k).y;
    if (!flags.isObstacle(i,j+1,k)) set -= vel(i,j+1,k).y;
    if (!flags.isObstacle(i,j,k-1)) set += vel(i,j,k).z;
    if (!flags.isObstacle(i,j,k+1)) set -= vel(i,j,k+1).z; */
    Real set = vel(i,j,k).x - vel(i+1,j,k).x + vel(i,j,k).y - vel(i,j+1,k).y + vel(i,j,k).z - vel(i,j,k+1).z;
    
    // per cell divergence correction
    if(perCellCorr) 
        set += perCellCorr->get(i,j,k);
    
    // obtain sum, cell count
    sum += set;
    cnt++;
    
    rhs(i,j,k) = set;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_rhs, m_vel, m_perCellCorr, cnt, sum); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_rhs; MACGrid& m_vel; Grid<Real >* m_perCellCorr; int cnt; double sum;  }; 

//! Kernel: Apply velocity update from poisson equation


struct CorrectVelocity : public KernelBase {  CorrectVelocity (FlagGrid& _flags, MACGrid& _vel, Grid<Real >& _pressure) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_pressure(_pressure) { run(); } CorrectVelocity (const CorrectVelocity& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_pressure(o.m_pressure) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Grid<Real >& pressure)  {
    // correct all faces between fluid-fluid and fluid-empty cells
	// skip everything with obstacles...
	if (flags.isObstacle(i,j,k))
        return;
    
    // skip faces between two empty cells
	const bool curEmpty = flags.isEmpty(i,j,k);
    const Real p = pressure(i,j,k);

    if (!curEmpty || !flags.isEmpty(i-1,j,k)) vel(i,j,k).x -= (p - pressure(i-1,j,k));
    if (!curEmpty || !flags.isEmpty(i,j-1,k)) vel(i,j,k).y -= (p - pressure(i,j-1,k));
    if (!curEmpty || !flags.isEmpty(i,j,k-1)) vel(i,j,k).z -= (p - pressure(i,j,k-1));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_pressure); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Grid<Real >& m_pressure;  }; 

//! Kernel: Set matrix stencils and velocities to enable open boundaries


struct SetOpenBound : public KernelBase {  SetOpenBound (Grid<Real >& _A0, Grid<Real >& _Ai, Grid<Real >& _Aj, Grid<Real >& _Ak, MACGrid& _vel, Vector3D<bool > _lowerBound, Vector3D<bool > _upperBound) : KernelBase(&_A0, 0), parent((_A0).getParent()), m_A0(_A0), m_Ai(_Ai), m_Aj(_Aj), m_Ak(_Ak), m_vel(_vel), m_lowerBound(_lowerBound), m_upperBound(_upperBound) { run(); } SetOpenBound (const SetOpenBound& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_A0(o.m_A0), m_Ai(o.m_Ai), m_Aj(o.m_Aj), m_Ak(o.m_Ak), m_vel(o.m_vel), m_lowerBound(o.m_lowerBound), m_upperBound(o.m_upperBound) {} inline void op(int i,int j,int k,Grid<Real >& A0, Grid<Real >& Ai, Grid<Real >& Aj, Grid<Real >& Ak, MACGrid& vel, Vector3D<bool > lowerBound, Vector3D<bool > upperBound)  {    
    // set velocity boundary conditions
    if (lowerBound.x && i == 0) vel(0,j,k) = vel(1,j,k);
    if (lowerBound.y && j == 0) vel(i,0,k) = vel(i,1,k);
    if (lowerBound.z && k == 0) vel(i,j,0) = vel(i,j,1);
    if (upperBound.x && i == maxX-1) vel(maxX-1,j,k) = vel(maxX-2,j,k);
    if (upperBound.y && j == maxY-1) vel(i,maxY-1,k) = vel(i,maxY-2,k);
    if (upperBound.z && k == maxZ-1) vel(i,j,maxZ-1) = vel(i,j,maxZ-2);
    
    // set matrix stencils at boundary
    if ((lowerBound.x && i<=1) || (upperBound.x && i>=maxX-2) ||
        (lowerBound.y && j<=1) || (upperBound.y && j>=maxY-2) ||
        (lowerBound.z && k<=1) || (upperBound.z && k>=maxZ-2)) {
        A0(i,j,k) = vel.is3D() ? 6.0 : 4.0;
        Ai(i,j,k) = -1.0;
        Aj(i,j,k) = -1.0;
        if (vel.is3D()) Ak(i,j,k) = -1.0;
    }
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k, m_A0, m_Ai, m_Aj, m_Ak, m_vel, m_lowerBound, m_upperBound); } FluidSolver* parent; Grid<Real >& m_A0; Grid<Real >& m_Ai; Grid<Real >& m_Aj; Grid<Real >& m_Ak; MACGrid& m_vel; Vector3D<bool > m_lowerBound; Vector3D<bool > m_upperBound;  }; 

//! Kernel: Set matrix rhs for outflow

struct SetOutflow : public KernelBase {  SetOutflow (Grid<Real >& _rhs, Vector3D<bool > _lowerBound, Vector3D<bool > _upperBound, int _height) : KernelBase(&_rhs, 0), parent((_rhs).getParent()), m_rhs(_rhs), m_lowerBound(_lowerBound), m_upperBound(_upperBound), m_height(_height) { run(); } SetOutflow (const SetOutflow& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_rhs(o.m_rhs), m_lowerBound(o.m_lowerBound), m_upperBound(o.m_upperBound), m_height(o.m_height) {} inline void op(int i,int j,int k,Grid<Real >& rhs, Vector3D<bool > lowerBound, Vector3D<bool > upperBound, int height)  {
    if ((lowerBound.x && i < height) || (upperBound.x && i >= maxX-1-height) ||
        (lowerBound.y && j < height) || (upperBound.y && j >= maxY-1-height) ||
        (lowerBound.z && k < height) || (upperBound.z && k >= maxZ-1-height))
        rhs(i,j,k) = 0;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k, m_rhs, m_lowerBound, m_upperBound, m_height); } FluidSolver* parent; Grid<Real >& m_rhs; Vector3D<bool > m_lowerBound; Vector3D<bool > m_upperBound; int m_height;  }; 


// *****************************************************************************
// Ghost fluid helpers
// TODO set sides individually?

// iso surface level, usually zero
static const int LEVELSET_ISOSURFACE = 0.;

static inline Real getGhostMatrixAddition(Real a, Real b, const Real accuracy) {
	Real ret = 0.f;

	if(a < 0 && b < 0)
		ret = 1;
	else if( (a >= 0) && (b < 0))
		ret = b / (b - a);
	else if ( (a < 0) && (b >= 0))
		ret = a / (a - b);
	else
		ret = 0;

	if(ret < accuracy)
		ret = accuracy;

	Real invret = 1./ret;
	return invret;
}

//! Kernel: Adapt A0 for ghost fluid


struct ApplyGhostFluid : public KernelBase {  ApplyGhostFluid (FlagGrid& _flags, Grid<Real >& _phi, Grid<Real >& _A0, Real _accuracy) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_phi(_phi), m_A0(_A0), m_accuracy(_accuracy) { run(); } ApplyGhostFluid (const ApplyGhostFluid& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_phi(o.m_phi), m_A0(o.m_A0), m_accuracy(o.m_accuracy) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<Real >& phi, Grid<Real >& A0, Real accuracy)  {
    if (!flags.isFluid(i,j,k))
        return;
    
    const Real curPhi = phi(i,j,k);
    
    if (flags.isEmpty(i-1,j,k)) A0(i,j,k) += getGhostMatrixAddition( phi(i-1,j,k), curPhi, accuracy);
    if (flags.isEmpty(i+1,j,k)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i+1,j,k), accuracy);
    if (flags.isEmpty(i,j-1,k)) A0(i,j,k) += getGhostMatrixAddition( phi(i,j-1,k), curPhi, accuracy);
    if (flags.isEmpty(i,j+1,k)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i,j+1,k), accuracy);
    if (flags.is3D() && flags.isEmpty(i,j,k-1)) A0(i,j,k) += getGhostMatrixAddition( phi(i,j,k-1), curPhi, accuracy);
    if (flags.is3D() && flags.isEmpty(i,j,k+1)) A0(i,j,k) += getGhostMatrixAddition( curPhi, phi(i,j,k+1), accuracy);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_phi, m_A0, m_accuracy); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_phi; Grid<Real >& m_A0; Real m_accuracy;  }; 

//! Kernel: Correct velocities for ghost fluids


struct CorrectVelGhostFluid : public KernelBase {  CorrectVelGhostFluid (FlagGrid& _flags, MACGrid& _vel, Grid<Real >& _pressure) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_pressure(_pressure) { run(); } CorrectVelGhostFluid (const CorrectVelGhostFluid& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_pressure(o.m_pressure) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Grid<Real >& pressure)  {
    bool curFluid = flags.isFluid(i,j,k);
    if (!curFluid && !flags.isEmpty(i,j,k))
        return;
    
    const Real curPress = pressure(i,j,k);

    //const Real curPhi = phi(i,j,k);
	// TODO - include ghost fluid factor  NT_DEBUG
    
    // in contrast to old implementation:
    // make sure to add gradient for all fluid-empty or fluid-fluid combinations
    // of neighbors...

    if (!flags.isObstacle(i-1,j,k) && (curFluid || flags.isFluid(i-1,j,k)))
        vel(i,j,k).x -= curPress - pressure(i-1,j,k);
    
    if (!flags.isObstacle(i,j-1,k) && (curFluid || flags.isFluid(i,j-1,k)))
        vel(i,j,k).y -= curPress - pressure(i,j-1,k);
    
    if (flags.is3D() && (!flags.isObstacle(i,j,k-1) && (curFluid || flags.isFluid(i,j,k-1))))
        vel(i,j,k).z -= curPress - pressure(i,j,k-1);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_pressure); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Grid<Real >& m_pressure;  }; 

inline void convertDescToVec(const string& desc, Vector3D<bool>& lo, Vector3D<bool>& up) {
    for(size_t i=0; i<desc.size(); i++) {
        if (desc[i] == 'x') lo.x = true;
        else if (desc[i] == 'y') lo.y = true;
        else if (desc[i] == 'z') lo.z = true;
        else if (desc[i] == 'X') up.x = true;
        else if (desc[i] == 'Y') up.y = true;
        else if (desc[i] == 'Z') up.z = true;
        else errMsg("invalid character in boundary description string. Only [xyzXYZ] allowed.");
    }
}

//! Perform pressure projection of the velocity grid












void solvePressure(MACGrid& vel,Grid<Real>& pressure,FlagGrid& flags,Grid<Real>* phi = 0,Grid<Real>* perCellCorr = 0,Real ghostAccuracy = 0,Real cgMaxIterFac = 1.5,Real cgAccuracy = 1e-3,string openBound = "",string outflow = "",int outflowHeight = 1,int precondition = 0,bool enforceCompatibility = false,bool useResNorm = true , FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    //assertMsg(vel.is3D(), "Only 3D grids supported so far");
    
    // parse strings
    Vector3D<bool> loOpenBound, upOpenBound, loOutflow, upOutflow;
    convertDescToVec(openBound, loOpenBound, upOpenBound);
    convertDescToVec(outflow, loOutflow, upOutflow);
    if (vel.is2D() && (loOpenBound.z || upOpenBound.z))
        errMsg("open boundaries for z specified for 2D grid");
    
    // reserve temp grids
    Grid<Real> rhs(parent);
    Grid<Real> residual(parent);
    Grid<Real> search(parent);
    Grid<Real> A0(parent);
    Grid<Real> Ai(parent);
    Grid<Real> Aj(parent);
    Grid<Real> Ak(parent);
    Grid<Real> tmp(parent);
    Grid<Real> pca0(parent);
    Grid<Real> pca1(parent);
    Grid<Real> pca2(parent);
    Grid<Real> pca3(parent);
        
    // setup matrix and boundaries
    MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);
    SetOpenBound (A0, Ai, Aj, Ak, vel, loOpenBound, upOpenBound);
    
    if (ghostAccuracy > 0) {
        if (!phi) errMsg("solve_pressure: if ghostAccuracy>0, need to specify levelset phi=xxx");
        ApplyGhostFluid (flags, A0, *phi, ghostAccuracy);
    }
    
    // compute divergence and init right hand side
    MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr);
    
    if (!outflow.empty())
        SetOutflow (rhs, loOutflow, upOutflow, outflowHeight);
    
    if (enforceCompatibility)
        rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
    
    // CG
    const int maxIter = (int)(cgMaxIterFac * flags.getSize().max());
    GridCgInterface *gcg;
    if (vel.is3D())
        gcg = new GridCg<ApplyMatrix>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
    else
        gcg = new GridCg<ApplyMatrix2D>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
    
    gcg->setAccuracy( cgAccuracy ); 
    gcg->setUseResNorm( useResNorm );

    // optional preconditioning
    gcg->setPreconditioner( (GridCgInterface::PreconditionType)precondition, &pca0, &pca1, &pca2, &pca3);

    for (int iter=0; iter<maxIter; iter++) {
        if (!gcg->iterate()) iter=maxIter;
    } 
    debMsg("FluidSolver::solvePressure iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
    delete gcg;
    
    if(ghostAccuracy<=0.) {
        // ghost fluid off, normal correction
        CorrectVelocity (flags, vel, pressure );
    } else {        
        CorrectVelGhostFluid (flags, vel, pressure);
    }    
}PyObject* _plugin_solvePressure (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "solvePressure"); PyObject *_retval = NULL; { ArgLocker _lock; MACGrid& vel = *_args.get< MACGrid* > (0,"vel", &_lock); Grid<Real>& pressure = *_args.get< Grid<Real>* > (1,"pressure", &_lock); FlagGrid& flags = *_args.get< FlagGrid* > (2,"flags", &_lock); Grid<Real>* phi = _args.getOpt< Grid<Real>* > (3,"phi", 0, &_lock); Grid<Real>* perCellCorr = _args.getOpt< Grid<Real>* > (4,"perCellCorr", 0, &_lock); Real ghostAccuracy = _args.getOpt< Real > (5,"ghostAccuracy", 0, &_lock); Real cgMaxIterFac = _args.getOpt< Real > (6,"cgMaxIterFac", 1.5, &_lock); Real cgAccuracy = _args.getOpt< Real > (7,"cgAccuracy", 1e-3, &_lock); string openBound = _args.getOpt< string > (8,"openBound", "", &_lock); string outflow = _args.getOpt< string > (9,"outflow", "", &_lock); int outflowHeight = _args.getOpt< int > (10,"outflowHeight", 1, &_lock); int precondition = _args.getOpt< int > (11,"precondition", 0, &_lock); bool enforceCompatibility = _args.getOpt< bool > (12,"enforceCompatibility", false, &_lock); bool useResNorm = _args.getOpt< bool > (13,"useResNorm", true, &_lock); _retval = getPyNone();solvePressure(vel, pressure, flags, phi, perCellCorr, ghostAccuracy, cgMaxIterFac, cgAccuracy, openBound, outflow, outflowHeight, precondition, enforceCompatibility, useResNorm, parent);_args.check(); } pbFinalizePlugin(parent,"solvePressure"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("solvePressure",e.what()); return 0; } } 

} // end namespace



