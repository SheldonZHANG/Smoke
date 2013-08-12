




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/advection.cpp"
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
#include "grid.h"
#include "kernel.h"

using namespace std;

namespace Manta { 

//! Semi-Lagrange interpolation kernel


template <class T> struct SemiLagrange : public KernelBase {  SemiLagrange (FlagGrid& _flags, MACGrid& _vel, Grid<T >& _dst, Grid<T >& _src, Real _dt, bool _isLevelset) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_dst(_dst), m_src(_src), m_dt(_dt), m_isLevelset(_isLevelset) { run(); } SemiLagrange (const SemiLagrange& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_dst(o.m_dst), m_src(o.m_src), m_dt(o.m_dt), m_isLevelset(o.m_isLevelset) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Grid<T >& dst, Grid<T >& src, Real dt, bool isLevelset)  {
    cout<<"sheldon's comment:enter SemiLagrange"<<endl;
    if (flags.isObstacle(i,j,k)) {
        dst(i,j,k) = 0;
        return;
    }
    if (!isLevelset && !flags.isFluid(i,j,k) && !flags.isFluid(i-1,j,k) && 
        !flags.isFluid(i,j-1,k) && !flags.isFluid(i,j,k-1)) {
        dst(i,j,k) = src(i,j,k);
        return;
    }
    
    // SL traceback
    Vec3 pos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getCentered(i,j,k) * dt;
    dst(i,j,k) = src.getInterpolated(pos);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_dst, m_src, m_dt, m_isLevelset); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Grid<T >& m_dst; Grid<T >& m_src; Real m_dt; bool m_isLevelset;  }; 

//! Semi-Lagrange interpolation kernel for MAC grids


struct SemiLagrangeMAC : public KernelBase {  SemiLagrangeMAC (FlagGrid& _flags, MACGrid& _vel, MACGrid& _dst, MACGrid& _src, Real _dt) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_dst(_dst), m_src(_src), m_dt(_dt) { run(); } SemiLagrangeMAC (const SemiLagrangeMAC& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_dst(o.m_dst), m_src(o.m_src), m_dt(o.m_dt) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, MACGrid& dst, MACGrid& src, Real dt)  {
    cout<<"sheldon's comment enter SemiLagrangeMAC"<<endl;
    if (flags.isObstacle(i,j,k)) {
        dst(i,j,k) = 0;
        return;
    }
    if (!flags.isFluid(i,j,k) && !flags.isFluid(i-1,j,k) && 
        !flags.isFluid(i,j-1,k) && !flags.isFluid(i,j,k-1)) {
        dst(i,j,k) = src(i,j,k);
        return;
    }
    
    // get currect velocity at MAC position
    // no need to shift xpos etc. as lookup field is also shifted
    Vec3 xpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACX(i,j,k) * dt;
    Real vx = src.getInterpolatedComponent<0>(xpos);
    Vec3 ypos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACY(i,j,k) * dt;
    Real vy = src.getInterpolatedComponent<1>(ypos);
    Vec3 zpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACZ(i,j,k) * dt;
    Real vz = src.getInterpolatedComponent<2>(zpos);
    
    dst(i,j,k) = Vec3(vx,vy,vz);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_dst, m_src, m_dt); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; MACGrid& m_dst; MACGrid& m_src; Real m_dt;  }; 

//! Kernel: Correct based on forward and backward SL steps



template <class T> struct MacCormackCorrect : public KernelBase {  MacCormackCorrect (FlagGrid& _flags, Grid<T >& _dst, Grid<T >& _old, Grid<T >& _fwd, Grid<T >& _bwd, Real _strength, bool _isLevelSet) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_dst(_dst), m_old(_old), m_fwd(_fwd), m_bwd(_bwd), m_strength(_strength), m_isLevelSet(_isLevelSet) { run(); } MacCormackCorrect (const MacCormackCorrect& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_dst(o.m_dst), m_old(o.m_old), m_fwd(o.m_fwd), m_bwd(o.m_bwd), m_strength(o.m_strength), m_isLevelSet(o.m_isLevelSet) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<T >& dst, Grid<T >& old, Grid<T >& fwd, Grid<T >& bwd, Real strength, bool isLevelSet)  {
    cout<<"sheldon's comment enter MacCormackCorrect"<<endl;
    const int idx = flags.index(i,j,k);
    
    if (!flags.isFluid(idx) && !flags.isFluid(i-1,j,k) && !flags.isFluid(i,j-1,k) && !flags.isFluid(i,j,k-1)) {
        dst[idx] = isLevelSet ? fwd[idx] : (T)0.0;
        return;
    }
    
    dst[idx] = fwd[idx] + 0.5*(old[idx] - bwd[idx]);
    
    // interpolate between SL and MC
    if (strength < 1.0)
        dst[idx] = (1.0-strength) * fwd[idx] + strength * dst[idx];
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_dst, m_old, m_fwd, m_bwd, m_strength, m_isLevelSet); } FluidSolver* parent; FlagGrid& m_flags; Grid<T >& m_dst; Grid<T >& m_old; Grid<T >& m_fwd; Grid<T >& m_bwd; Real m_strength; bool m_isLevelSet;  }; 

// Helper to collect min/max in a template
template<class T> inline void getMinMax(T& minv, T& maxv, const T& val) {
    if (val < minv) minv = val;
    if (val > maxv) maxv = val;
}
template<> inline void getMinMax<Vec3>(Vec3& minv, Vec3& maxv, const Vec3& val) {
    getMinMax(minv.x, maxv.x, val.x);
    getMinMax(minv.y, maxv.y, val.y);
    getMinMax(minv.z, maxv.z, val.z);
}

//! Kernel: Clamp obtained value to min/max in source area, and reset values that point out of grid or into boundaries


template <class T> struct MacCormackClamp : public KernelBase {  MacCormackClamp (FlagGrid& _flags, MACGrid& _vel, Grid<T >& _dst, Grid<T >& _orig, Grid<T >& _fwd, Real _dt) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_dst(_dst), m_orig(_orig), m_fwd(_fwd), m_dt(_dt) { run(); } MacCormackClamp (const MacCormackClamp& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_dst(o.m_dst), m_orig(o.m_orig), m_fwd(o.m_fwd), m_dt(o.m_dt) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Grid<T >& dst, Grid<T >& orig, Grid<T >& fwd, Real dt)  {
    cout<<"sheldon's comment enter MacCormackClamp"<<endl;
    if (flags.isObstacle(i,j,k))
        return;
    if ((!flags.isFluid(i,j,k) && !flags.isFluid(i-1,j,k) && !flags.isFluid(i,j-1,k) && !flags.isFluid(i,j,k-1))) {
        dst(i,j,k) = fwd(i,j,k);
        return;
    }
    
    // lookup forward/backward
    Vec3i posFwd = toVec3i( Vec3(i,j,k) - vel.getCentered(i,j,k) * dt );
    Vec3i posBwd = toVec3i( Vec3(i,j,k) + vel.getCentered(i,j,k) * dt );
    
    // clamp forward lookup to grid
    const int i0 = clamp(posFwd.x, 0, flags.getSizeX()-2);
    const int j0 = clamp(posFwd.y, 0, flags.getSizeY()-2);
    const int k0 = clamp(posFwd.z, 0, flags.getSizeZ()-2);
    const int i1 = i0+1, j1 = j0+1, k1=k0+1;
    
    if (orig.isInBounds(Vec3i(i0,j0,k0),1)) {           
        // find min/max around fwd pos
        T minv = orig(i0,j0,k0), maxv = minv;
        getMinMax(minv, maxv, orig(i1,j0,k0));
        getMinMax(minv, maxv, orig(i0,j1,k0));
        getMinMax(minv, maxv, orig(i1,j1,k0));
        getMinMax(minv, maxv, orig(i0,j0,k1));
        getMinMax(minv, maxv, orig(i1,j0,k1));
        getMinMax(minv, maxv, orig(i0,j1,k1));
        getMinMax(minv, maxv, orig(i1,j1,k1));
        
        // write clamped value
        dst(i,j,k) = clamp(dst(i,j,k), minv, maxv);
    }
    
    // test if lookups point out of grid or into obstacle
    if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
        posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
        posFwd.x >= flags.getSizeX()-1 || posFwd.y >= flags.getSizeY()-1 || posFwd.z >= flags.getSizeZ()-1 ||
        posBwd.x >= flags.getSizeX()-1 || posBwd.y >= flags.getSizeY()-1 || posBwd.z >= flags.getSizeZ()-1 ||
        flags.isObstacle(posFwd) || flags.isObstacle(posBwd) ) 
    {        
        dst(i,j,k) = fwd(i,j,k);
    }
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_dst, m_orig, m_fwd, m_dt); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Grid<T >& m_dst; Grid<T >& m_orig; Grid<T >& m_fwd; Real m_dt;  }; 

//! Helper function for clamping MAC grids
template<int c> 
inline Real doClampComponent(const Vec3i& upperClamp, MACGrid& orig, Real dst, const Vec3& posFwd) {
    cout<<"sheldon's comment enter doClampComponent"<<endl;
    // clamp forward lookup to grid
    const int i0 = clamp((int)posFwd.x, 0, upperClamp.x);
    const int j0 = clamp((int)posFwd.y, 0, upperClamp.y);
    const int k0 = clamp((int)posFwd.z, 0, upperClamp.z);
    const int i1 = i0+1, j1 = j0+1, k1=k0+1;
    if (!orig.isInBounds(Vec3i(i0,j0,k0),1)) 
        return dst;
    
    // find min/max around fwd pos
    Real minv = orig(i0,j0,k0)[c], maxv = minv;
    getMinMax(minv, maxv, orig(i1,j0,k0)[c]);
    getMinMax(minv, maxv, orig(i0,j1,k0)[c]);
    getMinMax(minv, maxv, orig(i1,j1,k0)[c]);
    getMinMax(minv, maxv, orig(i0,j0,k1)[c]);
    getMinMax(minv, maxv, orig(i1,j0,k1)[c]);
    getMinMax(minv, maxv, orig(i0,j1,k1)[c]);
    getMinMax(minv, maxv, orig(i1,j1,k1)[c]);
    
    return clamp(dst, minv, maxv);    
}

//! Kernel: Clamp obtained value to min/max in source area, and reset values that point out of grid or into boundaries. 
//! Specialized version for MAC grids


struct MacCormackClampMAC : public KernelBase {  MacCormackClampMAC (FlagGrid& _flags, MACGrid& _vel, MACGrid& _dst, MACGrid& _orig, MACGrid& _fwd, Real _dt) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_dst(_dst), m_orig(_orig), m_fwd(_fwd), m_dt(_dt) { run(); } MacCormackClampMAC (const MacCormackClampMAC& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_dst(o.m_dst), m_orig(o.m_orig), m_fwd(o.m_fwd), m_dt(o.m_dt) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, MACGrid& dst, MACGrid& orig, MACGrid& fwd, Real dt)  {
    cout<<"sheldon's comment enter MacCormackClampMac"<<endl;
    if (flags.isObstacle(i,j,k))
        return;
    if ((!flags.isFluid(i,j,k) && !flags.isFluid(i-1,j,k) && !flags.isFluid(i,j-1,k) && !flags.isFluid(i,j,k-1))) {
        dst(i,j,k) = fwd(i,j,k);
        return;
    }
    
    Vec3 pos(i,j,k);
    Vec3 dval = dst(i,j,k);
    Vec3i upperClamp = flags.getSize() - 1;
    
    // clamp individual components
    dval.x = doClampComponent<0>(upperClamp, orig, dval.x, pos - vel.getAtMACX(i,j,k) * dt);
    dval.y = doClampComponent<1>(upperClamp, orig, dval.y, pos - vel.getAtMACY(i,j,k) * dt);
    dval.z = doClampComponent<2>(upperClamp, orig, dval.z, pos - vel.getAtMACZ(i,j,k) * dt);
    
    // get total fwd lookup
    Vec3i posFwd = toVec3i( Vec3(i,j,k) - vel.getCentered(i,j,k) * dt );
    Vec3i posBwd = toVec3i( Vec3(i,j,k) + vel.getCentered(i,j,k) * dt );
    
    // test if lookups point out of grid or into obstacle
    if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
        posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
        posFwd.x >= upperClamp.x || posFwd.y >= upperClamp.y || posFwd.z >= upperClamp.z ||
        posBwd.x >= upperClamp.x || posBwd.y >= upperClamp.y || posBwd.z >= upperClamp.z ||
        flags.isObstacle(posFwd) || flags.isObstacle(posBwd) ) 
    {        
        dval = fwd(i,j,k);
    }
    
    // writeback
    dst(i,j,k) = dval;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_dst, m_orig, m_fwd, m_dt); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; MACGrid& m_dst; MACGrid& m_orig; MACGrid& m_fwd; Real m_dt;  }; 

//! template function for performing SL advection
template<class GridType> 
void fnAdvectSemiLagrange(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, GridType& orig, int order, Real strength) {
    cout<<"sheldon's comment enter fnAdvectSemiLagrange"<<endl;
    typedef typename GridType::BASETYPE T;
    
    Real dt = parent->getDt();
    bool levelset = orig.getType() & GridBase::TypeLevelset;
    
    // forward step
    GridType fwd(parent);
    SemiLagrange<T> (flags, vel, fwd, orig, dt, levelset);
    
    if (order == 1) {
        orig.swap(fwd);
    }
    else if (order == 2) { // MacCormack
        GridType bwd(parent);
        GridType newGrid(parent);
    
        // bwd <- backwards step
        SemiLagrange<T> (flags, vel, bwd, fwd, -dt, levelset);
        
        // newGrid <- compute correction
        MacCormackCorrect<T> (flags, newGrid, orig, fwd, bwd, strength, levelset);
        
        // clamp values
        MacCormackClamp<T> (flags, vel, newGrid, orig, fwd, dt);
        
        orig.swap(newGrid);
    }
}

//! template function for performing SL advection: specialized version for MAC grids
template<> 
void fnAdvectSemiLagrange<MACGrid>(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, MACGrid& orig, int order, Real strength) {
    cout<<"sheldon's comment fnAdvectSemiLagrange"<<endl;
    Real dt = parent->getDt();
    
    // forward step
    MACGrid fwd(parent);    
    SemiLagrangeMAC (flags, vel, fwd, orig, dt);
    
    if (order == 1) {
        orig.swap(fwd);
    }
    else if (order == 2) { // MacCormack
        MACGrid bwd(parent);
        MACGrid newGrid(parent);
        
        // bwd <- backwards step
        SemiLagrangeMAC (flags, vel, bwd, fwd, -dt);
        
        // newGrid <- compute correction
        MacCormackCorrect<Vec3> (flags, newGrid, orig, fwd, bwd, strength, false);
        
        // clamp values
        MacCormackClampMAC (flags, vel, newGrid, orig, fwd, dt);
        
        orig.swap(newGrid);
    }
}

//! Perform semi-lagrangian advection of target Real- or Vec3 grid


void advectSemiLagrange(FlagGrid* flags,MACGrid* vel,GridBase* grid,int order = 1,Real strength = 0.5, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {    
    cout<<"sheldon's comment advectSemiLagrange"<<endl;
    assertMsg(order==1 || order==2, "AdvectSemiLagrange: Only order 1 (regular SL) and 2 (MacCormack) supported");
    
    // determine type of grid    
    if (grid->getType() & GridBase::TypeReal) {
        fnAdvectSemiLagrange< Grid<Real> >(parent, *flags, *vel, *((Grid<Real>*) grid), order, strength);
    }
    else if (grid->getType() & GridBase::TypeMAC) {    
        fnAdvectSemiLagrange< MACGrid >(parent, *flags, *vel, *((MACGrid*) grid), order, strength);
    }
    else if (grid->getType() & GridBase::TypeVec3) {    
        fnAdvectSemiLagrange< Grid<Vec3> >(parent, *flags, *vel, *((Grid<Vec3>*) grid), order, strength);
    }
    else
        errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");    
}PyObject* _plugin_advectSemiLagrange (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "advectSemiLagrange"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid* flags = _args.get< FlagGrid* > (0,"flags", &_lock); MACGrid* vel = _args.get< MACGrid* > (1,"vel", &_lock); GridBase* grid = _args.get< GridBase* > (2,"grid", &_lock); int order = _args.getOpt< int > (3,"order", 1, &_lock); Real strength = _args.getOpt< Real > (4,"strength", 0.5, &_lock); _retval = getPyNone();advectSemiLagrange(flags, vel, grid, order, strength, parent);_args.check(); } pbFinalizePlugin(parent,"advectSemiLagrange"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("advectSemiLagrange",e.what()); return 0; } } 

} // end namespace DDF 



