




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/extforces.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Set boundary conditions, gravity
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"

using namespace std;

namespace Manta { 

//! add Forces between fl/fl and fl/em cells
struct KnAddForceField : public KernelBase {  KnAddForceField (FlagGrid& _flags, MACGrid& _vel, Grid<Vec3 >& _force) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_force(_force) { run(); } KnAddForceField (const KnAddForceField& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_force(o.m_force) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Grid<Vec3 >& force)  {
    cout<<"sheldon's comment KnAddForce"<<endl;
    bool curFluid = flags.isFluid(i,j,k);
    bool curEmpty = flags.isEmpty(i,j,k);
    if (!curFluid && !curEmpty) return;
    
    if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
        vel(i,j,k).x += 0.5*(force(i-1,j,k).x + force(i,j,k).x);
    if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
        vel(i,j,k).y += 0.5*(force(i,j-1,k).y + force(i,j,k).y);
    if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
        vel(i,j,k).z += 0.5*(force(i,j,k-1).z + force(i,j,k).z);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_force); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Grid<Vec3 >& m_force;  }; 

//! add Forces between fl/fl and fl/em cells
struct KnAddForce : public KernelBase {  KnAddForce (FlagGrid& _flags, MACGrid& _vel, Vec3 _force) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel), m_force(_force) { run(); } KnAddForce (const KnAddForce& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel), m_force(o.m_force) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel, Vec3 force)  {
    cout<<"sheldon comment KnAddForce"<<endl;
    bool curFluid = flags.isFluid(i,j,k);
    bool curEmpty = flags.isEmpty(i,j,k);
    if (!curFluid && !curEmpty) return;
    
    if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
        vel(i,j,k).x += force.x;
    if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
        vel(i,j,k).y += force.y;
    if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
        vel(i,j,k).z += force.z;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel, m_force); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel; Vec3 m_force;  }; 

//! add gravity forces to all fluid cells
void addGravity(FlagGrid& flags,MACGrid& vel,Vec3 gravity, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {    
    cout<<"sheldon's comment addGravity"<<endl;
    Vec3 f = gravity * parent->getDt() / flags.getDx();
    KnAddForce(flags, vel, f);
}PyObject* _plugin_addGravity (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addGravity"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *_args.get< MACGrid* > (1,"vel", &_lock); Vec3 gravity = _args.get< Vec3 > (2,"gravity", &_lock); _retval = getPyNone();addGravity(flags, vel, gravity, parent);_args.check(); } pbFinalizePlugin(parent,"addGravity"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("addGravity",e.what()); return 0; } } 

//! add Buoyancy force based on smoke density
struct KnAddBuoyancy : public KernelBase {  KnAddBuoyancy (FlagGrid& _flags, Grid<Real >& _density, MACGrid& _vel, Vec3 _strength) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_density(_density), m_vel(_vel), m_strength(_strength) { run(); } KnAddBuoyancy (const KnAddBuoyancy& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_density(o.m_density), m_vel(o.m_vel), m_strength(o.m_strength) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<Real >& density, MACGrid& vel, Vec3 strength)  {    
    cout<<"sheldon comment KnAddBuoyancy"<<endl;
    if (!flags.isFluid(i,j,k)) return;
    if (flags.isFluid(i-1,j,k))
        vel(i,j,k).x += (0.5 * strength.x) * (density(i,j,k)+density(i-1,j,k));
    if (flags.isFluid(i,j-1,k))
        vel(i,j,k).y += (0.5 * strength.y) * (density(i,j,k)+density(i,j-1,k));
    if (vel.is3D() && flags.isFluid(i,j,k-1))
        vel(i,j,k).z += (0.5 * strength.z) * (density(i,j,k)+density(i,j,k-1));    
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_density, m_vel, m_strength); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_density; MACGrid& m_vel; Vec3 m_strength;  }; 

//! add Buoyancy force based on smoke density
void addBuoyancy(FlagGrid& flags,Grid<Real>& density,MACGrid& vel,Vec3 gravity, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    cout<<"sheldon comment addBuoyancy"<<endl;
    Vec3 f = - gravity * parent->getDt() / parent->getDx();
    KnAddBuoyancy(flags,density, vel, f);
}PyObject* _plugin_addBuoyancy (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "addBuoyancy"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); Grid<Real>& density = *_args.get< Grid<Real>* > (1,"density", &_lock); MACGrid& vel = *_args.get< MACGrid* > (2,"vel", &_lock); Vec3 gravity = _args.get< Vec3 > (3,"gravity", &_lock); _retval = getPyNone();addBuoyancy(flags, density, vel, gravity, parent);_args.check(); } pbFinalizePlugin(parent,"addBuoyancy"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("addBuoyancy",e.what()); return 0; } } 

        
//! set no-stick wall boundary condition between ob/fl and ob/ob cells
struct KnSetWallBcs : public KernelBase {  KnSetWallBcs (FlagGrid& _flags, MACGrid& _vel) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel) { run(); } KnSetWallBcs (const KnSetWallBcs& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel)  {
    bool curFluid = flags.isFluid(i,j,k);
    bool curObstacle = flags.isObstacle(i,j,k);
    if (!curFluid && !curObstacle) return;
    
    // we use i>0 instead of bnd=1 to check outer wall
    if (i>0 && (flags.isObstacle(i-1,j,k) || (curObstacle && flags.isFluid(i-1,j,k))))
        vel(i,j,k).x = 0;
    if (j>0 && (flags.isObstacle(i,j-1,k) || (curObstacle && flags.isFluid(i,j-1,k))))
        vel(i,j,k).y = 0;
    if (vel.is2D() || (k>0 && (flags.isObstacle(i,j,k-1) || (curObstacle && flags.isFluid(i,j,k-1)))))
        vel(i,j,k).z = 0;
		
	if (curFluid) {
		if ((i>0 && flags.isStick(i-1,j,k)) || (i<flags.getSizeX()-1 && flags.isStick(i+1,j,k)))
			vel(i,j,k).y = vel(i,j,k).z = 0;
		if ((j>0 && flags.isStick(i,j-1,k)) || (j<flags.getSizeY()-1 && flags.isStick(i,j+1,k)))
			vel(i,j,k).x = vel(i,j,k).z = 0;
		if (vel.is3D() && ((k>0 && flags.isStick(i,j,k-1)) || (k<flags.getSizeZ()-1 && flags.isStick(i,j,k+1))))
			vel(i,j,k).x = vel(i,j,k).y = 0;
	}
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k, m_flags, m_vel); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel;  }; 

//! set no-stick boundary condition on walls
void setWallBcs(FlagGrid& flags,MACGrid& vel, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    KnSetWallBcs(flags, vel);
}PyObject* _plugin_setWallBcs (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setWallBcs"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *_args.get< MACGrid* > (1,"vel", &_lock); _retval = getPyNone();setWallBcs(flags, vel, parent);_args.check(); } pbFinalizePlugin(parent,"setWallBcs"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("setWallBcs",e.what()); return 0; } }  

//! set boundary conditions at empty cells
struct KnSetLiquidBcs : public KernelBase {  KnSetLiquidBcs (FlagGrid& _flags, MACGrid& _vel) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_vel(_vel) { run(); } KnSetLiquidBcs (const KnSetLiquidBcs& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_vel(o.m_vel) {} inline void op(int i,int j,int k,FlagGrid& flags, MACGrid& vel)  {
    if (!flags.isFluid(i,j,k)) return;
    
    // init empty cells from fluid
    if (flags.isEmpty(i+1,j,k)) 
        vel(i+1,j,k).x = vel(i,j,k).x;
    if (flags.isEmpty(i,j+1,k)) 
        vel(i,j+1,k).y = vel(i,j,k).y;
    if (flags.isEmpty(i,j,k+1)) 
        vel(i,j,k+1).z = vel(i,j,k).z;
    
    // "left" sides of fluid - fluid cells neighboring
    // empty cells along neg. dir, get velocities from fluid 
    if (flags.isEmpty(i-1,j,k) && flags.isFluid(i+1,j,k)) 
        vel(i,j,k).x = vel(i+1,j,k).x;
    if (flags.isEmpty(i,j-1,k) && flags.isFluid(i,j+1,k)) 
        vel(i,j,k).y = vel(i,j+1,k).y;
    if (flags.isEmpty(i,j,k-1) && flags.isFluid(i,j,k+1)) 
        vel(i,j,k).z = vel(i,j,k+1).z;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_vel); } FluidSolver* parent; FlagGrid& m_flags; MACGrid& m_vel;  }; 

//! set boundary conditions at empty cells
void setLiquidBcs(FlagGrid& flags,MACGrid& vel, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    FOR_IDX(flags) {
        if (flags.isEmpty(idx)) vel[idx]=_0;
    }
    KnSetLiquidBcs(flags, vel);
}PyObject* _plugin_setLiquidBcs (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setLiquidBcs"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *_args.get< MACGrid* > (1,"vel", &_lock); _retval = getPyNone();setLiquidBcs(flags, vel, parent);_args.check(); } pbFinalizePlugin(parent,"setLiquidBcs"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("setLiquidBcs",e.what()); return 0; } }  

//! Kernel: gradient norm operator
struct KnConfForce : public KernelBase {  KnConfForce (Grid<Vec3 >& _force, const Grid<Real >& _grid, const Grid<Vec3 >& _curl, Real _str) : KernelBase(&_force, 1), parent((_force).getParent()), m_force(_force), m_grid(_grid), m_curl(_curl), m_str(_str) { run(); } KnConfForce (const KnConfForce& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_force(o.m_force), m_grid(o.m_grid), m_curl(o.m_curl), m_str(o.m_str) {} inline void op(int i,int j,int k,Grid<Vec3 >& force, const Grid<Real >& grid, const Grid<Vec3 >& curl, Real str)  {
    Vec3 grad = 0.5 * Vec3( grid(i+1,j,k)-grid(i-1,j,k), grid(i,j+1,k)-grid(i,j-1,k), grid(i,j,k+1)-grid(i,j,k-1));
    normalize(grad);
    force(i,j,k) = str*cross(grad, curl(i,j,k));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_force, m_grid, m_curl, m_str); } FluidSolver* parent; Grid<Vec3 >& m_force; const Grid<Real >& m_grid; const Grid<Vec3 >& m_curl; Real m_str;  }; 

void vorticityConfinement(MACGrid& vel,FlagGrid& flags,Real strength, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    assertMsg(vel.is3D(), "Only 3D grids supported so far");
    Grid<Vec3> velCenter(parent), curl(parent), force(parent);
    Grid<Real> norm(parent);
    
    GetCentered(velCenter, vel);
    CurlOp(velCenter, curl);
    GridNorm(norm, curl);
    KnConfForce(force, norm, curl, strength);
    KnAddForceField(flags, vel, force);
}PyObject* _plugin_vorticityConfinement (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "vorticityConfinement"); PyObject *_retval = NULL; { ArgLocker _lock; MACGrid& vel = *_args.get< MACGrid* > (0,"vel", &_lock); FlagGrid& flags = *_args.get< FlagGrid* > (1,"flags", &_lock); Real strength = _args.get< Real > (2,"strength", &_lock); _retval = getPyNone();vorticityConfinement(vel, flags, strength, parent);_args.check(); } pbFinalizePlugin(parent,"vorticityConfinement"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("vorticityConfinement",e.what()); return 0; } } 


} // namespace


