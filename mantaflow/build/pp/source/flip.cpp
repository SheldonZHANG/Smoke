




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/flip.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * FLIP (fluid implicit particles)
 *
 ******************************************************************************/

#include "flip.h"

using namespace std;
namespace Manta {

//! Set velocities from grid with given PIC/FLIP mixture

struct CopyVelocitiesFromGrid : public ParticleKernelBase {  CopyVelocitiesFromGrid (FlipSystem& _p, FlagGrid& _flags, MACGrid& _vel, MACGrid& _oldVel, Real _flipRatio) : ParticleKernelBase((_p).size()), parent((_p).getParent()), m_p(_p), m_flags(_flags), m_vel(_vel), m_oldVel(_oldVel), m_flipRatio(_flipRatio) { run(); } CopyVelocitiesFromGrid (const CopyVelocitiesFromGrid& o) : ParticleKernelBase(o.size), parent(o.parent), m_p(o.m_p), m_flags(o.m_flags), m_vel(o.m_vel), m_oldVel(o.m_oldVel), m_flipRatio(o.m_flipRatio) {} inline void op(int i,FlipSystem& p, FlagGrid& flags, MACGrid& vel, MACGrid& oldVel, Real flipRatio)  {
    unusedParameter(flags);
    
    if (!p.isActive(i)) return;
    
    /*if (!flags.isFluid(p[i].pos)) {
        p[i].flag |= ParticleBase::PDELETE;
        return;
    }*/
    
    Vec3 v = vel.getInterpolated(p[i].pos);
    Vec3 delta = v - oldVel.getInterpolated(p[i].pos);
    
    p[i].vel = flipRatio * (p[i].vel + delta) + (1.0f - flipRatio) * v;    
} void run() { const int _sz = size; for (int i=0; i < _sz; i++) op(i, m_p, m_flags, m_vel, m_oldVel, m_flipRatio); } FluidSolver* parent; FlipSystem& m_p; FlagGrid& m_flags; MACGrid& m_vel; MACGrid& m_oldVel; Real m_flipRatio;  }; 

//! Set velocities on the grid from the particle system

struct CopyVelocitiesToGrid : public ParticleKernelBase {  CopyVelocitiesToGrid (FlipSystem& _p, FlagGrid& _flags, MACGrid& _vel, Grid<Vec3 >& _tmp) : ParticleKernelBase((_p).size()), parent((_p).getParent()), m_p(_p), m_flags(_flags), m_vel(_vel), m_tmp(_tmp) { run(); } CopyVelocitiesToGrid (const CopyVelocitiesToGrid& o) : ParticleKernelBase(o.size), parent(o.parent), m_p(o.m_p), m_flags(o.m_flags), m_vel(o.m_vel), m_tmp(o.m_tmp) {} inline void op(int i,FlipSystem& p, FlagGrid& flags, MACGrid& vel, Grid<Vec3 >& tmp)  {
    unusedParameter(flags);
    
    if (!p.isActive(i)) return;
    
    vel.setInterpolated(p[i].pos, p[i].vel, &tmp[0]);
} void run() { const int _sz = size; for (int i=0; i < _sz; i++) op(i, m_p, m_flags, m_vel, m_tmp); } FluidSolver* parent; FlipSystem& m_p; FlagGrid& m_flags; MACGrid& m_vel; Grid<Vec3 >& m_tmp;  }; 

void FlipSystem::velocitiesFromGrid(FlagGrid& flags, MACGrid& vel, Real flipRatio) {
    assertMsg(vel.is3D(), "Only 3D grids supported so far");
    CopyVelocitiesFromGrid(*this, flags, vel, mOldVel, flipRatio);
}

void FlipSystem::velocitiesToGrid(FlagGrid& flags, MACGrid& vel) {
    assertMsg(vel.is3D(), "Only 3D grids supported so far");
    
    // interpol -> grid. tmpgrid for counting
    Grid<Vec3> tmp(mParent);
    vel.clear();
    CopyVelocitiesToGrid(*this, flags, vel, tmp);
    vel.safeDivide(tmp);
    
    // store diff
    mOldVel = vel;
}

void FlipSystem::adjustNumber(MACGrid& vel, FlagGrid& flags, int minParticles, int maxParticles) {
    Grid<int> tmp(mParent);
    
    // count particles in cells, and delete excess particles
    for (size_t i=0; i<mData.size(); i++) {
        if (isActive(i)) {
            Vec3i p = toVec3i(mData[i].pos);
            int num = tmp(p);
            
            if (!flags.isFluid(p) || num > maxParticles)
                mData[i].flag |= PDELETE;
            else
                tmp(p) = num+1;
        }
    }
    
    compress();
    
    // seed new particles
    FOR_IJK(tmp) {
        int cnt = tmp(i,j,k);
        if (flags.isFluid(i,j,k) && cnt < minParticles) {
            for (int m=cnt; m < minParticles; m++) { 
                Vec3 rndPos (i + mRand.getReal(), j + mRand.getReal(), k + mRand.getReal());
                add(FlipData(rndPos, vel.getInterpolated(rndPos.x)));                
            }
        }
    }
}

void FlipSystem::markFluidCells(FlagGrid& flags) {
    // kill all fluid cells
    flags.fillGrid(FlagGrid::TypeEmpty);
    
    // mark all particles in flaggrid
    for(size_t i=0;i<mData.size();i++) {
        const Vec3i p = toVec3i(mData[i].pos);
        if (flags.isInBounds(p) && flags.isEmpty(p))
            flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
    }
}

ParticleBase* FlipSystem::clone() {
    FlipSystem* nm = new FlipSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
    return nm;
}



} // namespace


