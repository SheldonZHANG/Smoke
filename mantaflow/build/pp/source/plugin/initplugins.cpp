




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/initplugins.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "noisefield.h"

using namespace std;

namespace Manta {
    
//! Apply noise to grid


struct KnApplyNoise : public KernelBase {  KnApplyNoise (FlagGrid& _flags, Grid<Real >& _dens, WaveletNoiseField& _noise, Grid<Real >& _sdf, Real _scale, Real _sigma) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_dens(_dens), m_noise(_noise), m_sdf(_sdf), m_scale(_scale), m_sigma(_sigma) { run(); } KnApplyNoise (const KnApplyNoise& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_dens(o.m_dens), m_noise(o.m_noise), m_sdf(o.m_sdf), m_scale(o.m_scale), m_sigma(o.m_sigma) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<Real >& dens, WaveletNoiseField& noise, Grid<Real >& sdf, Real scale, Real sigma)  {
    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
    Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
    
    Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
    if (dens(i,j,k) < target)
        dens(i,j,k) = target;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k, m_flags, m_dens, m_noise, m_sdf, m_scale, m_sigma); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_dens; WaveletNoiseField& m_noise; Grid<Real >& m_sdf; Real m_scale; Real m_sigma;  }; 

//! Init noise-moduled density inside shape

void densityInflow(FlagGrid& flags,Grid<Real>& density,WaveletNoiseField& noise,Shape* shape,Real scale=1.0,Real sigma=0, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Grid<Real> sdf = shape->computeLevelset();    
    KnApplyNoise(flags, density, noise, sdf, scale, sigma);
}PyObject* _plugin_densityInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "densityInflow"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); Grid<Real>& density = *_args.get< Grid<Real>* > (1,"density", &_lock); WaveletNoiseField& noise = *_args.get< WaveletNoiseField* > (2,"noise", &_lock); Shape* shape = _args.get< Shape* > (3,"shape", &_lock); Real scale = _args.getOpt< Real > (4,"scale", 1.0, &_lock); Real sigma = _args.getOpt< Real > (5,"sigma", 0, &_lock); _retval = getPyNone();densityInflow(flags, density, noise, shape, scale, sigma, parent);_args.check(); } pbFinalizePlugin(parent,"densityInflow"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("densityInflow",e.what()); return 0; } } 

    
} // namespace

