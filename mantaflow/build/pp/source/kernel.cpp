




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/kernel.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Function and macros for defining compution kernels over grids
 *
 ******************************************************************************/

#include "kernel.h"
#include "grid.h"
#include "particle.h"

namespace Manta {

KernelBase::KernelBase(const GridBase* base, int bnd) :    
    maxX (base->getSizeX()-bnd),
    maxY (base->getSizeY()-bnd),
    maxZ (base->is3D() ? (base->getSizeZ()-bnd) : 1),
    minZ (base->is3D() ? bnd : 0),
    X (base->getStrideX()),
    Y (base->getStrideY()),
    Z (base->getStrideZ()),
    maxCells (base->getSizeX() * base->getSizeY() * base->getSizeZ()) {}
    
KernelBase::KernelBase(int _maxX, int _maxY, int _maxZ, int _maxC, int _minZ, int _X, int _Y, int _Z) :
    maxX(_maxX), maxY(_maxY), maxZ(_maxZ), maxCells(_maxC), minZ(_minZ), X(_X), Y(_Y), Z(_Z) {}
    
ParticleKernelBase::ParticleKernelBase(int sz) :
    size(sz) {}
    
    
} // namespace

