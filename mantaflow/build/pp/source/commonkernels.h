




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/commonkernels.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Common grid kernels
 *
 ******************************************************************************/

#ifndef _COMMONKERNELS_H
#define _COMMONKERNELS_H

#include "general.h"
#include "kernel.h"
#include "grid.h"

namespace Manta {
   
//! Kernel: Invert real values, if positive and fluid


struct InvertCheckFluid : public KernelBase {  InvertCheckFluid (FlagGrid& _flags, Grid<Real >& _grid) : KernelBase(&_flags, 0), parent((_flags).getParent()), m_flags(_flags), m_grid(_grid) { run(); } InvertCheckFluid (const InvertCheckFluid& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_grid(o.m_grid) {} inline void op(int idx,FlagGrid& flags, Grid<Real >& grid)  {
    if (flags.isFluid(idx) && grid[idx] > 0)
        grid[idx] = 1.0 / grid[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_flags, m_grid); } FluidSolver* parent; FlagGrid& m_flags; Grid<Real >& m_grid;  }; 

//! Kernel: Squared sum over grid

struct GridSumSqr : public KernelBase {  GridSumSqr (Grid<Real >& _grid) : KernelBase(&_grid, 0), parent((_grid).getParent()), m_grid(_grid), sum(0) { run(); } GridSumSqr (const GridSumSqr& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), sum(0) {} inline void op(int idx,Grid<Real >& grid, double& sum)  {
    sum += square((double)grid[idx]);    
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_grid, sum); } operator double  () { return sum; } const double& getRet() const { return sum; } Grid<Real >& getArg0() { return m_grid; } typedef Grid<Real > type0; FluidSolver* parent; Grid<Real >& m_grid; double sum;  }; 

//! Kernel: rotation operator \nabla x v for centered vector fields

struct CurlOp : public KernelBase {  CurlOp (const Grid<Vec3 >& _grid, Grid<Vec3 >& _dst) : KernelBase(&_grid, 1), parent((_grid).getParent()), m_grid(_grid), m_dst(_dst) { run(); } CurlOp (const CurlOp& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), m_dst(o.m_dst) {} inline void op(int i,int j,int k,const Grid<Vec3 >& grid, Grid<Vec3 >& dst)  {
    dst(i,j,k) = Vec3(0.5*((grid(i,j+1,k).z - grid(i,j-1,k).z) - (grid(i,j,k+1).y - grid(i,j,k-1).y)),
                      0.5*((grid(i,j,k+1).x - grid(i,j,k-1).x) - (grid(i+1,j,k).z - grid(i-1,j,k).z)),
                      0.5*((grid(i+1,j,k).y - grid(i-1,j,k).y) - (grid(i,j+1,k).x - grid(i,j-1,k).x)));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_grid, m_dst); } FluidSolver* parent; const Grid<Vec3 >& m_grid; Grid<Vec3 >& m_dst;  }; ;

//! Kernel: divergence operator (from MAC grid)

struct DivergenceOpMAC : public KernelBase {  DivergenceOpMAC (Grid<Real >& _div, const MACGrid& _grid) : KernelBase(&_div, 1), parent((_div).getParent()), m_div(_div), m_grid(_grid) { run(); } DivergenceOpMAC (const DivergenceOpMAC& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_div(o.m_div), m_grid(o.m_grid) {} inline void op(int i,int j,int k,Grid<Real >& div, const MACGrid& grid)  {
    Vec3 del = Vec3(grid(i+1,j,k).x, grid(i,j+1,k).y, grid(i,j,k+1).z) - grid(i,j,k);
    div(i,j,k) = del.x + del.y + del.z;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_div, m_grid); } FluidSolver* parent; Grid<Real >& m_div; const MACGrid& m_grid;  }; 

//! Kernel: gradient operator (create MAC grid)
struct GradientOpMAC : public KernelBase {  GradientOpMAC (MACGrid& _gradient, const Grid<Real >& _grid) : KernelBase(&_gradient, 1), parent((_gradient).getParent()), m_gradient(_gradient), m_grid(_grid) { run(); } GradientOpMAC (const GradientOpMAC& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_gradient(o.m_gradient), m_grid(o.m_grid) {} inline void op(int i,int j,int k,MACGrid& gradient, const Grid<Real >& grid)  {
    gradient(i,j,k) = (Vec3(grid(i,j,k)) - Vec3(grid(i-1,j,k), grid(i,j-1,k), grid(i,j,k-1)));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_gradient, m_grid); } FluidSolver* parent; MACGrid& m_gradient; const Grid<Real >& m_grid;  }; 

//! Kernel: gradient operator 
struct GradientOp : public KernelBase {  GradientOp (Grid<Vec3 >& _gradient, const Grid<Real >& _grid) : KernelBase(&_gradient, 1), parent((_gradient).getParent()), m_gradient(_gradient), m_grid(_grid) { run(); } GradientOp (const GradientOp& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_gradient(o.m_gradient), m_grid(o.m_grid) {} inline void op(int i,int j,int k,Grid<Vec3 >& gradient, const Grid<Real >& grid)  {
    gradient(i,j,k) = 0.5 * Vec3( grid(i+1,j,k)-grid(i-1,j,k), grid(i,j+1,k)-grid(i,j-1,k), grid(i,j,k+1)-grid(i,j,k-1));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_gradient, m_grid); } FluidSolver* parent; Grid<Vec3 >& m_gradient; const Grid<Real >& m_grid;  }; 

//! Kernel: get component at MAC positions
struct GetShiftedComponent : public KernelBase {  GetShiftedComponent (const Grid<Vec3 >& _grid, Grid<Real >& _comp, int _dim) : KernelBase(&_grid, 1), parent((_grid).getParent()), m_grid(_grid), m_comp(_comp), m_dim(_dim) { run(); } GetShiftedComponent (const GetShiftedComponent& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), m_comp(o.m_comp), m_dim(o.m_dim) {} inline void op(int i,int j,int k,const Grid<Vec3 >& grid, Grid<Real >& comp, int dim)  {
    Vec3i ishift(i,j,k);
    ishift[dim]--;
    comp(i,j,k) = 0.5*(grid(i,j,k)[dim] + grid(ishift)[dim]);
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_grid, m_comp, m_dim); } FluidSolver* parent; const Grid<Vec3 >& m_grid; Grid<Real >& m_comp; int m_dim;  }; ;

//! Kernel: get component (not shifted)
struct GetComponent : public KernelBase {  GetComponent (const Grid<Vec3 >& _grid, Grid<Real >& _comp, int _dim) : KernelBase(&_grid, 0), parent((_grid).getParent()), m_grid(_grid), m_comp(_comp), m_dim(_dim) { run(); } GetComponent (const GetComponent& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), m_comp(o.m_comp), m_dim(o.m_dim) {} inline void op(int idx,const Grid<Vec3 >& grid, Grid<Real >& comp, int dim)  {
    comp[idx] = grid[idx][dim];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_grid, m_comp, m_dim); } FluidSolver* parent; const Grid<Vec3 >& m_grid; Grid<Real >& m_comp; int m_dim;  }; ;

//! Kernel: get norm of centered grid
struct GridNorm : public KernelBase {  GridNorm (Grid<Real >& _n, const Grid<Vec3 >& _grid) : KernelBase(&_n, 0), parent((_n).getParent()), m_n(_n), m_grid(_grid) { run(); } GridNorm (const GridNorm& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_n(o.m_n), m_grid(o.m_grid) {} inline void op(int idx,Grid<Real >& n, const Grid<Vec3 >& grid)  {
    n[idx] = norm(grid[idx]);
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_n, m_grid); } FluidSolver* parent; Grid<Real >& m_n; const Grid<Vec3 >& m_grid;  }; ;

//! Kernel: set component (not shifted)
struct SetComponent : public KernelBase {  SetComponent (Grid<Vec3 >& _grid, const Grid<Real >& _comp, int _dim) : KernelBase(&_grid, 0), parent((_grid).getParent()), m_grid(_grid), m_comp(_comp), m_dim(_dim) { run(); } SetComponent (const SetComponent& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), m_comp(o.m_comp), m_dim(o.m_dim) {} inline void op(int idx,Grid<Vec3 >& grid, const Grid<Real >& comp, int dim)  {
    grid[idx][dim] = comp[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_grid, m_comp, m_dim); } FluidSolver* parent; Grid<Vec3 >& m_grid; const Grid<Real >& m_comp; int m_dim;  }; ;

//! Kernel: compute centered velocity field from MAC
struct GetCentered : public KernelBase {  GetCentered (Grid<Vec3 >& _center, const MACGrid& _vel) : KernelBase(&_center, 1), parent((_center).getParent()), m_center(_center), m_vel(_vel) { run(); } GetCentered (const GetCentered& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_center(o.m_center), m_vel(o.m_vel) {} inline void op(int i,int j,int k,Grid<Vec3 >& center, const MACGrid& vel)  {
    center(i,j,k) = 0.5*(vel(i,j,k)+Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, vel(i,j,k+1).z));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_center, m_vel); } FluidSolver* parent; Grid<Vec3 >& m_center; const MACGrid& m_vel;  }; ;

//! Kernel: compute MAC from centered velocity field
struct GetMAC : public KernelBase {  GetMAC (MACGrid& _vel, const Grid<Vec3 >& _center) : KernelBase(&_vel, 1), parent((_vel).getParent()), m_vel(_vel), m_center(_center) { run(); } GetMAC (const GetMAC& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_vel(o.m_vel), m_center(o.m_center) {} inline void op(int i,int j,int k,MACGrid& vel, const Grid<Vec3 >& center)  {
    vel(i,j,k) = 0.5*(center(i,j,k)+Vec3(center(i-1,j,k).x, center(i,j-1,k).y, center(i,j,k-1).z));
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_vel, m_center); } FluidSolver* parent; MACGrid& m_vel; const Grid<Vec3 >& m_center;  }; ;

} // namespace
#endif

