




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/grid.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Grid representation
 *
 ******************************************************************************/

#include "grid.h"
#include "levelset.h"
#include "kernel.h"
#include <limits>
#include <sstream>
#include <cstring>
#include "fileio.h"

using namespace std;
namespace Manta {

//******************************************************************************
// GridBase members

GridBase::GridBase (FluidSolver* parent) 
    : PbClass(parent), mType(TypeNone)
{
    checkParent();
    m3D = getParent()->is3D();
}

//******************************************************************************
// Grid<T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>() { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>() { return GridBase::TypeInt; }
template<> inline GridBase::GridType typeList<Vec3>() { return GridBase::TypeVec3; }

template<class T>
Grid<T>::Grid(FluidSolver* parent, bool show)
    : GridBase(parent)
{
    mType = typeList<T>();
    mSize = parent->getGridSize();
    mData = parent->getGridPointer<T>();
    
    mStrideZ = parent->is2D() ? 0 : (mSize.x * mSize.y);
    mDx = 1.0 / mSize.max();
    clear();
    setHidden(!show);
}

template<class T>
Grid<T>::Grid(const Grid<T>& a) : GridBase(a.getParent()) {
    mSize = a.mSize;
    mType = a.mType;
    mStrideZ = a.mStrideZ;
    mDx = a.mDx;
    FluidSolver *gp = a.getParent();
    mData = gp->getGridPointer<T>();
    memcpy(mData, a.mData, sizeof(T) * a.mSize.x * a.mSize.y * a.mSize.z);
}

template<class T>
Grid<T>::~Grid() {
    mParent->freeGridPointer<T>(mData);
}

template<class T>
void Grid<T>::clear() {
    memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z);    
}

template<class T>
void Grid<T>::swap(Grid<T>& other) {
    if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ())
        errMsg("Grid::swap(): Grid dimensions mismatch.");
    
    T* dswap = other.mData;
    other.mData = mData;
    mData = dswap;
}

template<class T>
void Grid<T>::load(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".raw")
        readGridRaw(name, this);
    else if (ext == ".uni")
        readGridUni(name, this);
    else
        errMsg("file '" + name +"' filetype not supported");
}

template<class T>
void Grid<T>::save(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".raw")
        writeGridRaw(name, this);
    else if (ext == ".uni")
        writeGridUni(name, this);
    else
        errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<T> operators

//! Kernel: Compute min value of Real grid

struct CompMinReal : public KernelBase {  CompMinReal (Grid<Real >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), minVal(std::numeric_limits<Real>::max()) { run(); } CompMinReal (const CompMinReal& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), minVal(std::numeric_limits<Real>::max()) {} inline void op(int idx,Grid<Real >& val, Real& minVal)  {
    if (val[idx] < minVal)
        minVal = val[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, minVal); } operator Real  () { return minVal; } const Real& getRet() const { return minVal; } Grid<Real >& getArg0() { return m_val; } typedef Grid<Real > type0; FluidSolver* parent; Grid<Real >& m_val; Real minVal;  }; 

//! Kernel: Compute max value of Real grid

struct CompMaxReal : public KernelBase {  CompMaxReal (Grid<Real >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), maxVal(-std::numeric_limits<Real>::max()) { run(); } CompMaxReal (const CompMaxReal& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), maxVal(-std::numeric_limits<Real>::max()) {} inline void op(int idx,Grid<Real >& val, Real& maxVal)  {
    if (val[idx] > maxVal)
        maxVal = val[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, maxVal); } operator Real  () { return maxVal; } const Real& getRet() const { return maxVal; } Grid<Real >& getArg0() { return m_val; } typedef Grid<Real > type0; FluidSolver* parent; Grid<Real >& m_val; Real maxVal;  }; 

//! Kernel: Compute min value of int grid

struct CompMinInt : public KernelBase {  CompMinInt (Grid<int >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), minVal(std::numeric_limits<int>::max()) { run(); } CompMinInt (const CompMinInt& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), minVal(std::numeric_limits<int>::max()) {} inline void op(int idx,Grid<int >& val, int& minVal)  {
    if (val[idx] < minVal)
        minVal = val[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, minVal); } operator int  () { return minVal; } const int& getRet() const { return minVal; } Grid<int >& getArg0() { return m_val; } typedef Grid<int > type0; FluidSolver* parent; Grid<int >& m_val; int minVal;  }; 

//! Kernel: Compute max value of int grid

struct CompMaxInt : public KernelBase {  CompMaxInt (Grid<int >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), maxVal(-std::numeric_limits<int>::min()) { run(); } CompMaxInt (const CompMaxInt& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), maxVal(-std::numeric_limits<int>::min()) {} inline void op(int idx,Grid<int >& val, int& maxVal)  {
    if (val[idx] > maxVal)
        maxVal = val[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, maxVal); } operator int  () { return maxVal; } const int& getRet() const { return maxVal; } Grid<int >& getArg0() { return m_val; } typedef Grid<int > type0; FluidSolver* parent; Grid<int >& m_val; int maxVal;  }; 

//! Kernel: Compute min norm of vec grid

struct CompMinVec : public KernelBase {  CompMinVec (Grid<Vec3 >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), minVal(std::numeric_limits<Real>::max()) { run(); } CompMinVec (const CompMinVec& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), minVal(std::numeric_limits<Real>::max()) {} inline void op(int idx,Grid<Vec3 >& val, Real& minVal)  {
    const Real s = normSquare(val[idx]);
    if (s < minVal)
        minVal = s;
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, minVal); } operator Real  () { return minVal; } const Real& getRet() const { return minVal; } Grid<Vec3 >& getArg0() { return m_val; } typedef Grid<Vec3 > type0; FluidSolver* parent; Grid<Vec3 >& m_val; Real minVal;  }; 

//! Kernel: Compute max norm of vec grid

struct CompMaxVec : public KernelBase {  CompMaxVec (Grid<Vec3 >& _val) : KernelBase(&_val, 0), parent((_val).getParent()), m_val(_val), maxVal(0) { run(); } CompMaxVec (const CompMaxVec& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_val(o.m_val), maxVal(0) {} inline void op(int idx,Grid<Vec3 >& val, Real& maxVal)  {
    const Real s = normSquare(val[idx]);
    if (s > maxVal)
        maxVal = s;
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_val, maxVal); } operator Real  () { return maxVal; } const Real& getRet() const { return maxVal; } Grid<Vec3 >& getArg0() { return m_val; } typedef Grid<Vec3 > type0; FluidSolver* parent; Grid<Vec3 >& m_val; Real maxVal;  }; 


template<class T> Grid<T>& Grid<T>::safeDivide (const Grid<T>& a) {
    gridSafeDiv<T> (*this, a);
    return *this;
}
template<class T> Grid<T>& Grid<T>::operator= (const Grid<T>& a) {
    assertMsg (a.mSize.x == mSize.x && a.mSize.y == mSize.y && a.mSize.z == mSize.z, "different grid resolutions");
    memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
    mType = a.mType; // copy type marker
    return *this;
}
template<class T> Grid<T>& Grid<T>::operator= (const T& a) {
    FOR_IDX(*this) { mData[idx] = a; }
    return *this;
}
template<class T> void Grid<T>::scaledAdd(const Grid<T>& a, const T& factor) {
    gridScaleAdd<T> (*this, a, factor);
}
template<> Real Grid<Real>::getMaxValue() {
    return CompMaxReal (*this);
}
template<> Real Grid<Real>::getMinValue() {
    return CompMinReal (*this);
}
template<> Real Grid<Real>::getMaxAbsValue() {
    Real amin = CompMinReal (*this);
    Real amax = CompMaxReal (*this);
    return max( fabs(amin), fabs(amax));
}
template<> Real Grid<Vec3>::getMaxValue() {
    return sqrt(CompMaxVec (*this));
}
template<> Real Grid<Vec3>::getMinValue() { 
    return sqrt(CompMinVec (*this));
}
template<> Real Grid<Vec3>::getMaxAbsValue() {
    return sqrt(CompMaxVec (*this));
}
template<> Real Grid<int>::getMaxValue() {
    return (Real) CompMaxInt (*this);
}
template<> Real Grid<int>::getMinValue() {
    return (Real) CompMinInt (*this);
}
template<> Real Grid<int>::getMaxAbsValue() {
    int amin = CompMinInt (*this);
    int amax = CompMaxInt (*this);
    return max( fabs((Real)amin), fabs((Real)amax));
}
template<class T> void Grid<T>::add(const Grid<T>& a, const Grid<T>& b) {
    gridAdd2<T>(*this, a, b);
}

//******************************************************************************
// Specialization classes

void FlagGrid::initDomain(int boundaryWidth) {
    FOR_IDX(*this)
        mData[idx] = TypeEmpty;
    initBoundaries(boundaryWidth);
}

void FlagGrid::initBoundaries(int boundaryWidth) {
    const int w = boundaryWidth;
    FOR_IJK(*this) {
        bool bnd = (i<=w || i>=mSize.x-1-w || j<=w || j>=mSize.y-1-w || (is3D() && (k<=w || k>=mSize.z-1-w)));
        if (bnd) 
            mData[index(i,j,k)] = TypeObstacle;
    }
}

void FlagGrid::updateFromLevelset(LevelsetGrid& levelset) {
    FOR_IDX(*this) {
        if (!isObstacle(idx)) {
            const Real phi = levelset[idx];
            if (phi <= levelset.invalidTimeValue()) continue;
            
            mData[idx] &= ~(TypeEmpty | TypeFluid); // clear empty/fluid flags
            mData[idx] |= (phi <= 0) ? TypeFluid : TypeEmpty; // set resepctive flag
        }
    }
}   

void FlagGrid::fillGrid(int type) {
    FOR_IDX(*this) {
        if ((mData[idx] & TypeObstacle)==0)
            mData[idx] = (mData[idx] & ~(TypeEmpty | TypeFluid)) | type;
    }
}

// explicit instantiation
template class Grid<int>;
template class Grid<Real>;
template class Grid<Vec3>;

} //namespace

