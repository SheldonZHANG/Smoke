




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/levelset.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Levelset
 *
 ******************************************************************************/

#include "levelset.h"
#include "fastmarch.h"
#include "kernel.h"
#include "mcubes.h"
#include "mesh.h"

using namespace std;
namespace Manta {

//************************************************************************
// Helper functions and kernels for marching

static const int FlagInited = FastMarch<FmHeapEntryOut, +1>::FlagInited;

// neighbor lookup vectors
static const Vec3i neighbors[6] = { Vec3i(-1,0,0), Vec3i(1,0,0), Vec3i(0,-1,0), Vec3i(0,1,0), Vec3i(0,0,-1), Vec3i(0,0,1) };
    

struct InitFmIn : public KernelBase {  InitFmIn (FlagGrid& _flags, Grid<int >& _fmFlags, LevelsetGrid& _phi, bool _ignoreWalls) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_fmFlags(_fmFlags), m_phi(_phi), m_ignoreWalls(_ignoreWalls) { run(); } InitFmIn (const InitFmIn& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_fmFlags(o.m_fmFlags), m_phi(o.m_phi), m_ignoreWalls(o.m_ignoreWalls) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<int >& fmFlags, LevelsetGrid& phi, bool ignoreWalls)  {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (v>=0 && (!ignoreWalls || !flags.isObstacle(idx)))
        fmFlags[idx] = FlagInited;
    else
        fmFlags[idx] = 0;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_fmFlags, m_phi, m_ignoreWalls); } FluidSolver* parent; FlagGrid& m_flags; Grid<int >& m_fmFlags; LevelsetGrid& m_phi; bool m_ignoreWalls;  }; 


struct InitFmOut : public KernelBase {  InitFmOut (FlagGrid& _flags, Grid<int >& _fmFlags, LevelsetGrid& _phi, bool _ignoreWalls) : KernelBase(&_flags, 1), parent((_flags).getParent()), m_flags(_flags), m_fmFlags(_fmFlags), m_phi(_phi), m_ignoreWalls(_ignoreWalls) { run(); } InitFmOut (const InitFmOut& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_flags(o.m_flags), m_fmFlags(o.m_fmFlags), m_phi(o.m_phi), m_ignoreWalls(o.m_ignoreWalls) {} inline void op(int i,int j,int k,FlagGrid& flags, Grid<int >& fmFlags, LevelsetGrid& phi, bool ignoreWalls)  {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (ignoreWalls) {
        fmFlags[idx] = (v<0) ? FlagInited : 0;
        if (flags.isObstacle(idx)) {
            fmFlags[idx] = 0;
            phi[idx] = 0;
        }
    }
    else
        fmFlags[idx] = (v<0) ? FlagInited : 0;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_flags, m_fmFlags, m_phi, m_ignoreWalls); } FluidSolver* parent; FlagGrid& m_flags; Grid<int >& m_fmFlags; LevelsetGrid& m_phi; bool m_ignoreWalls;  }; 


struct SetUninitialized : public KernelBase {  SetUninitialized (Grid<int >& _fmFlags, LevelsetGrid& _phi, const Real _val) : KernelBase(&_fmFlags, 1), parent((_fmFlags).getParent()), m_fmFlags(_fmFlags), m_phi(_phi), m_val(_val) { run(); } SetUninitialized (const SetUninitialized& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_fmFlags(o.m_fmFlags), m_phi(o.m_phi), m_val(o.m_val) {} inline void op(int i,int j,int k,Grid<int >& fmFlags, LevelsetGrid& phi, const Real val)  {
    if (fmFlags(i,j,k) != FlagInited)
        phi(i,j,k) = val;
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_fmFlags, m_phi, m_val); } FluidSolver* parent; Grid<int >& m_fmFlags; LevelsetGrid& m_phi; const Real m_val;  }; 

template<bool inward>
inline bool isAtInterface(Grid<int>& fmFlags, LevelsetGrid& phi, const Vec3i& p) {
    // check for interface
    for (int nb=0; nb<6; nb++) {
        const Vec3i pn(p + neighbors[nb]);
        if (!fmFlags.isInBounds(pn)) continue;
        
        if (fmFlags(pn) != FlagInited) continue;
        if ((inward && phi(pn) >= 0) || 
            (!inward && phi(pn) < 0)) return true;
    }
    return false;
}

//************************************************************************
// Levelset class def

LevelsetGrid::LevelsetGrid(FluidSolver* parent, bool show) 
    : Grid<Real>(parent, show) 
{ 
    mType = (GridType)(TypeLevelset | TypeReal);    
}    

Real LevelsetGrid::invalidTimeValue() {
    return FastMarch<FmHeapEntryOut, 1>::InvalidTime();
}

//! Kernel: perform levelset union
struct KnJoin : public KernelBase {  KnJoin (Grid<Real >& _a, const Grid<Real >& _b) : KernelBase(&_a, 0), parent((_a).getParent()), m_a(_a), m_b(_b) { run(); } KnJoin (const KnJoin& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_a(o.m_a), m_b(o.m_b) {} inline void op(int idx,Grid<Real >& a, const Grid<Real >& b)  {
    a[idx] = min(a[idx], b[idx]);
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_a, m_b); } FluidSolver* parent; Grid<Real >& m_a; const Grid<Real >& m_b;  }; 

void LevelsetGrid::join(const LevelsetGrid& o) {
    KnJoin(*this, o);
}

void LevelsetGrid::reinitMarching(FlagGrid& flags, Real maxTime, MACGrid* velTransport, bool ignoreWalls, bool correctOuterLayer)
{
    assertMsg(is3D(), "Only 3D grids supported so far");
    
    Grid<int> fmFlags(mParent);
    LevelsetGrid& phi = *this;
    
    FastMarch<FmHeapEntryOut, +1> marchOut(flags, fmFlags, phi, maxTime, velTransport);
    FastMarch<FmHeapEntryIn, -1> marchIn(flags, fmFlags, phi, maxTime, NULL);
    
    // march inside
    InitFmIn (flags, fmFlags, phi, ignoreWalls);
    
    FOR_IJK_BND(flags, 1) {
        if (fmFlags(i,j,k) == FlagInited) continue;
        if (flags.isObstacle(i,j,k)) continue;
        const Vec3i p(i,j,k);
                
        if(isAtInterface<true>(fmFlags, phi, p)) {
            // set value
            fmFlags(p) = FlagInited;
            
            // add neighbors that are not at the interface
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                if (flags.isObstacle(pn)) continue;
                
                // check neighbors of neighbor
                if (phi(pn) < 0 && !isAtInterface<true>(fmFlags, phi, pn)) {
                    marchIn.addToList(pn, p); 
                }
            }            
        }
    }
    marchIn.performMarching();     
    
    // set un initialized regions
    SetUninitialized (fmFlags, phi, -maxTime); 
    
    // done with inwards marching, now march out...    
    InitFmOut (flags, fmFlags, phi, ignoreWalls);
    
    // by default, correctOuterLayer is on
    if (correctOuterLayer) {    
        // normal version, inwards march is done, now add all outside values (0..2] to list
        // note, this might move the interface a bit! but keeps a nice signed distance field...        
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            const Vec3i p(i,j,k);
            
            // check nbs
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                
                if (fmFlags(pn) != FlagInited) continue;
                if (flags.isObstacle(pn)) continue;
                
                const Real nbPhi = phi(pn);
                
                // only add nodes near interface, not e.g. outer boundary vs. invalid region                
                if (nbPhi < 0 && nbPhi >= -2)
                    marchOut.addToList(p, pn); 
            }
        }         
    } else {
        // alternative version, keep interface, do not distort outer cells
        // add all ouside values, but not those at the IF layer
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            
            // only look at ouside values
            const Vec3i p(i,j,k);
            if (phi(p) < 0) continue;
            
            if (isAtInterface<false>(fmFlags, phi, p)) {
                // now add all non, interface neighbors
                fmFlags(p) = FlagInited;
                
                // add neighbors that are not at the interface
                for (int nb=0; nb<6; nb++) {
                    const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                    if (flags.isObstacle(pn)) continue;
                
                    // check neighbors of neighbor
                    if (phi(pn) > 0 && !isAtInterface<false>(fmFlags, phi, pn))
                        marchOut.addToList(pn, p);
                }            
            }
        }
    }    
    marchOut.performMarching();
    
    // set un initialized regions
    SetUninitialized (fmFlags, phi, +maxTime);    
    
}


void LevelsetGrid::initFromFlags(FlagGrid& flags, bool ignoreWalls) {
    FOR_IDX(*this) {
        if (flags.isFluid(idx) || (ignoreWalls && flags.isObstacle(idx)))
            mData[idx] = -0.5;
        else
            mData[idx] = 0.5;
    }
    reinitMarching(flags, getSize().max(), NULL, ignoreWalls, true);
}
    

// helper function
inline Vec3 getNormal(const Grid<Real>& data, int i, int j, int k) {
    if (i > data.getSizeX()-2) i= data.getSizeX()-2;
    if (j > data.getSizeY()-2) j= data.getSizeY()-2;
    if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
    return Vec3 (data(i-1,j  ,k  ) - data(i+1,j  ,k  ),
                 data(i  ,j-1,k  ) - data(i  ,j+1,k  ),
                 data(i  ,j  ,k-1) - data(i  ,j  ,k+1));
}


void LevelsetGrid::createMesh(Mesh& mesh) {
    assertMsg(is3D(), "Only 3D grids supported so far");
    
    mesh.clear();
        
    const Real invalidTime = invalidTimeValue();
    const Real isoValue = 1e-4;
    
    // create some temp grids
    Grid<int> edgeVX(mParent);
    Grid<int> edgeVY(mParent);
    Grid<int> edgeVZ(mParent);
    
    for(int i=1; i<mSize.x-1; i++)
    for(int j=1; j<mSize.y-1; j++)
    for(int k=1; k<mSize.z-1; k++) {
         Real value[8] = { get(i,j,k),   get(i+1,j,k),   get(i+1,j+1,k),   get(i,j+1,k),
                                get(i,j,k+1), get(i+1,j,k+1), get(i+1,j+1,k+1), get(i,j+1,k+1) };
        
        // build lookup index, check for invalid times
        bool skip = false;
        int cubeIdx = 0;
        for (int l=0;l<8;l++) {
            value[l] *= -1;
            if (-value[l] <= invalidTime)
                skip = true;
            if (value[l] < isoValue) 
                cubeIdx |= 1<<l;
        }
        if (skip || (mcEdgeTable[cubeIdx] == 0)) continue;
        
        // where to look up if this point already exists
        int triIndices[12];
        int *eVert[12] = { &edgeVX(i,j,k),   &edgeVY(i+1,j,k),   &edgeVX(i,j+1,k),   &edgeVY(i,j,k), 
                           &edgeVX(i,j,k+1), &edgeVY(i+1,j,k+1), &edgeVX(i,j+1,k+1), &edgeVY(i,j,k+1), 
                           &edgeVZ(i,j,k),   &edgeVZ(i+1,j,k),   &edgeVZ(i+1,j+1,k), &edgeVZ(i,j+1,k) };
        
        const Vec3 pos[9] = { Vec3(i,j,k),   Vec3(i+1,j,k),   Vec3(i+1,j+1,k),   Vec3(i,j+1,k),
                        Vec3(i,j,k+1), Vec3(i+1,j,k+1), Vec3(i+1,j+1,k+1), Vec3(i,j+1,k+1) };
        
        for (int e=0; e<12; e++) {
            if (mcEdgeTable[cubeIdx] & (1<<e)) {
                // vertex already calculated ?
                if (*eVert[e] == 0) {
                    // interpolate edge
                    const int e1 = mcEdges[e*2  ];
                    const int e2 = mcEdges[e*2+1];
                    const Vec3 p1 = pos[ e1  ];    // scalar field pos 1
                    const Vec3 p2 = pos[ e2  ];    // scalar field pos 2
                    const float valp1  = value[ e1  ];  // scalar field val 1
                    const float valp2  = value[ e2  ];  // scalar field val 2
                    const float mu = (isoValue - valp1) / (valp2 - valp1);

                    // init isolevel vertex
                    Node vertex;
                    vertex.pos = p1 + (p2-p1)*mu;
                    vertex.normal = getNormalized( 
                                        getNormal( *this, i+cubieOffsetX[e1], j+cubieOffsetY[e1], k+cubieOffsetZ[e1]) * (1.0-mu) +
                                        getNormal( *this, i+cubieOffsetX[e2], j+cubieOffsetY[e2], k+cubieOffsetZ[e2]) * (    mu)) ;
                    
                    triIndices[e] = mesh.addNode(vertex) + 1;
                    
                    // store vertex 
                    *eVert[e] = triIndices[e];
                } else {
                    // retrieve  from vert array
                    triIndices[e] = *eVert[e];
                }
            }
        }
        
        // Create the triangles... 
        for(int e=0; mcTriTable[cubeIdx][e]!=-1; e+=3) {
            mesh.addTri( Triangle( triIndices[ mcTriTable[cubeIdx][e+0]] - 1,
                                        triIndices[ mcTriTable[cubeIdx][e+1]] - 1,
                                        triIndices[ mcTriTable[cubeIdx][e+2]] - 1));
        }
    }
    
    //mesh.rebuildCorners();
    //mesh.rebuildLookup();
}


} //namespace

