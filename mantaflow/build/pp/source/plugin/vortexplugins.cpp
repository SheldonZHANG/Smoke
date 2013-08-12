




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/vortexplugins.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for using vortex sheet meshes 
 *
 ******************************************************************************/
 
#include <iostream>
#include "vortexsheet.h"
#include "vortexpart.h"
#include "shapes.h"
#include "commonkernels.h"
#include "conjugategrad.h"
#include "randomstream.h"
#include "levelset.h"

using namespace std;

namespace Manta {
    
//! Mark area of mesh inside shape as fixed nodes. 
//! Remove all other fixed nodes if 'exclusive' is set

void markAsFixed(Mesh& mesh,Shape* shape,bool exclusive=true, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    cout<<"Enter markAsFixed"<<endl;
    for (int i=0; i<mesh.numNodes(); i++) {
        if (shape->isInside(mesh.nodes(i).pos))
            mesh.nodes(i).flags |= Mesh::NfFixed;
        else if (exclusive)
            mesh.nodes(i).flags &= ~Mesh::NfFixed;
    }
    cout<<"Ext markAsFixed"<<endl;
}PyObject* _plugin_markAsFixed (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "markAsFixed"); PyObject *_retval = NULL; { ArgLocker _lock; Mesh& mesh = *_args.get< Mesh* > (0,"mesh", &_lock); Shape* shape = _args.get< Shape* > (1,"shape", &_lock); bool exclusive = _args.getOpt< bool > (2,"exclusive", true, &_lock); _retval = getPyNone();markAsFixed(mesh, shape, exclusive, parent);_args.check(); } pbFinalizePlugin(parent,"markAsFixed"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("markAsFixed",e.what()); return 0; } } 

//! Adapt texture coordinates of mesh inside shape
//! to obtain an effective inflow effect

void texcoordInflow(VortexSheetMesh& mesh,Shape* shape,MACGrid& vel, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    static Vec3 t0 = Vec3::Zero;
    
    // get mean velocity
    int cnt=0;
    Vec3 meanV(_0);
    FOR_IJK(vel) {
        if (shape->isInsideGrid(i,j,k)) {
            cnt++;
            meanV += vel.getCentered(i,j,k);
        }
    }
    meanV /= (Real) cnt;
    t0 -= parent->getDt() * meanV;
    mesh.setReferenceTexOffset(t0);

    // apply mean velocity
    for (int i=0; i<mesh.numNodes(); i++) {
        if (shape->isInside(mesh.nodes(i).pos)) {
            Vec3 tc = mesh.nodes(i).pos + t0;
            mesh.tex1(i) = tc;
            mesh.tex2(i) = tc;
        }
    }
}PyObject* _plugin_texcoordInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "texcoordInflow"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); Shape* shape = _args.get< Shape* > (1,"shape", &_lock); MACGrid& vel = *_args.get< MACGrid* > (2,"vel", &_lock); _retval = getPyNone();texcoordInflow(mesh, shape, vel, parent);_args.check(); } pbFinalizePlugin(parent,"texcoordInflow"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("texcoordInflow",e.what()); return 0; } } ;

//! Init smoke density values of the mesh surface inside source shape

void meshSmokeInflow(VortexSheetMesh& mesh,Shape* shape,Real amount, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    for (int t=0; t<mesh.numTris(); t++) {
        if (shape->isInside(mesh.getFaceCenter(t)))
            mesh.sheet(t).smokeAmount = amount;
    }    
}PyObject* _plugin_meshSmokeInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "meshSmokeInflow"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); Shape* shape = _args.get< Shape* > (1,"shape", &_lock); Real amount = _args.get< Real > (2,"amount", &_lock); _retval = getPyNone();meshSmokeInflow(mesh, shape, amount, parent);_args.check(); } pbFinalizePlugin(parent,"meshSmokeInflow"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("meshSmokeInflow",e.what()); return 0; } } 


struct KnAcceleration : public KernelBase {  KnAcceleration (MACGrid& _a, const MACGrid& _v1, const MACGrid& _v0, const Real _idt) : KernelBase(&_a, 0), parent((_a).getParent()), m_a(_a), m_v1(_v1), m_v0(_v0), m_idt(_idt) { run(); } KnAcceleration (const KnAcceleration& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_a(o.m_a), m_v1(o.m_v1), m_v0(o.m_v0), m_idt(o.m_idt) {} inline void op(int idx,MACGrid& a, const MACGrid& v1, const MACGrid& v0, const Real idt)  { 
    a[idx] = (v1[idx]-v0[idx])*idt; 
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_a, m_v1, m_v0, m_idt); } FluidSolver* parent; MACGrid& m_a; const MACGrid& m_v1; const MACGrid& m_v0; const Real m_idt;  }; 

//! Add vorticity to vortex sheets based on buoyancy



void vorticitySource(VortexSheetMesh& mesh,Vec3 gravity,MACGrid* vel=NULL,MACGrid* velOld=NULL,Real scale = 0.1,Real maxAmount = 0,Real mult = 1.0, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Real dt = parent->getDt();
    Real dx = parent->getDx();
    MACGrid acceleration(parent);
    if (vel)
        KnAcceleration(acceleration, *vel, *velOld, 1.0/dt);
    const Real A= -1.0;
    Real maxV = 0, meanV = 0;
    
    for (int t=0; t<mesh.numTris(); t++) {
        Vec3 fn = mesh.getFaceNormal(t);
        Vec3 source;
        if (vel) {
            Vec3 a = acceleration.getInterpolated(mesh.getFaceCenter(t));        
            source = A*cross(fn, a-gravity) * scale;
        } else {
            source = A*cross(fn, -gravity) * scale;
        }
        
        if (mesh.isTriangleFixed(t)) source = 0;
    
        mesh.sheet(t).vorticity *= mult;
        mesh.sheet(t).vorticity += dt * source / dx;
        // upper limit
        Real v = norm(mesh.sheet(t).vorticity);
        if (maxAmount>0 && v > maxAmount)
            mesh.sheet(t).vorticity *= maxAmount/v;
        
        //stats
        if (v > maxV) maxV = v;
        meanV += v;
    }
    
    cout << "vorticity: max " << maxV << " / mean " << meanV/mesh.numTris() << endl;
}PyObject* _plugin_vorticitySource (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "vorticitySource"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); Vec3 gravity = _args.get< Vec3 > (1,"gravity", &_lock); MACGrid* vel = _args.getOpt< MACGrid* > (2,"vel", NULL, &_lock); MACGrid* velOld = _args.getOpt< MACGrid* > (3,"velOld", NULL, &_lock); Real scale = _args.getOpt< Real > (4,"scale", 0.1, &_lock); Real maxAmount = _args.getOpt< Real > (5,"maxAmount", 0, &_lock); Real mult = _args.getOpt< Real > (6,"mult", 1.0, &_lock); _retval = getPyNone();vorticitySource(mesh, gravity, vel, velOld, scale, maxAmount, mult, parent);_args.check(); } pbFinalizePlugin(parent,"vorticitySource"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("vorticitySource",e.what()); return 0; } } 


void smoothVorticity(VortexSheetMesh& mesh,int iter=1,Real sigma=0.2,Real alpha=0.8, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    const Real mult = -0.5 / sigma / sigma;
    
    // pre-calculate positions and weights
    vector<Vec3> vort(mesh.numTris()), pos(mesh.numTris());    
    vector<Real> weights(3*mesh.numTris());
    vector<int> index(3*mesh.numTris());
    for(int i=0; i<mesh.numTris(); i++) {
        pos[i] = mesh.getFaceCenter(i);
        mesh.sheet(i).vorticitySmoothed = mesh.sheet(i).vorticity;
    }
    for(int i=0; i<mesh.numTris(); i++) {
        for (int c=0; c<3; c++) {
            int oc = mesh.corners(i,c).opposite;
            if (oc>=0) {
                int t = mesh.corners(oc).tri;
                weights[3*i+c] = exp(normSquare(pos[t]-pos[i])*mult);
                index[3*i+c] = t;
            }
            else {
                weights[3*i+c] = 0;
                index[3*i+c] = 0;
            }
        }        
    }
        
    for (int it=0; it<iter; ++it) {
        // first, preload
        for(int i=0; i<mesh.numTris(); i++) vort[i] = mesh.sheet(i).vorticitySmoothed;
            
        for(int i=0,idx=0; i<mesh.numTris(); i++) {            
            // loop over adjacent tris
            Real sum=1.0f;
            Vec3 v=vort[i];
            for (int c=0;c<3;c++,idx++) {
                Real w = weights[index[idx]];
                v += w*vort[index[idx]];
                sum += w;
            }
            mesh.sheet(i).vorticitySmoothed = v/sum;
        }
    }
    for(int i=0; i<mesh.numTris(); i++) mesh.sheet(i).vorticitySmoothed *= alpha;
}PyObject* _plugin_smoothVorticity (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "smoothVorticity"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); int iter = _args.getOpt< int > (1,"iter", 1, &_lock); Real sigma = _args.getOpt< Real > (2,"sigma", 0.2, &_lock); Real alpha = _args.getOpt< Real > (3,"alpha", 0.8, &_lock); _retval = getPyNone();smoothVorticity(mesh, iter, sigma, alpha, parent);_args.check(); } pbFinalizePlugin(parent,"smoothVorticity"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("smoothVorticity",e.what()); return 0; } } 

//! Seed Vortex Particles inside shape with K41 characteristics
void VPseedK41(VortexParticleSystem& system,Shape* shape,Real strength=0,Real sigma0=0.2,Real sigma1=1.0,Real probability=1.0,Real N=3.0, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Grid<Real> temp(parent);
    const Real dt = parent->getDt();
    static RandomStream rand(3489572);
    Real s0 = pow( (Real)sigma0, (Real)(-N+1.0) );
    Real s1 = pow( (Real)sigma1, (Real)(-N+1.0) );
    
    FOR_IJK(temp) {
        if (shape->isInsideGrid(i,j,k)) {
            if (rand.getReal() < probability*dt) {
                Real p = rand.getReal();
                Real sigma = pow( (1.0-p)*s0 + p*s1, 1./(-N+1.0) );
                Vec3 randDir (rand.getReal(), rand.getReal(), rand.getReal());
                Vec3 posUpd (i+rand.getReal(), j+rand.getReal(), k+rand.getReal());
                normalize(randDir);
                Vec3 vorticity = randDir * strength * pow( (Real)sigma, (Real)(-10./6.+N/2.0) );
                system.add(VortexParticleData(posUpd, vorticity, sigma)); 
            }
        }
    }    
}PyObject* _plugin_VPseedK41 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "VPseedK41"); PyObject *_retval = NULL; { ArgLocker _lock; VortexParticleSystem& system = *_args.get< VortexParticleSystem* > (0,"system", &_lock); Shape* shape = _args.get< Shape* > (1,"shape", &_lock); Real strength = _args.getOpt< Real > (2,"strength", 0, &_lock); Real sigma0 = _args.getOpt< Real > (3,"sigma0", 0.2, &_lock); Real sigma1 = _args.getOpt< Real > (4,"sigma1", 1.0, &_lock); Real probability = _args.getOpt< Real > (5,"probability", 1.0, &_lock); Real N = _args.getOpt< Real > (6,"N", 3.0, &_lock); _retval = getPyNone();VPseedK41(system, shape, strength, sigma0, sigma1, probability, N, parent);_args.check(); } pbFinalizePlugin(parent,"VPseedK41"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("VPseedK41",e.what()); return 0; } } 
        
//! Vortex-in-cell integration

void VICintegration(VortexSheetMesh& mesh,Real sigma,Grid<Vec3>& vel,FlagGrid& flags,Grid<Vec3>* vorticity=NULL,Real cgMaxIterFac=1.5,Real cgAccuracy=1e-3,Real scale = 0.01,int precondition=0, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    
    MuTime t0;
    const Real fac = 16.0; // experimental factor to balance out regularization
    
    // if no vort grid is given, use a temporary one
    Grid<Vec3> vortTemp(parent);    
    Grid<Vec3>& vort = (vorticity) ? (*vorticity) : (vortTemp);
    vort.clear();
    
    // map vorticity to grid using Peskin kernel
    int sgi = ceil(sigma);
    Real pkfac=M_PI/sigma;
    const int numTris = mesh.numTris();
    for (int t=0; t<numTris; t++) {
        Vec3 pos = mesh.getFaceCenter(t);
        Vec3 v = mesh.sheet(t).vorticity * mesh.getFaceArea(t) * fac;
                    
        // inner kernel
        // first, summate 
        Real sum=0;
        for (int i=-sgi; i<sgi; i++) {
            if (pos.x+i < 0 || (int)pos.x+i >= vort.getSizeX()) continue;
            for (int j=-sgi; j<sgi; j++) {
                if (pos.y+j < 0 || (int)pos.y+j >= vort.getSizeY()) continue;            
                for (int k=-sgi; k<sgi; k++) {
                    if (pos.z+k < 0 || (int)pos.z+k >= vort.getSizeZ()) continue;                                
                    Vec3i cell(pos.x+i, pos.y+j, pos.z+k);
                    if (!flags.isFluid(cell)) continue;
                    Vec3 d = pos - Vec3(i+0.5+floor(pos.x), j+0.5+floor(pos.y), k+0.5+floor(pos.z));
                    Real dl = norm(d);
                    if (dl > sigma) continue;
                    // precalc Peskin kernel
                    sum += 1.0 + cos(dl * pkfac);
                }
            }
        }
        // then, apply normalized kernel
        Real wnorm = 1.0/sum;
        for (int i=-sgi; i<sgi; i++) {
            if (pos.x+i < 0 || (int)pos.x+i >= vort.getSizeX()) continue;
            for (int j=-sgi; j<sgi; j++) {
                if (pos.y+j < 0 || (int)pos.y+j >= vort.getSizeY()) continue;            
                for (int k=-sgi; k<sgi; k++) {
                    if (pos.z+k < 0 || (int)pos.z+k >= vort.getSizeZ()) continue;                                
                    Vec3i cell(pos.x+i, pos.y+j, pos.z+k);  
                    if (!flags.isFluid(cell)) continue;                    
                    Vec3 d = pos - Vec3(i+0.5+floor(pos.x), j+0.5+floor(pos.y), k+0.5+floor(pos.z));
                    Real dl = norm(d);
                    if (dl > sigma) continue;
                    Real w = (1.0 + cos(dl * pkfac))*wnorm;
                    vort(cell) += v * w;
                }
            }
        }
    }
    
    // Prepare grids for poisson solve
    Grid<Vec3> vortexCurl(parent);
    Grid<Real> rhs(parent);
    Grid<Real> solution(parent);
    Grid<Real> residual(parent);
    Grid<Real> search(parent);
    Grid<Real> temp1(parent);
    Grid<Real> A0(parent);
    Grid<Real> Ai(parent);
    Grid<Real> Aj(parent);
    Grid<Real> Ak(parent);
    Grid<Real> pca0(parent);
    Grid<Real> pca1(parent);
    Grid<Real> pca2(parent);
    Grid<Real> pca3(parent);
    
    MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak);    
    CurlOp(vort, vortexCurl);    
    
    // Solve vector poisson equation
    for (int c=0; c<3; c++) {
        // construct rhs    
        if (vel.getType() & GridBase::TypeMAC)
            GetShiftedComponent(vortexCurl, rhs, c);
        else
            GetComponent(vortexCurl, rhs, c);
                
        // prepare CG solver
        const int maxIter = (int)(cgMaxIterFac * vel.getSize().max());    
        GridCgInterface *gcg = new GridCg<ApplyMatrix>(solution, rhs, residual, search, flags, temp1, &A0, &Ai, &Aj, &Ak );
        gcg->setAccuracy(cgAccuracy); 
        gcg->setUseResNorm(true);
        gcg->setPreconditioner( (GridCgInterface::PreconditionType)precondition, &pca0, &pca1, &pca2, &pca3); 
        
        // iterations
        for (int iter=0; iter<maxIter; iter++) {
            if (!gcg->iterate()) iter=maxIter;
        } 
        debMsg("VICintegration CG iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
        delete gcg;
    
        // copy back
        solution *= scale;
        SetComponent(vel, solution, c);
    }
}PyObject* _plugin_VICintegration (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "VICintegration"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); Real sigma = _args.get< Real > (1,"sigma", &_lock); Grid<Vec3>& vel = *_args.get< Grid<Vec3>* > (2,"vel", &_lock); FlagGrid& flags = *_args.get< FlagGrid* > (3,"flags", &_lock); Grid<Vec3>* vorticity = _args.getOpt< Grid<Vec3>* > (4,"vorticity", NULL, &_lock); Real cgMaxIterFac = _args.getOpt< Real > (5,"cgMaxIterFac", 1.5, &_lock); Real cgAccuracy = _args.getOpt< Real > (6,"cgAccuracy", 1e-3, &_lock); Real scale = _args.getOpt< Real > (7,"scale", 0.01, &_lock); int precondition = _args.getOpt< int > (8,"precondition", 0, &_lock); _retval = getPyNone();VICintegration(mesh, sigma, vel, flags, vorticity, cgMaxIterFac, cgAccuracy, scale, precondition, parent);_args.check(); } pbFinalizePlugin(parent,"VICintegration"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("VICintegration",e.what()); return 0; } } 

//! Obtain density field from levelset with linear gradient of size sigma over the interface
void densityFromLevelset(LevelsetGrid& phi,Grid<Real>& density,Real value=1.0,Real sigma=1.0, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    FOR_IJK(phi) {
        // remove boundary
        if (i<2 || j<2 || k<2 || i>=phi.getSizeX()-2 || j>=phi.getSizeY()-2 || k>=phi.getSizeZ()-2)
            density(i,j,k) = 0;
        else if (phi(i,j,k) < -sigma)
            density(i,j,k) = value;
        else if (phi(i,j,k) > sigma)
            density(i,j,k) = 0;
        else
            density(i,j,k) = clamp((Real)(0.5*value/sigma*(1.0-phi(i,j,k))), _0, value);
    }    
}PyObject* _plugin_densityFromLevelset (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "densityFromLevelset"); PyObject *_retval = NULL; { ArgLocker _lock; LevelsetGrid& phi = *_args.get< LevelsetGrid* > (0,"phi", &_lock); Grid<Real>& density = *_args.get< Grid<Real>* > (1,"density", &_lock); Real value = _args.getOpt< Real > (2,"value", 1.0, &_lock); Real sigma = _args.getOpt< Real > (3,"sigma", 1.0, &_lock); _retval = getPyNone();densityFromLevelset(phi, density, value, sigma, parent);_args.check(); } pbFinalizePlugin(parent,"densityFromLevelset"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("densityFromLevelset",e.what()); return 0; } } 

} // namespace


