




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/plugin/turbulence.cpp"
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
 
#include "grid.h"
#include "commonkernels.h"
#include "vortexsheet.h"

using namespace std;

namespace Manta {

// k-epsilon model constants
const Real keCmu = 0.09;
const Real keC1 = 1.44;
const Real keC2 = 1.92;
const Real keS1 = 1.0;
const Real keS2 = 1.3;

// k-epsilon limiters
const Real keU0 = 1.0;
const Real keImin = 2e-3;
const Real keImax = 1.0;
const Real keNuMin = 1e-3;
const Real keNuMax = 4.0;

//! clamp k and epsilon to limits    

struct KnTurbulenceClamp : public KernelBase {  KnTurbulenceClamp (Grid<Real >& _kgrid, Grid<Real >& _egrid, Real _minK, Real _maxK, Real _minNu, Real _maxNu) : KernelBase(&_kgrid, 0), parent((_kgrid).getParent()), m_kgrid(_kgrid), m_egrid(_egrid), m_minK(_minK), m_maxK(_maxK), m_minNu(_minNu), m_maxNu(_maxNu) { run(); } KnTurbulenceClamp (const KnTurbulenceClamp& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_kgrid(o.m_kgrid), m_egrid(o.m_egrid), m_minK(o.m_minK), m_maxK(o.m_maxK), m_minNu(o.m_minNu), m_maxNu(o.m_maxNu) {} inline void op(int idx,Grid<Real >& kgrid, Grid<Real >& egrid, Real minK, Real maxK, Real minNu, Real maxNu)  {
    Real eps = egrid[idx];
    Real ke = clamp(kgrid[idx],minK,maxK);
    Real nu = keCmu*square(ke)/eps;
    if (nu > maxNu) 
        eps = keCmu*square(ke)/maxNu;
    if (nu < minNu) 
        eps = keCmu*square(ke)/minNu;

    kgrid[idx] = ke;
    egrid[idx] = eps;
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_kgrid, m_egrid, m_minK, m_maxK, m_minNu, m_maxNu); } FluidSolver* parent; Grid<Real >& m_kgrid; Grid<Real >& m_egrid; Real m_minK; Real m_maxK; Real m_minNu; Real m_maxNu;  }; 

//! clamp k and epsilon to limits    
void TurbulenceClampMesh(VortexSheetMesh& mesh, Real minK, Real maxK, Real minNu, Real maxNu) {
    for (int idx=0; idx<mesh.numNodes(); idx++) {
        Real eps = mesh.turb(idx).epsilon;
        Real ke = clamp(mesh.turb(idx).k,minK,maxK);
        Real nu = keCmu*square(ke)/eps;
        if (nu > maxNu) 
            eps = keCmu*square(ke)/maxNu;
        if (nu < minNu) 
            eps = keCmu*square(ke)/minNu;

        mesh.turb(idx).k = ke;
        mesh.turb(idx).epsilon = eps;
    }
}

//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps



struct KnComputeProductionStrain : public KernelBase {  KnComputeProductionStrain (const MACGrid& _vel, const Grid<Vec3 >& _velCenter, const Grid<Real >& _ke, const Grid<Real >& _eps, Grid<Real >& _prod, Grid<Real >& _nuT, Real _pscale = 1.0f) : KernelBase(&_vel, 1), parent((_vel).getParent()), m_vel(_vel), m_velCenter(_velCenter), m_ke(_ke), m_eps(_eps), m_prod(_prod), m_nuT(_nuT), m_pscale(_pscale) { run(); } KnComputeProductionStrain (const KnComputeProductionStrain& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_vel(o.m_vel), m_velCenter(o.m_velCenter), m_ke(o.m_ke), m_eps(o.m_eps), m_prod(o.m_prod), m_nuT(o.m_nuT), m_pscale(o.m_pscale) {} inline void op(int i,int j,int k,const MACGrid& vel, const Grid<Vec3 >& velCenter, const Grid<Real >& ke, const Grid<Real >& eps, Grid<Real >& prod, Grid<Real >& nuT, Real pscale)  {
    Real curEps = eps(i,j,k);
    if (curEps > 0) {
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(ke(i,j,k)) / curEps;
        
        // compute Sij = 1/2 * (dU_i/dx_j + dU_j/dx_i)
        Vec3 diag = Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, vel(i,j,k+1).z) - vel(i,j,k);
        Vec3 ux = 0.5*(velCenter(i+1,j,k)-velCenter(i-1,j,k));
        Vec3 uy = 0.5*(velCenter(i,j+1,k)-velCenter(i,j-1,k));
        Vec3 uz = 0.5*(velCenter(i,j,k+1)-velCenter(i,j,k-1));
        Real S12 = 0.5*(ux.y+uy.x);
        Real S13 = 0.5*(ux.z+uz.x);
        Real S23 = 0.5*(uy.z+uz.y);
        Real S2 = square(diag.x) + square(diag.y) + square(diag.z) +
                  2.0*square(S12) + 2.0*square(S13) + 2.0*square(S23);
        
        // P = 2*nu_T*sum_ij(Sij^2)
        prod(i,j,k) = 2.0 * curNu * S2 * pscale;
        nuT(i,j,k) = curNu;        
    } 
    else {
        prod(i,j,k) = 0;
        nuT(i,j,k) = 0;
    }
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_vel, m_velCenter, m_ke, m_eps, m_prod, m_nuT, m_pscale); } FluidSolver* parent; const MACGrid& m_vel; const Grid<Vec3 >& m_velCenter; const Grid<Real >& m_ke; const Grid<Real >& m_eps; Grid<Real >& m_prod; Grid<Real >& m_nuT; Real m_pscale;  }; 

//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Omegaij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps



struct KnComputeProductionCurl : public KernelBase {  KnComputeProductionCurl (const Grid<Vec3 >& _curl, const Grid<Vec3 >& _bcurl, const Grid<Real >& _ke, const Grid<Real >& _eps, Grid<Real >& _prod, Grid<Real >& _nuT, Real _pscale = 1.0f, Grid<Vec3 >* _debug = NULL) : KernelBase(&_curl, 1), parent((_curl).getParent()), m_curl(_curl), m_bcurl(_bcurl), m_ke(_ke), m_eps(_eps), m_prod(_prod), m_nuT(_nuT), m_pscale(_pscale), m_debug(_debug) { run(); } KnComputeProductionCurl (const KnComputeProductionCurl& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_curl(o.m_curl), m_bcurl(o.m_bcurl), m_ke(o.m_ke), m_eps(o.m_eps), m_prod(o.m_prod), m_nuT(o.m_nuT), m_pscale(o.m_pscale), m_debug(o.m_debug) {} inline void op(int i,int j,int k,const Grid<Vec3 >& curl, const Grid<Vec3 >& bcurl, const Grid<Real >& ke, const Grid<Real >& eps, Grid<Real >& prod, Grid<Real >& nuT, Real pscale, Grid<Vec3 >* debug)  {
    Real curEps = eps(i,j,k);
    if (curEps > 0) {
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(ke(i,j,k)) / curEps;
        
        // diff with clamping
        Vec3 actual = curl(i,j,k), buoyant=bcurl(i,j,k);
        Vec3 diff = actual - buoyant;
        for (int c=0; c<3; c++) {
            if (actual[c]*diff[c] < 0) 
                diff[c] = 0; // clamp to 0 if buoyant is bigger than actual
            if (actual[c]*buoyant[c] < 0) 
                diff[c] = actual[c]; // ignore if initial dirs point in different direction            
            if (fabs(diff[c]) > fabs(actual[c]))
                diff[c] = actual[c];
        }
        if (debug)
            (*debug)(i,j,k) = diff;
        
        prod(i,j,k) = 2.0 * curNu * 2.0 * normSquare(diff) * pscale;
        nuT(i,j,k) = curNu;                
    } 
    else {
        prod(i,j,k) = 0;
        nuT(i,j,k) = 0;
    }
} void run() { const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ; for (int k=minZ; k < _maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k, m_curl, m_bcurl, m_ke, m_eps, m_prod, m_nuT, m_pscale, m_debug); } FluidSolver* parent; const Grid<Vec3 >& m_curl; const Grid<Vec3 >& m_bcurl; const Grid<Real >& m_ke; const Grid<Real >& m_eps; Grid<Real >& m_prod; Grid<Real >& m_nuT; Real m_pscale; Grid<Vec3 >* m_debug;  }; 
    
//! Compute k-epsilon production term P = 2*nu_T*sum_ij(Sij^2) and the turbulent viscosity nu_T=C_mu*k^2/eps

void KEpsilonComputeProduction(MACGrid& vel,Grid<Real>& k,Grid<Real>& eps,Grid<Real>& prod,Grid<Real>& nuT,Real pscale = 1.0f,Grid<Vec3>* bcurl=NULL,Grid<Vec3>* debug=NULL, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    // get centered velocity grid
    Grid<Vec3> vcenter(parent);
    GetCentered(vcenter, vel);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);    
    KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);
    
    if (bcurl) {
        Grid<Vec3> curl(parent);
        CurlOp(vcenter, curl);
        
        // compute production field
        KnComputeProductionCurl(curl,*bcurl, k, eps, prod, nuT, pscale, debug);
    } else {
        KnComputeProductionStrain(vel, vcenter, k, eps, prod, nuT, pscale);
    }
}PyObject* _plugin_KEpsilonComputeProduction (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonComputeProduction"); PyObject *_retval = NULL; { ArgLocker _lock; MACGrid& vel = *_args.get< MACGrid* > (0,"vel", &_lock); Grid<Real>& k = *_args.get< Grid<Real>* > (1,"k", &_lock); Grid<Real>& eps = *_args.get< Grid<Real>* > (2,"eps", &_lock); Grid<Real>& prod = *_args.get< Grid<Real>* > (3,"prod", &_lock); Grid<Real>& nuT = *_args.get< Grid<Real>* > (4,"nuT", &_lock); Real pscale = _args.getOpt< Real > (5,"pscale", 1.0f, &_lock); Grid<Vec3>* bcurl = _args.getOpt< Grid<Vec3>* > (6,"bcurl", NULL, &_lock); Grid<Vec3>* debug = _args.getOpt< Grid<Vec3>* > (7,"debug", NULL, &_lock); _retval = getPyNone();KEpsilonComputeProduction(vel, k, eps, prod, nuT, pscale, bcurl, debug, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonComputeProduction"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonComputeProduction",e.what()); return 0; } } 

//! Integrate source terms of k-epsilon equation

struct KnAddTurbulenceSource : public KernelBase {  KnAddTurbulenceSource (Grid<Real >& _kgrid, Grid<Real >& _egrid, const Grid<Real >& _pgrid, Real _dt) : KernelBase(&_kgrid, 0), parent((_kgrid).getParent()), m_kgrid(_kgrid), m_egrid(_egrid), m_pgrid(_pgrid), m_dt(_dt) { run(); } KnAddTurbulenceSource (const KnAddTurbulenceSource& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_kgrid(o.m_kgrid), m_egrid(o.m_egrid), m_pgrid(o.m_pgrid), m_dt(o.m_dt) {} inline void op(int idx,Grid<Real >& kgrid, Grid<Real >& egrid, const Grid<Real >& pgrid, Real dt)  {
    Real eps = egrid[idx], prod = pgrid[idx], ke = kgrid[idx];
    if (ke <= 0) ke = 1e-3; // pre-clamp to avoid nan
    
    Real newK = ke + dt*(prod - eps);
    Real newEps = eps + dt*(prod * keC1 - eps * keC2) * (eps / ke);
    if (newEps <= 0) newEps = 1e-4; // pre-clamp to avoid nan

    kgrid[idx] = newK;
    egrid[idx] = newEps;
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_kgrid, m_egrid, m_pgrid, m_dt); } FluidSolver* parent; Grid<Real >& m_kgrid; Grid<Real >& m_egrid; const Grid<Real >& m_pgrid; Real m_dt;  }; 


//! Integrate source terms of k-epsilon equation
void KEpsilonSources(Grid<Real>& k,Grid<Real>& eps,Grid<Real>& prod, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Real dt = parent->getDt();
        
    KnAddTurbulenceSource(k, eps, prod, dt);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);
    KnTurbulenceClamp(k, eps, minK, maxK, keNuMin, keNuMax);    
}PyObject* _plugin_KEpsilonSources (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonSources"); PyObject *_retval = NULL; { ArgLocker _lock; Grid<Real>& k = *_args.get< Grid<Real>* > (0,"k", &_lock); Grid<Real>& eps = *_args.get< Grid<Real>* > (1,"eps", &_lock); Grid<Real>& prod = *_args.get< Grid<Real>* > (2,"prod", &_lock); _retval = getPyNone();KEpsilonSources(k, eps, prod, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonSources"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonSources",e.what()); return 0; } } 

//! Integrate source terms of k-epsilon equation
void AddTurbulenceSourceMesh(VortexSheetMesh& mesh, Grid<Vec3>& curl, Grid<Vec3>& bcurl, Real dt) {
    for (int idx=0; idx<mesh.numNodes(); idx++) {
        Real eps = mesh.turb(idx).epsilon;
        Real k = mesh.turb(idx).k;
        const Vec3& pos = mesh.nodes(idx).pos;
        
        // turbulent viscosity: nu_T = C_mu * k^2/eps
        Real curNu = keCmu * square(k) / eps;
        
        // diff with clamping
        Vec3 actual = curl.getInterpolated(pos), buoyant=bcurl.getInterpolated(pos);
        Vec3 diff = actual - buoyant;
        for (int c=0; c<3; c++) {
            if (actual[c]*diff[c] < 0) 
                diff[c] = 0; // clamp to 0 if buoyant is bigger than actual
            if (actual[c]*buoyant[c] < 0) 
                diff[c] = actual[c]; // ignore if initial dirs point in different direction            
            if (fabs(diff[c]) > fabs(actual[c]))
                diff[c] = actual[c];
        }
        
        // production
        Real prod = 2.0 * curNu * 2.0 * normSquare(diff);
        Real newK = k + dt*(prod - eps);
        Real newEps = eps + dt*(prod * keC1 - eps * keC2) * (eps / k);
        
        mesh.turb(idx).k = newK;
        mesh.turb(idx).epsilon = newEps;
    }    
}

//! Integrate source terms of k-epsilon equation
void KEpsilonSourcesMesh(VortexSheetMesh& mesh,MACGrid& vel,Grid<Vec3>& bcurl,Grid<Real>& kgrid,Grid<Real>& epsGrid, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Real dt = parent->getDt();
    
    // get centered velocity grid
    Grid<Vec3> vcenter(parent);
    GetCentered(vcenter, vel);
    
    // compute limits
    const Real minK = 1.5*square(keU0)*square(keImin);
    const Real maxK = 1.5*square(keU0)*square(keImax);    
    
    Grid<Vec3> curl(parent);
    CurlOp(vcenter, curl);
    
    TurbulenceClampMesh(mesh, minK, maxK, keNuMin, keNuMax);
    AddTurbulenceSourceMesh(mesh, curl, bcurl, dt);
    TurbulenceClampMesh(mesh, minK, maxK, keNuMin, keNuMax);
    
    Grid<Real> sum(parent);
    kgrid.clear();
    epsGrid.clear();
    for (int i=0; i<mesh.numNodes(); i++) {
        const Vec3& p = mesh.nodes(i).pos;
        kgrid.setInterpolated(p, mesh.turb(i).k, sum);
        epsGrid.setInterpolated(p, mesh.turb(i).epsilon, sum);
    }
    sum *= 0.5;
    kgrid.safeDivide(sum);
    epsGrid.safeDivide(sum);
}PyObject* _plugin_KEpsilonSourcesMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonSourcesMesh"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); MACGrid& vel = *_args.get< MACGrid* > (1,"vel", &_lock); Grid<Vec3>& bcurl = *_args.get< Grid<Vec3>* > (2,"bcurl", &_lock); Grid<Real>& kgrid = *_args.get< Grid<Real>* > (3,"kgrid", &_lock); Grid<Real>& epsGrid = *_args.get< Grid<Real>* > (4,"epsGrid", &_lock); _retval = getPyNone();KEpsilonSourcesMesh(mesh, vel, bcurl, kgrid, epsGrid, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonSourcesMesh"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonSourcesMesh",e.what()); return 0; } } 

void KEpsilonInit(Grid<Real>& k,Grid<Real>& eps,Real intensity,Real nu, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    // compute limits
    const Real vk = 1.5*square(keU0)*square(intensity);
    const Real ve = keCmu*square(vk) / nu;
    
    FOR_IDX(k) {
        k[idx] = vk;
        eps[idx] = ve;
    }
}PyObject* _plugin_KEpsilonInit (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonInit"); PyObject *_retval = NULL; { ArgLocker _lock; Grid<Real>& k = *_args.get< Grid<Real>* > (0,"k", &_lock); Grid<Real>& eps = *_args.get< Grid<Real>* > (1,"eps", &_lock); Real intensity = _args.get< Real > (2,"intensity", &_lock); Real nu = _args.get< Real > (3,"nu", &_lock); _retval = getPyNone();KEpsilonInit(k, eps, intensity, nu, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonInit"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonInit",e.what()); return 0; } } 

void KEpsilonInitMesh(VortexSheetMesh& mesh,Real intensity,Real nu, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    // compute limits
    const Real vk = 1.5*square(keU0)*square(intensity);
    const Real ve = keCmu*square(vk) / nu;
    
    for(int i=0; i<mesh.numNodes(); i++) {
        mesh.turb(i).k = vk;
        mesh.turb(i).epsilon = ve;
    }
}PyObject* _plugin_KEpsilonInitMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonInitMesh"); PyObject *_retval = NULL; { ArgLocker _lock; VortexSheetMesh& mesh = *_args.get< VortexSheetMesh* > (0,"mesh", &_lock); Real intensity = _args.get< Real > (1,"intensity", &_lock); Real nu = _args.get< Real > (2,"nu", &_lock); _retval = getPyNone();KEpsilonInitMesh(mesh, intensity, nu, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonInitMesh"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonInitMesh",e.what()); return 0; } } 

//! Compute k-epsilon turbulent viscosity
void KEpsilonGradientDiffusion(Grid<Real>& k,Grid<Real>& eps,Grid<Real>& nuT, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Real dt = parent->getDt();
    MACGrid grad(parent);
    Grid<Real> div(parent);
    
    // gradient diffusion of k
    GradientOpMAC(grad, k);
    grad *= nuT;
    DivergenceOpMAC(div, grad);
    div *= dt/keC1;
    k += div;

    // gradient diffusion of epsilon
    GradientOpMAC(grad, eps);
    grad *= nuT;
    DivergenceOpMAC(div, grad);
    div *= dt/keC2;
    eps += div;
}PyObject* _plugin_KEpsilonGradientDiffusion (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "KEpsilonGradientDiffusion"); PyObject *_retval = NULL; { ArgLocker _lock; Grid<Real>& k = *_args.get< Grid<Real>* > (0,"k", &_lock); Grid<Real>& eps = *_args.get< Grid<Real>* > (1,"eps", &_lock); Grid<Real>& nuT = *_args.get< Grid<Real>* > (2,"nuT", &_lock); _retval = getPyNone();KEpsilonGradientDiffusion(k, eps, nuT, parent);_args.check(); } pbFinalizePlugin(parent,"KEpsilonGradientDiffusion"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("KEpsilonGradientDiffusion",e.what()); return 0; } } 

} // namespace

