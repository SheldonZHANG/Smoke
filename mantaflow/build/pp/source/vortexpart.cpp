




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/vortexpart.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex particles
 *
 ******************************************************************************/

#include "vortexpart.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {

// vortex particle effect: (cyl coord around wp)
// u = -|wp|*rho*exp( (-rho^2-z^2)/(2sigma^2) ) e_phi
inline Vec3 VortexKernel(const Vec3& p, const vector<VortexParticleData>& vp, Real scale) {
    Vec3 u(_0);
    for (size_t i=0; i<vp.size(); i++) {
        if (vp[i].flag & ParticleBase::PDELETE) continue;
        
        // cutoff radius
        const Vec3 r = p - vp[i].pos;
        const Real rlen2 = normSquare(r);   
        const Real sigma2 = square(vp[i].sigma);
        if (rlen2 > 6.0 * sigma2 || rlen2 < 1e-8) continue;
        
        // split vortex strength
        Vec3 vortNorm = vp[i].vorticity;
        Real strength = normalize(vortNorm) * scale;
    
        // transform in cylinder coordinate system
        const Real rlen = sqrt(rlen2);
        const Real z = dot(r, vortNorm);
        const Vec3 ePhi = cross(r, vortNorm) / rlen;
        const Real rho2 = rlen2 - z*z;
    
        Real vortex = 0;
        if (rho2 > 1e-10) {
            // evaluate Kernel      
            vortex = strength * sqrt(rho2) * exp (rlen2 * -0.5/sigma2);  
        }
        u += vortex * ePhi;
    }
    return u;
}


struct KnVpAdvectMesh : public ParticleKernelBase {  KnVpAdvectMesh (vector<Node >& _nodes, const vector<VortexParticleData >& _vp, Real _scale) : ParticleKernelBase((_nodes).size()), m_nodes(_nodes), m_vp(_vp), m_scale(_scale), u(size) { run(); } KnVpAdvectMesh (const KnVpAdvectMesh& o) : ParticleKernelBase(o.size), m_nodes(o.m_nodes), m_vp(o.m_vp), m_scale(o.m_scale), u(size) {} inline void op(int i,vector<Node >& nodes, const vector<VortexParticleData >& vp, Real scale, vector<Vec3 >& u)  {
    if (nodes[i].flags & Mesh::NfFixed)
        u[i] = _0;
    else
        u[i] = VortexKernel(nodes[i].pos, vp, scale);
} void run() { const int _sz = size; for (int i=0; i < _sz; i++) op(i, m_nodes, m_vp, m_scale, u); } operator vector<Vec3>  () { return u; } const vector<Vec3 >& getRet() const { return u; } vector<Node >& getArg0() { return m_nodes; } typedef vector<Node > type0; const vector<VortexParticleData >& getArg1() { return m_vp; } typedef vector<VortexParticleData > type1; Real& getArg2() { return m_scale; } typedef Real type2; vector<Node >& m_nodes; const vector<VortexParticleData >& m_vp; Real m_scale; vector<Vec3 > u;  }; 


struct KnVpAdvectSelf : public ParticleKernelBase {  KnVpAdvectSelf (vector<VortexParticleData >& _vp, Real _scale) : ParticleKernelBase((_vp).size()), m_vp(_vp), m_scale(_scale), u(size) { run(); } KnVpAdvectSelf (const KnVpAdvectSelf& o) : ParticleKernelBase(o.size), m_vp(o.m_vp), m_scale(o.m_scale), u(size) {} inline void op(int i,vector<VortexParticleData >& vp, Real scale, vector<Vec3 >& u)  {
    if (vp[i].flag & ParticleBase::PDELETE) 
        u[i] = _0;
    else
        u[i] = VortexKernel(vp[i].pos, vp, scale);
} void run() { const int _sz = size; for (int i=0; i < _sz; i++) op(i, m_vp, m_scale, u); } operator vector<Vec3>  () { return u; } const vector<Vec3 >& getRet() const { return u; } vector<VortexParticleData >& getArg0() { return m_vp; } typedef vector<VortexParticleData > type0; Real& getArg1() { return m_scale; } typedef Real type1; vector<VortexParticleData >& m_vp; Real m_scale; vector<Vec3 > u;  }; 
    
VortexParticleSystem::VortexParticleSystem(FluidSolver* parent) :
    ParticleSystem<VortexParticleData>(parent)
{ 
}

void VortexParticleSystem::advectSelf(Real scale, int integrationMode) {
    KnVpAdvectSelf kernel(mData, scale* getParent()->getDt());
    integratePointSet( kernel, integrationMode);    
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
    KnVpAdvectMesh kernel(mesh.getNodeData(), mData, scale* getParent()->getDt());
    integratePointSet( kernel, integrationMode);    
}

ParticleBase* VortexParticleSystem::clone() {
    VortexParticleSystem* nm = new VortexParticleSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
    return nm;
}

    

} // namespace


