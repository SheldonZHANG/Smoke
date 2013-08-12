




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/vortexpart.h"
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

#ifndef _VORTEXPART_H
#define _VORTEXPART_H

#include "particle.h"

namespace Manta {
class Mesh;
    
struct VortexParticleData {
    VortexParticleData() : pos(_0),vorticity(_0),sigma(0),flag(0) {}
    VortexParticleData(const Vec3& p, const Vec3& v, Real sig) : pos(p),vorticity(v),sigma(sig),flag(0) {}
    Vec3 pos, vorticity;
    Real sigma;
    int flag;    
    static ParticleBase::SystemType getType() { return ParticleBase::VORTEX; }
};

//! Vortex particles
class VortexParticleSystem : public ParticleSystem<VortexParticleData> {
public:
     VortexParticleSystem(FluidSolver* parent) ;
  
    void advectSelf(Real scale=1.0,int integrationMode=IntRK4) ;PyObject* _advectSelf (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexParticleSystem::advectSelf"); { ArgLocker _lock; Real scale = __args.getOpt< Real > (0,"scale", 1.0, &_lock); int integrationMode = __args.getOpt< int > (1,"integrationMode", IntRK4, &_lock); this->_args.copy(__args); _retval = getPyNone();advectSelf(scale, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexParticleSystem::advectSelf"); return _retval; } catch(std::exception& e) { pbSetError("VortexParticleSystem::advectSelf",e.what()); return 0; } } 
    void applyToMesh(Mesh& mesh,Real scale=1.0,int integrationMode=IntRK4) ;PyObject* _applyToMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexParticleSystem::applyToMesh"); { ArgLocker _lock; Mesh& mesh = *__args.get< Mesh* > (0,"mesh", &_lock); Real scale = __args.getOpt< Real > (1,"scale", 1.0, &_lock); int integrationMode = __args.getOpt< int > (2,"integrationMode", IntRK4, &_lock); this->_args.copy(__args); _retval = getPyNone();applyToMesh(mesh, scale, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexParticleSystem::applyToMesh"); return _retval; } catch(std::exception& e) { pbSetError("VortexParticleSystem::applyToMesh",e.what()); return 0; } } 
    
    virtual ParticleBase* clone();
protected:PbArgs _args;};;

} // namespace


#endif

