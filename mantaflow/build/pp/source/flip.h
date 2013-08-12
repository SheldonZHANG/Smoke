




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/flip.h"
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

#ifndef _FLIP_H
#define _FLIP_H

#include "particle.h"
#include "grid.h"
#include "randomstream.h"

namespace Manta {
    
struct FlipData {
    FlipData() : pos(_0),vel(_0),flag(0) {}
    FlipData(const Vec3& p, const Vec3& v) : pos(p),vel(v),flag(0) {}
    Vec3 pos, vel;
    int flag;
    static ParticleBase::SystemType getType() { return ParticleBase::FLIP; }
};

//! FLIP particle system
class FlipSystem : public ParticleSystem<FlipData> {
public:
     FlipSystem(FluidSolver* parent) : ParticleSystem<FlipData>(parent), mOldVel(parent), mRand(1238943) {}
  
    //! Copy velocities from grid with given PIC/FLIP ratio
    void velocitiesFromGrid(FlagGrid& flags,MACGrid& vel,Real flipRatio=0.95) ;PyObject* _velocitiesFromGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FlipSystem::velocitiesFromGrid"); { ArgLocker _lock; FlagGrid& flags = *__args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *__args.get< MACGrid* > (1,"vel", &_lock); Real flipRatio = __args.getOpt< Real > (2,"flipRatio", 0.95, &_lock); this->_args.copy(__args); _retval = getPyNone();velocitiesFromGrid(flags, vel, flipRatio);this->_args.check(); } pbFinalizePlugin(this->mParent,"FlipSystem::velocitiesFromGrid"); return _retval; } catch(std::exception& e) { pbSetError("FlipSystem::velocitiesFromGrid",e.what()); return 0; } } 
	//! Write back velocities to grid
    void velocitiesToGrid(FlagGrid& flags,MACGrid& vel) ;PyObject* _velocitiesToGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FlipSystem::velocitiesToGrid"); { ArgLocker _lock; FlagGrid& flags = *__args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *__args.get< MACGrid* > (1,"vel", &_lock); this->_args.copy(__args); _retval = getPyNone();velocitiesToGrid(flags, vel);this->_args.check(); } pbFinalizePlugin(this->mParent,"FlipSystem::velocitiesToGrid"); return _retval; } catch(std::exception& e) { pbSetError("FlipSystem::velocitiesToGrid",e.what()); return 0; } } 
	//! Ensure minimum/maximum number of particles per cell
    void adjustNumber(MACGrid& vel,FlagGrid& flags,int minParticles=8,int maxParticles=12) ;PyObject* _adjustNumber (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FlipSystem::adjustNumber"); { ArgLocker _lock; MACGrid& vel = *__args.get< MACGrid* > (0,"vel", &_lock); FlagGrid& flags = *__args.get< FlagGrid* > (1,"flags", &_lock); int minParticles = __args.getOpt< int > (2,"minParticles", 8, &_lock); int maxParticles = __args.getOpt< int > (3,"maxParticles", 12, &_lock); this->_args.copy(__args); _retval = getPyNone();adjustNumber(vel, flags, minParticles, maxParticles);this->_args.check(); } pbFinalizePlugin(this->mParent,"FlipSystem::adjustNumber"); return _retval; } catch(std::exception& e) { pbSetError("FlipSystem::adjustNumber",e.what()); return 0; } } 
    //! Mark cells with particles as fluid, otherwise empty
    void markFluidCells(FlagGrid& flags) ;PyObject* _markFluidCells (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FlipSystem::markFluidCells"); { ArgLocker _lock; FlagGrid& flags = *__args.get< FlagGrid* > (0,"flags", &_lock); this->_args.copy(__args); _retval = getPyNone();markFluidCells(flags);this->_args.check(); } pbFinalizePlugin(this->mParent,"FlipSystem::markFluidCells"); return _retval; } catch(std::exception& e) { pbSetError("FlipSystem::markFluidCells",e.what()); return 0; } } 
        
    virtual ParticleBase* clone();
private:
	MACGrid mOldVel;
    RandomStream mRand;
protected:PbArgs _args;};;

} // namespace

#endif

