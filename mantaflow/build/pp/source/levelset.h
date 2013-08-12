




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/levelset.h"
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

#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "grid.h"

namespace Manta {
class Mesh;

//! Special function for levelsets
class LevelsetGrid : public Grid<Real> {
public:
     LevelsetGrid(FluidSolver* parent,bool show = true) ;
    
    //! reconstruct the levelset using fast marching
    void reinitMarching(FlagGrid& flags,Real maxTime=4.0,MACGrid* velTransport=NULL,bool ignoreWalls=false,bool correctOuterLayer=true) ;PyObject* _reinitMarching (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "LevelsetGrid::reinitMarching"); { ArgLocker _lock; FlagGrid& flags = *__args.get< FlagGrid* > (0,"flags", &_lock); Real maxTime = __args.getOpt< Real > (1,"maxTime", 4.0, &_lock); MACGrid* velTransport = __args.getOpt< MACGrid* > (2,"velTransport", NULL, &_lock); bool ignoreWalls = __args.getOpt< bool > (3,"ignoreWalls", false, &_lock); bool correctOuterLayer = __args.getOpt< bool > (4,"correctOuterLayer", true, &_lock); this->_args.copy(__args); _retval = getPyNone();reinitMarching(flags, maxTime, velTransport, ignoreWalls, correctOuterLayer);this->_args.check(); } pbFinalizePlugin(this->mParent,"LevelsetGrid::reinitMarching"); return _retval; } catch(std::exception& e) { pbSetError("LevelsetGrid::reinitMarching",e.what()); return 0; } } 
    //! create a triangle mesh from the levelset isosurface
    void createMesh(Mesh& mesh) ;PyObject* _createMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "LevelsetGrid::createMesh"); { ArgLocker _lock; Mesh& mesh = *__args.get< Mesh* > (0,"mesh", &_lock); this->_args.copy(__args); _retval = getPyNone();createMesh(mesh);this->_args.check(); } pbFinalizePlugin(this->mParent,"LevelsetGrid::createMesh"); return _retval; } catch(std::exception& e) { pbSetError("LevelsetGrid::createMesh",e.what()); return 0; } } 
    
    //! union with another levelset
    void join(const LevelsetGrid& o) ;PyObject* _join (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "LevelsetGrid::join"); { ArgLocker _lock; const LevelsetGrid& o = *__args.get< LevelsetGrid* > (0,"o", &_lock); this->_args.copy(__args); _retval = getPyNone();join(o);this->_args.check(); } pbFinalizePlugin(this->mParent,"LevelsetGrid::join"); return _retval; } catch(std::exception& e) { pbSetError("LevelsetGrid::join",e.what()); return 0; } } 
    
    //! completely init levelset from flags
    void initFromFlags(FlagGrid& flags,bool ignoreWalls=false) ;PyObject* _initFromFlags (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "LevelsetGrid::initFromFlags"); { ArgLocker _lock; FlagGrid& flags = *__args.get< FlagGrid* > (0,"flags", &_lock); bool ignoreWalls = __args.getOpt< bool > (1,"ignoreWalls", false, &_lock); this->_args.copy(__args); _retval = getPyNone();initFromFlags(flags, ignoreWalls);this->_args.check(); } pbFinalizePlugin(this->mParent,"LevelsetGrid::initFromFlags"); return _retval; } catch(std::exception& e) { pbSetError("LevelsetGrid::initFromFlags",e.what()); return 0; } } 
    
    static Real invalidTimeValue();
protected:PbArgs _args;};;

} //namespace
#endif

