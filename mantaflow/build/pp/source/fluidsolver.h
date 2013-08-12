




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/fluidsolver.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Main class for the fluid solver
 *
 ******************************************************************************/

#ifndef _FLUIDSOLVER_H
#define _FLUIDSOLVER_H

#include "pclass.h"
#include "vectorbase.h"
#include <vector>
#include <map>

namespace Manta { 
    
//! Encodes grid size, timstep etc.

class FluidSolver : public PbClass {
public:
     FluidSolver(Vec3i gridSize,int dim=3,int is_debug=1) ;
    virtual ~FluidSolver();
    
    // accessors
    Vec3i getGridSize() { return mGridSize; }PyObject* _getGridSize (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FluidSolver::getGridSize"); { ArgLocker _lock; this->_args.copy(__args); _retval = d_toPy(getGridSize() );this->_args.check(); } pbFinalizePlugin(this->mParent,"FluidSolver::getGridSize"); return _retval; } catch(std::exception& e) { pbSetError("FluidSolver::getGridSize",e.what()); return 0; } } 
    inline Real getDt() { return mDt; }
    inline Real getTime() { return mTimeTotal; }
    inline Real getDx() { return 1.0 / mGridSize.max(); }
    inline Real getScale() { return mScale; }
    //! Check dimensionality
    inline bool is2D() const { return mDim==2; }
    //! Check dimensionality
    inline bool is3D() const { return mDim==3; }
    
    // Python callable methods    
    //! output performace statistics
    void printTimings() ;PyObject* _printTimings (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FluidSolver::printTimings"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();printTimings();this->_args.check(); } pbFinalizePlugin(this->mParent,"FluidSolver::printTimings"); return _retval; } catch(std::exception& e) { pbSetError("FluidSolver::printTimings",e.what()); return 0; } } 
    void saveMeanTimings(std::string filename) ;PyObject* _saveMeanTimings (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FluidSolver::saveMeanTimings"); { ArgLocker _lock; std::string filename = __args.get< std::string > (0,"filename", &_lock); this->_args.copy(__args); _retval = getPyNone();saveMeanTimings(filename);this->_args.check(); } pbFinalizePlugin(this->mParent,"FluidSolver::saveMeanTimings"); return _retval; } catch(std::exception& e) { pbSetError("FluidSolver::saveMeanTimings",e.what()); return 0; } } 
    
    //! Advance the solver one timestep, update GUI if present
    void step() ;PyObject* _step (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FluidSolver::step"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();step();this->_args.check(); } pbFinalizePlugin(this->mParent,"FluidSolver::step"); return _retval; } catch(std::exception& e) { pbSetError("FluidSolver::step",e.what()); return 0; } } 
    
    //! create a object with the solver as its parent
    PbClass* create(PbType type,const std::string& name = "") ;PyObject* _create (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "FluidSolver::create"); { ArgLocker _lock; PbType type = __args.get< PbType > (0,"type", &_lock); const std::string& name = __args.getOpt< std::string > (1,"name", "", &_lock); this->_args.copy(__args); _retval = d_toPy(create(type, name) );this->_args.check(); } pbFinalizePlugin(this->mParent,"FluidSolver::create"); return _retval; } catch(std::exception& e) { pbSetError("FluidSolver::create",e.what()); return 0; } }  //struct PbType {std::string str;}  sample:
    //create(FlagGrid)  create(RealGrid,false)
    
    // temp grid and plugin stuff: you shouldn't call this manually
    template<class T> T* getGridPointer();
    template<class T> void freeGridPointer(T* ptr);    
    void pluginStart(const std::string& name);
    void pluginStop(const std::string& name);        
    int to_debug;
protected:
    //! subclass for managing grid memory
    //! stored as a stack to allow fast allocation
    template<class T> struct GridStorage {
        GridStorage() : used(0) {}
        T* get(Vec3i size);
        void free();
        void release(T* ptr);
        
        std::vector<T*> grids;
        int used;
    };
    
    Vec3i mGridSize;
    const int mDim;
    Real mDt;friend PyObject* _get_FluidSolver_mDt(PyObject* self, void* cl);friend int _set_FluidSolver_mDt(PyObject* self, PyObject* val, void* cl);
    Real mTimeTotal, mScale;
    int mFrame;
        
    GridStorage<int> mGridsInt;
    GridStorage<Real> mGridsReal;
    GridStorage<Vec3> mGridsVec;

    // for timing plugins
    MuTime mPluginTimer;
    std::string mLastPlugin;
    std::vector<std::pair<std::string, MuTime> > mTimings;
    std::map<std::string, std::pair<int,MuTime> > mTimingsTotal;
protected:PbArgs _args;};;

}

#endif


