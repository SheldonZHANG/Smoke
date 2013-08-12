




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/vortexsheet.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex sheets
 *
 ******************************************************************************/

#ifndef _VORTEXSHEET_H
#define _VORTEXSHEET_H

#include "mesh.h"

namespace Manta {

//! Stores vortex sheet info
struct VortexSheetInfo {
    VortexSheetInfo() : vorticity(0.0), vorticitySmoothed(0.0), circulation(0.0), smokeAmount(1.0), smokeParticles(0.0) {}
    
    Vec3 vorticity;
    Vec3 vorticitySmoothed;
    Vec3 circulation;
    Real smokeAmount, smokeParticles;
};

//! Manages vortex sheet info
struct VorticityChannel : public SimpleTriChannel<VortexSheetInfo> {
    virtual TriChannel* clone() { VorticityChannel* vc = new VorticityChannel(); *vc = *this; return vc;}    
};

//! Manages 3D texture coordinates
struct TexCoord3Channel : public SimpleNodeChannel<Vec3> {
    virtual NodeChannel* clone() { TexCoord3Channel* tc = new TexCoord3Channel(); *tc = *this; return tc; }
    
    void addInterpol(int a, int b, Real alpha) { data.push_back((1.0-alpha)*data[a] + alpha*data[b]);}
    void mergeWith(int node, int delnode, Real alpha) { data[node] = 0.5*(data[node]+data[delnode]); }
};

struct TurbulenceInfo {
    TurbulenceInfo() : k(0.0), epsilon(0.0) {}
    TurbulenceInfo(const TurbulenceInfo& a, const TurbulenceInfo& b, Real alpha) : k((1.0-alpha)*a.k+alpha*b.k), epsilon((1.0-alpha)*a.epsilon+alpha*b.epsilon) {}
    Real k, epsilon;    
};

//! Manages k-epsilon information
struct TurbulenceChannel : public SimpleNodeChannel<TurbulenceInfo> {
    virtual NodeChannel* clone() { TurbulenceChannel* tc = new TurbulenceChannel(); *tc = *this; return tc; }
    
    void addInterpol(int a, int b, Real alpha) { data.push_back(TurbulenceInfo(data[a], data[b], alpha)); }
    void mergeWith(int node, int delnode, Real alpha) { data[node] = TurbulenceInfo(data[node], data[delnode], 0.5); }
};

//! Typed Mesh with a vorticity and 2 texcoord3 channels
class VortexSheetMesh : public Mesh {
public:
     VortexSheetMesh(FluidSolver* parent) ;
    virtual Mesh* clone();
    
    virtual MeshType getType() { return TypeVortexSheet; }    
    
    inline VortexSheetInfo& sheet(int i) { return mVorticity.data[i]; };
    inline Vec3& tex1(int i) { return mTex1.data[i]; }
    inline Vec3& tex2(int i) { return mTex2.data[i]; }
    inline TurbulenceInfo& turb(int i) { return mTurb.data[i]; }
    void setReferenceTexOffset(const Vec3& ref) { mTexOffset = ref; }
    void resetTex1();
    void resetTex2();
    
    void calcCirculation() ;PyObject* _calcCirculation (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexSheetMesh::calcCirculation"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();calcCirculation();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexSheetMesh::calcCirculation"); return _retval; } catch(std::exception& e) { pbSetError("VortexSheetMesh::calcCirculation",e.what()); return 0; } } 
    void calcVorticity() ;PyObject* _calcVorticity (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexSheetMesh::calcVorticity"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();calcVorticity();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexSheetMesh::calcVorticity"); return _retval; } catch(std::exception& e) { pbSetError("VortexSheetMesh::calcVorticity",e.what()); return 0; } } 
    void reinitTexCoords() ;PyObject* _reinitTexCoords (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexSheetMesh::reinitTexCoords"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();reinitTexCoords();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexSheetMesh::reinitTexCoords"); return _retval; } catch(std::exception& e) { pbSetError("VortexSheetMesh::reinitTexCoords",e.what()); return 0; } } 
    
protected:
    Vec3 mTexOffset;
    VorticityChannel mVorticity;
    TexCoord3Channel mTex1, mTex2;
    TurbulenceChannel mTurb;
protected:PbArgs _args;};;

}; // namespace

#endif

