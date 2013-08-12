




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/vortexfilament.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex filament
 *
 ******************************************************************************/

#ifndef _VORTEXFIL_H
#define _VORTEXFIL_H

#include "particle.h"
#include<set>
#include<map>

namespace Manta {
class Mesh;
struct Segment {
    Vec3 p1,p2;
    int index;
    Segment(Vec3 pp1,Vec3 pp2,int ring_num)
    {
        p1=pp1;
        p2=pp2;
        index=ring_num;
    }
};
    
struct VortexRing {
    VortexRing() : circulation(0.),flag(0),isClosed(false) {}
    VortexRing(Real c, bool closed=false) : circulation(c),flag(0),isClosed(closed) {}
    void renumber(int* _renumber);
    inline int size() const { return indices.size(); }
    inline int idx(int i) const { return indices[(i+indices.size()) % indices.size()]; }
    inline int idx0(int i) const { return indices[i]; }
    inline int idx1(int i) const { return indices[ (i+1) % indices.size() ]; }
    
    bool isClosed;
    int flag;
    Real circulation;
    Real total_circulation;
    Real total_length;
    std::vector<int> indices;
};

//! Vortex filaments
class VortexFilamentSystem : public ConnectedParticleSystem<BasicParticleData,VortexRing> {
public:
    virtual SystemType getType() const { return ParticleBase::FILAMENT; };
        
     VortexFilamentSystem(FluidSolver* parent) ;
  
    //! self-advect the filament system
    void advectSelf(Real scale=1.0,Real regularization=0.1,int integrationMode=IntRK4) ;PyObject* _advectSelf (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::advectSelf"); { ArgLocker _lock; Real scale = __args.getOpt< Real > (0,"scale", 1.0, &_lock); Real regularization = __args.getOpt< Real > (1,"regularization", 0.1, &_lock); int integrationMode = __args.getOpt< int > (2,"integrationMode", IntRK4, &_lock); this->_args.copy(__args); _retval = getPyNone();advectSelf(scale, regularization, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::advectSelf"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::advectSelf",e.what()); return 0; } } 
    //! advect a particle system 
    void advectParticles(TracerParticleSystem& sys,Real scale=1.0,Real regularization=0.1,int integrationMode=IntRK2) ;PyObject* _advectParticles (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::advectParticles"); { ArgLocker _lock; TracerParticleSystem& sys = *__args.get< TracerParticleSystem* > (0,"sys", &_lock); Real scale = __args.getOpt< Real > (1,"scale", 1.0, &_lock); Real regularization = __args.getOpt< Real > (2,"regularization", 0.1, &_lock); int integrationMode = __args.getOpt< int > (3,"integrationMode", IntRK2, &_lock); this->_args.copy(__args); _retval = getPyNone();advectParticles(sys, scale, regularization, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::advectParticles"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::advectParticles",e.what()); return 0; } } 
    //! advect triangle mesh using filaments
    void advectMesh(Mesh& mesh,Real scale=1.0,Real regularization=0.1,int integrationMode=IntRK4) ;PyObject* _advectMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::advectMesh"); { ArgLocker _lock; Mesh& mesh = *__args.get< Mesh* > (0,"mesh", &_lock); Real scale = __args.getOpt< Real > (1,"scale", 1.0, &_lock); Real regularization = __args.getOpt< Real > (2,"regularization", 0.1, &_lock); int integrationMode = __args.getOpt< int > (3,"integrationMode", IntRK4, &_lock); this->_args.copy(__args); _retval = getPyNone();advectMesh(mesh, scale, regularization, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::advectMesh"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::advectMesh",e.what()); return 0; } } 
    //! perform doubly-discrete smoke ring flow update
    //! as in [Weissmann,Pinkall 2009]
    void doublyDiscreteUpdate(Real regularization=0.1) ;PyObject* _doublyDiscreteUpdate (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::doublyDiscreteUpdate"); { ArgLocker _lock; Real regularization = __args.getOpt< Real > (0,"regularization", 0.1, &_lock); this->_args.copy(__args); _retval = getPyNone();doublyDiscreteUpdate(regularization);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::doublyDiscreteUpdate"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::doublyDiscreteUpdate",e.what()); return 0; } } 
    //! remesh long or strongly-curved segments
    void remesh(Real maxLen=3.0,Real minLen=1.0) ;PyObject* _remesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::remesh"); { ArgLocker _lock; Real maxLen = __args.getOpt< Real > (0,"maxLen", 3.0, &_lock); Real minLen = __args.getOpt< Real > (1,"minLen", 1.0, &_lock); this->_args.copy(__args); _retval = getPyNone();remesh(maxLen, minLen);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::remesh"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::remesh",e.what()); return 0; } } 
    //! split the filaments 
    void split_ring(Real cosine_threshold=-1.0,Real dist_threshold=0.0) ;PyObject* _split_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::split_ring"); { ArgLocker _lock; Real cosine_threshold = __args.getOpt< Real > (0,"cosine_threshold", -1.0, &_lock); Real dist_threshold = __args.getOpt< Real > (1,"dist_threshold", 0.0, &_lock); this->_args.copy(__args); _retval = getPyNone();split_ring(cosine_threshold, dist_threshold);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::split_ring"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::split_ring",e.what()); return 0; } } 
    //! reconnect the filaments
    void reconnect_ring(const Real cosine_threshold,const Real dist_threshold) ;PyObject* _reconnect_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::reconnect_ring"); { ArgLocker _lock; const Real cosine_threshold = __args.get< Real > (0,"cosine_threshold", &_lock); const Real dist_threshold = __args.get< Real > (1,"dist_threshold", &_lock); this->_args.copy(__args); _retval = getPyNone();reconnect_ring(cosine_threshold, dist_threshold);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::reconnect_ring"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::reconnect_ring",e.what()); return 0; } } 
    //! hairpin removal
    void merge_adj_edge(Real cosine_threshold=0.8) ;PyObject* _merge_adj_edge (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::merge_adj_edge"); { ArgLocker _lock; Real cosine_threshold = __args.getOpt< Real > (0,"cosine_threshold", 0.8, &_lock); this->_args.copy(__args); _retval = getPyNone();merge_adj_edge(cosine_threshold);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::merge_adj_edge"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::merge_adj_edge",e.what()); return 0; } } 
    void divide_ring(Real edge_length,Real dist_threshold) ;PyObject* _divide_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::divide_ring"); { ArgLocker _lock; Real edge_length = __args.get< Real > (0,"edge_length", &_lock); Real dist_threshold = __args.get< Real > (1,"dist_threshold", &_lock); this->_args.copy(__args); _retval = getPyNone();divide_ring(edge_length, dist_threshold);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::divide_ring"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::divide_ring",e.what()); return 0; } } 
    //!merge two rings filaments
    bool merge_ring(int fir,int sec,const Real cosine_threshold,const Real dist_threshold);
    //! reset the dirty 
    void reset_dirty() { 
        printf("Enter reset_dirty\n");
        dirty.erase(dirty.begin(),dirty.end());
        int size=dirty_edge.size();
        for(int i=0;i<size;i++)
            dirty_edge[i].erase(dirty_edge[i].begin(),dirty_edge[i].end());
        printf("Exit reset_dirty\n");
    }PyObject* _reset_dirty (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::reset_dirty"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();reset_dirty();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::reset_dirty"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::reset_dirty",e.what()); return 0; } } 
    //!Debug 
    void Debug_fun() ;PyObject* _Debug_fun (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::Debug_fun"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();Debug_fun();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::Debug_fun"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::Debug_fun",e.what()); return 0; } } 
    //!Revise the circulation after stretch
    void revise_circulation() ;PyObject* _revise_circulation (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::revise_circulation"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();revise_circulation();this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::revise_circulation"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::revise_circulation",e.what()); return 0; } } 
    //!Decimate rings with number of edges less than num_edge_threshold
    void Decimate_ring(int num_edge_threshold,double min_circum_threshold,double max_circum_threshold) ;PyObject* _Decimate_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::Decimate_ring"); { ArgLocker _lock; int num_edge_threshold = __args.get< int > (0,"num_edge_threshold", &_lock); double min_circum_threshold = __args.get< double > (1,"min_circum_threshold", &_lock); double max_circum_threshold = __args.get< double > (2,"max_circum_threshold", &_lock); this->_args.copy(__args); _retval = getPyNone();Decimate_ring(num_edge_threshold, min_circum_threshold, max_circum_threshold);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::Decimate_ring"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::Decimate_ring",e.what()); return 0; } } 
    //!Direct the motion of filament
    void Direct_motion(Real scale,Real regularization,int integrationMode) ;PyObject* _Direct_motion (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::Direct_motion"); { ArgLocker _lock; Real scale = __args.get< Real > (0,"scale", &_lock); Real regularization = __args.get< Real > (1,"regularization", &_lock); int integrationMode = __args.get< int > (2,"integrationMode", &_lock); this->_args.copy(__args); _retval = getPyNone();Direct_motion(scale, regularization, integrationMode);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::Direct_motion"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::Direct_motion",e.what()); return 0; } } 

    
    //! add a filament ring to the system
    void addRing(const Vec3& position,Real circulation,Real radius,Vec3 normal,int number) ;PyObject* _addRing (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::addRing"); { ArgLocker _lock; const Vec3& position = __args.get< Vec3 > (0,"position", &_lock); Real circulation = __args.get< Real > (1,"circulation", &_lock); Real radius = __args.get< Real > (2,"radius", &_lock); Vec3 normal = __args.get< Vec3 > (3,"normal", &_lock); int number = __args.get< int > (4,"number", &_lock); this->_args.copy(__args); _retval = getPyNone();addRing(position, circulation, radius, normal, number);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::addRing"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::addRing",e.what()); return 0; } } 
    //! add a line filament to the system
    void addLine(const Vec3& p0,const Vec3& p1,Real circulation) ;PyObject* _addLine (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "VortexFilamentSystem::addLine"); { ArgLocker _lock; const Vec3& p0 = __args.get< Vec3 > (0,"p0", &_lock); const Vec3& p1 = __args.get< Vec3 > (1,"p1", &_lock); Real circulation = __args.get< Real > (2,"circulation", &_lock); this->_args.copy(__args); _retval = getPyNone();addLine(p0, p1, circulation);this->_args.check(); } pbFinalizePlugin(this->mParent,"VortexFilamentSystem::addLine"); return _retval; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::addLine",e.what()); return 0; } } 
    
        
    virtual ParticleBase* clone();
protected:
    
    //! Biot-Savart line integration
    void integrate(const std::vector<Vec3>& nodesOld, std::vector<Vec3>& nodesNew, Real scale, Real reg, int integrationMode);
    int to_debug;
    Real max_length;
    Real min_length;
    std::set<int> dirty;
    std::vector< std::set<int> > dirty_edge;
protected:PbArgs _args;};;

} // namespace


#endif


