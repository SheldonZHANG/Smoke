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
PYTHON class VortexFilamentSystem : public ConnectedParticleSystem<BasicParticleData, VortexRing> {
public:
    virtual SystemType getType() const { return ParticleBase::FILAMENT; };
        
    PYTHON VortexFilamentSystem(FluidSolver* parent);
  
    //! self-advect the filament system
    PYTHON void advectSelf(Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK4);
    //! advect a particle system 
    PYTHON void advectParticles(TracerParticleSystem& sys, Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK2);
    //! advect triangle mesh using filaments
    PYTHON void advectMesh(Mesh& mesh, Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK4);
    //! perform doubly-discrete smoke ring flow update
    //! as in [Weissmann,Pinkall 2009]
    PYTHON void doublyDiscreteUpdate(Real regularization=0.1);
    //! remesh long or strongly-curved segments
    PYTHON void remesh(Real maxLen=3.0, Real minLen=1.0);
    //! split the filaments 
    PYTHON void split_ring(Real cosine_threshold,Real dist_threshold,int min_num_of_edge);
    //! reconnect the filaments
    PYTHON void reconnect_ring(const Real cosine_threshold,const Real dist_threshold);
    //! hairpin removal
    PYTHON void merge_adj_edge(Real cosine_threshold,double max_edge_length);
    PYTHON void divide_ring(Real avg_edge_length,Real dist_threshold);
    //!merge two rings filaments
    bool merge_ring(int fir,int sec,const Real cosine_threshold,const Real dist_threshold);
    //! reset the dirty 
    PYTHON void reset_dirty() { 
        printf("Enter reset_dirty\n");
        dirty.erase(dirty.begin(),dirty.end());
        int size=dirty_edge.size();
        for(int i=0;i<size;i++)
            dirty_edge[i].erase(dirty_edge[i].begin(),dirty_edge[i].end());
        printf("Exit reset_dirty\n");
    }
    //!Debug 
    PYTHON void Debug_fun();
    //!Revise the circulation after stretch
    PYTHON void revise_circulation();
    //!Decimate rings with number of edges less than num_edge_threshold
    PYTHON void Decimate_ring(int num_edge_threshold,double min_circum_threshold,double max_circum_threshold);
    //!cyclic interploation
    void cyclic_interpolate(int index,double avg_edge_length,int max_num_of_edge);
    //!resample ring 
    PYTHON void resample_ring(Real max_edge_length,Real avg_edge_length,int max_num_of_edge);
    //!Direct the motion of filament
    PYTHON void Direct_motion(Real scale,Real regularization,int integrationMode);

    
    //! add a filament ring to the system
    PYTHON void addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number);
    //! add a line filament to the system
    PYTHON void addLine(const Vec3& p0, const Vec3& p1, Real circulation);
    
        
    virtual ParticleBase* clone();
protected:
    
    //! Biot-Savart line integration
    void integrate(const std::vector<Vec3>& nodesOld, std::vector<Vec3>& nodesNew, Real scale, Real reg, int integrationMode);
    int to_debug;
    Real max_length;
    Real min_length;
    std::set<int> dirty;
    std::vector< std::set<int> > dirty_edge;
};

} // namespace


#endif
