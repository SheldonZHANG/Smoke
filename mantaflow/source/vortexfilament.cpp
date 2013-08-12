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

#include "vortexfilament.h"
#include "integrator.h"
#include "interpol.h"
#include "mesh.h"
#include "quaternion.h"

#include "vectorbase.h"
#include<memory.h>

using namespace std;
namespace Manta {

void VortexRing::renumber(int *_renumber) {
    for (size_t i=0; i<indices.size(); i++)
        indices[i] = _renumber[indices[i]];
}

inline Vec3 FilamentKernel(const Vec3& pos, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) {
    const Real strength = 0.25 / M_PI * scale; //In one paper, it seems that strength is a user-input
    //parameter---sheldon
    const Real a2 = square(reg);
    const Real cutoff2 = square(cutoff);
    //const Real mindist = 1e-6;
    const Real mindist=square(reg);
    Vec3 u(_0);
    //cout<<"initial "<<u<<endl;
    
    for (size_t i=0; i<rings.size(); i++) {
        const VortexRing& r = rings[i];
        if (r.flag & ParticleBase::PDELETE) continue;
        
        const int N = r.isClosed ? (r.size()) : (r.size()-1);
        const Real str = strength * r.circulation;
        for (int j=0; j<N; j++) {
            const Vec3 r0 = fp[r.idx0(j)].pos-pos; //in the vortexring data structure, only circulation and flag and
            //the indices of these position are stored
            const Vec3 r1 = fp[r.idx1(j)].pos-pos;
           
            const Real r0_2 = normSquare(r0), r1_2 = normSquare(r1);
            //cout<<r0_2<<"  "<<r1_2<<"  "<<cutoff2<<"  "<<mindist<<endl;
            if (r0_2 > cutoff2 || r1_2 > cutoff2 || r0_2 < mindist || r1_2 < mindist)
                continue;
            
            //u += A * cp; //what we should care is u is a vector rather than a scalar, so when we integrate we should
            //take care of this

            /*const Real res1=r1_2/sqrt(a2+r1_2);
            const Real res0=r0_2/sqrt(a2+r0_2);
            const Real denom=str*dot(r0,r1)/(a2*normSquare(r1-r0)+normSquare(cross(r1,r0)));
            const Vec3 res=(res0-res1)*denom*cross(r1,r0);
            u+=res;*/
            const Real denom1=(normSquare(r1-r0)*a2+normSquare(cross(r1,r1-r0)));
            const Real denom0=(normSquare(r1-r0)*a2+normSquare(cross(r0,r1-r0)));
            if(denom1==0||denom0==0)
            {
                continue;
            }
            const Vec3 res1=dot(r1,r1-r0)*cross(r1,r1-r0)/sqrt(a2+r1_2)/(normSquare(r1-r0)*a2+normSquare(cross(r1,r1-r0)));
            const Vec3 res0=dot(r0,r1-r0)*cross(r0,r1-r0)/sqrt(a2+r0_2)/(normSquare(r1-r0)*a2+normSquare(cross(r0,r1-r0)));
            u+=str*(res1-res0);
            //const Real res1=dot(r1,r0-r1)/r1n/(a2*normSquare(r1-r0)+normSquare(cross(r1,(r0-r1))));
            //const Real res0=dot(r0,r0-r1)/r0n/(a2*normSquare(r1-r0)+normSquare(cross(r0,(r0-r1))));
            //const Vec3 res=str*(res1*cross(r1,r0-r1)-res0*cross(r0,r0-r1));
            /*const Vec3 e = getNormalized(r1-r0);
            const Real r0n = 1.0f/sqrt(a2+r0_2);
            const Real r1n = 1.0f/sqrt(a2+r1_2);
            const Vec3 cp = cross(r0,e);
            const Real A = str * (dot(r1,e)*r1n - dot(r0,e)*r0n) / (a2 + normSquare(cp));
            u += A * cp;*/
        }
    }
    return u;
}
inline Vec3 Direct_FilamentKernel(const Vec3& pos, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) {
    return Vec3(0.0,1.0,0.0);
}

KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnFilamentAdvectParts(vector<BasicParticleData>& nodes, vector<BasicParticleData>& fp, const vector<VortexRing>& rings, Real reg, Real cutoff, Real scale) {
    //cout<<"Enter KnFilamentAdvectParts i="<<i<<endl;
    if (nodes[i].flag & ParticleBase::PDELETE)
        u[i] = _0;
    else
        u[i] = FilamentKernel(nodes[i].pos, rings, fp, reg, cutoff, scale);
}


KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> KnFilamentAdvectMesh(vector<Node>& nodes, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) {
    if (nodes[i].flags & Mesh::NfFixed)
        u[i] = _0;
    else
        u[i] = FilamentKernel(nodes[i].pos, rings, fp, reg, cutoff, scale);
}

/*void VortexFilamentSystem::addRing(Mesh& mesh,Vec3 center, Real radius,int num)
{
    std::set<int> mark;
    std::vector<int> indices; 
    int mesh_node_size=mesh.numNodes();
    for(int i=0;i<num;i++)
    {
        Vec3 pos=center+Vec3(radius*sin(2*M_PI/num),radius*cos(2*M_PI/num,0)); 
        Real min_dist=10000000000.0;
        int index=-1;
        for(int j=0;j<mesh_node_size;j++)
        {
            if(mark.find(j)!=mark.end())
                continue;
            Vec3 node=mesh.getNode(j);
            Real dist=normSquare(pos-node);
            if(dist<min_dist)
            {
                min_dist=dist,index=j;
            }
        }
        if(index!=-1)
        {
            indices.push_back(index);
            mark.insert(index);
        }
    }
}*/
KERNEL(pts) returns(vector<Vec3> u(size))
vector<Vec3> Direct_KnFilamentAdvectParts(vector<BasicParticleData>& nodes, vector<BasicParticleData>& fp, const vector<VortexRing>& rings, Real reg, Real cutoff, Real scale) {
    //cout<<"Enter KnFilamentAdvectParts i="<<i<<endl;
    if (nodes[i].flag & ParticleBase::PDELETE)
        u[i] = _0;
    else
        u[i] = Direct_FilamentKernel(nodes[i].pos, rings, fp, reg, cutoff, scale);
}
void VortexFilamentSystem::Direct_motion(Real scale,Real regularization,int integrationMode)
{
    if(to_debug)
        cerr<<"********VortexFilamentSystem::direct_motion "<<endl;
    Direct_KnFilamentAdvectParts kernel(mData,mData,mSegments,regularization,1e10,scale*getParent()->getDt());
    integratePointSet(kernel,integrationMode);
    if(to_debug)
    cerr<<"*****Exit VorexFilamentSystem::direct_motion"<<endl; 
}
void VortexFilamentSystem::advectSelf(Real scale, Real regularization, int integrationMode) {
    if(to_debug)
    cerr<<"******VortexFilamentSystem::advectSelf   "<<"  mData.size()= "<<mData.size()<<"   Segments.size()= "<<mSegments.size()<<endl;
    KnFilamentAdvectParts kernel(mData, mData, mSegments , regularization, 1e10, scale * getParent()->getDt()); //time step
    //mSegmetns is used to store the rings
    //mData stores the position
    integratePointSet( kernel, integrationMode);
    if(to_debug)
    cerr<<"*****Exit VorexFilamentSystem::advectSelf"<<endl;
}

void VortexFilamentSystem::advectMesh(Mesh& mesh, Real scale, Real regularization, int integrationMode) {
    if(to_debug)
    cerr<<"******VortexFilamentSystem::advectMesh"<<endl; 
    cout<<"******the number of nodes in mesh "<<mesh.numNodes()<<endl;
    KnFilamentAdvectMesh kernel(mesh.getNodeData(), mSegments, mData, regularization, 1e10, scale * getParent()->getDt());
    integratePointSet( kernel, integrationMode);   
    if(to_debug)
    cerr<<"******Exit VortexFilamentSystem::advectMesh"<<endl;
}

void VortexFilamentSystem::advectParticles(TracerParticleSystem& sys, Real scale, Real regularization, int integrationMode) {
    if(to_debug)
    cerr<<"******VortexFilamentSystem::advectParticles(TracerPatricleSystem "<<"size of tracer "<<(sys.getData()).size()<<"The size of mData"<<mData.size()<<"  "<<sys.infoString()<<"  "<<sys.getType()<<endl;
    std::vector<BasicParticleData> old_pos=sys.getData();
    KnFilamentAdvectParts kernel(sys.getData(), mData, mSegments, regularization, 1e10, scale * getParent()->getDt());
    //sys.getData() {return mData;} // here mData may be exactly the same as mData, make a copy of the mData
    integratePointSet( kernel, integrationMode);   
    if(to_debug)
    cerr<<"******Exit VortexFilamentSystem::advectParticles(TracerParticleSystem "<<endl;
    std::vector<BasicParticleData> particle_pos=sys.getData();
    int size=old_pos.size();
    /*if(to_debug)
    for(int i=0;i<size;i++)
        {
            //cout<<sys.getType()<<endl; //ERROR From here we could get that the Type is PARTICLE rather than TRACER
            //cout<<particle_pos[i].pos<<"  "<<old_pos[i].pos<<endl;
            Vec3 a=particle_pos[i].pos-old_pos[i].pos;
            float x=abs(a.max());
            //if(x>0.0002||abs(a.min())>0.0002)
                //cout<<"motion"<<endl;
             //   else
                    //cout<<"nomotion"<<endl;
        }
    cout<<endl;*/
}


void VortexFilamentSystem::split_ring(Real cosine_threshold,Real dist_threshold){
    cout<<"Enter split function"<<endl;
    int seg_size=segSize();
    for(int i=0;i<seg_size;i++){
        VortexRing& r=mSegments[i];
        int size=r.size();
        for(int j=0;j<size;j++){
            bool mark=false;
            for(int offset=1;offset<size;offset++){
                int k=j+offset; 
                if(k%size==j) continue;
                const Vec3 fir=mData[r.idx1(j)].pos-mData[r.idx0(j)].pos;
                const Vec3 mid1=(mData[r.idx1(j)].pos+mData[r.idx0(j)].pos)/2;
                const Vec3 sec=mData[r.idx1(k)].pos-mData[r.idx0(k)].pos;
                const Vec3 mid2=(mData[r.idx1(k)].pos+mData[r.idx0(k)].pos)/2;
                const Real cosine_res=dot(getNormalized(fir),getNormalized(sec));
                if(cosine_res>cosine_threshold)
                    continue;
                const Real dist1=normSquare(cross(mid1-mData[r.idx0(k)].pos,sec))/normSquare(sec); 
                const Real dist2=normSquare(cross(mid2-mData[r.idx0(j)].pos,fir))/normSquare(fir);
                Real dist=(sqrt(dist1)+sqrt(dist2))/2;
                if(normSquare(mid1-mid2)<dist_threshold*dist_threshold)
                   dist=sqrt(normSquare(mid1-mid2)); 
                if(dist>dist_threshold)
                    continue;
                //if(((k-j)>5&&(size-k+j)>5)==false)
                //    continue;
                vector<int> index1;
                for(int sub=j+1;sub<=k;sub++)
                    index1.push_back(r.indices[sub]);
                VortexRing new_ring;
                new_ring.isClosed=true;
                for(int sub=k+1;sub<=size+j;sub++)
                    new_ring.indices.push_back(r.indices[sub%size]);
                new_ring.indices.push_back(add(BasicParticleData(mid1)));
                new_ring.indices.push_back(add(BasicParticleData(mid2)));
                
                for(int sub=0;sub<=(k-j-1);sub++)
                    r.indices[sub]=index1[sub];
                if(r.indices.size()<k-j+2)
                    r.indices.resize(k-j+2);
                r.indices[k-j]=add(BasicParticleData(mid2));
                r.indices[k-j+1]=add(BasicParticleData(mid1));
                r.indices.resize(k-j);
                mSegments.push_back(new_ring); 
                dirty.insert(mSegments.size());
                dirty.insert(i); 
                mark=true;
                break;

            }
            if(mark==true)
                break;
        }
        }
        cout<<"Exit split function"<<endl;
}
bool VortexFilamentSystem::merge_ring(int fir,int sec,const Real cosine_threshold,const Real
dist_threshold)
{
    if(dirty.find(fir)!=dirty.end()||dirty.find(sec)!=dirty.end())
        return false;
    //cout<<"Enter merge_ring"<<endl;
    int len1=mSegments[fir].size();
    int len2=mSegments[sec].size();
    for(int i=0;i<len1;i++)
    {
        const Vec3 edge1=mData[mSegments[fir].idx1(i)].pos-mData[mSegments[fir].idx0(i)].pos;
        const Vec3 mid1=(mData[mSegments[fir].idx1(i)].pos+mData[mSegments[fir].idx0(i)].pos)/2;
        for(int j=0;j<len2;j++)
        {
            const Vec3 edge2=mData[mSegments[sec].idx1(j)].pos-mData[mSegments[sec].idx0(j)].pos;
            const Vec3 mid2=(mData[mSegments[sec].idx1(j)].pos+mData[mSegments[sec].idx0(j)].pos)/2;
            const Real cosine_res=dot(getNormalized(edge1),getNormalized(edge2));
            if(cosine_res>cosine_threshold)
                continue;
            //cout<<"normSquare(edge2) "<<normSquare(edge2)<<" normSquare(edge1)="<<normSquare(edge1)<<endl;
            Real dist1=normSquare(cross(mid1-mData[mSegments[sec].idx0(j)].pos,edge2))/normSquare(edge2);
            Real dist2=normSquare(cross(mid2-mData[mSegments[fir].idx0(i)].pos,edge1))/normSquare(edge1); 
            Real dist=(sqrt(dist1)+sqrt(dist2))/2;
            Real dist3=sqrt(normSquare(mid2-mid1));
            if(dist3<dist)
               dist=dist3; 
            if(dist<dist_threshold)
            {
                cout<<"enter if(dist<dist_threshold"<<endl;
                VortexRing new_ring;
                for(int sub=0;sub<=i;sub++)
                    new_ring.indices.push_back(mSegments[fir].indices[sub]); 
                for(int sub=0;sub<len2;sub++)
                    new_ring.indices.push_back(mSegments[sec].indices[(j+1+sub)%len2]);
                for(int sub=i+1;sub<len1;sub++)
                    new_ring.indices.push_back(mSegments[fir].indices[sub]);
                new_ring.isClosed=mSegments[fir].isClosed;
                new_ring.flag=mSegments[fir].flag;
                new_ring.circulation=mSegments[fir].circulation;
                mSegments.push_back(new_ring);
                dirty.insert(fir);
                dirty.insert(sec);
                dirty.insert(segSize());
                return true;

                /*std::vector<int> new_index;
                for(int sub=0;sub<=i;sub++)
                  new_index.push_back(mSegments[fir].indices[sub]); 
                for(int sub=0;sub<len2;sub++)
                    new_index.push_back(mSegments[sec].indices[(j+1+sub)%len2]);
                for(int sub=i+1;sub<len1;sub++)
                    new_index.push_back(mSegments[sec].indices[sub]);
                int new_size=new_index.size();
                cout<<"Before acessing mSegments"<<endl;
                for(int i=0;i<new_size;i++)
                    mSegments[fir].indices[i]=new_index[i];
                cout<<"before resize mSegments new_size "<<new_size<<" oldsize "<<len1<<endl; 
                VortexRing& refer=mSegments[fir];
                for(int i=0;i<len2;i++)
                    refer.indices.push_back(i);
                cout<<"len1 "<<len1<<" len2 "<<len2<<endl;
                refer.indices.resize(len1+len2);
                for(int sub=0;sub<len1+len2;sub++)
                    refer.indices[sub]=new_index[sub];
                cout<<"Exit merge_ring true"<<endl; 
                return true;*/
            }
                
        }
    }
    //cout<<"Exit merge_ring false"<<endl;
    return false;
}
void VortexFilamentSystem::Debug_fun()
{
    int seg_size=segSize();
    for(int i=0;i<seg_size;i++)
    {
        VortexRing ring=mSegments[i];
        int index_size=ring.indices.size();
        for(int j=0;j<index_size;j++)
        {
            cout<<ring.indices[j]<<",";
            
        }
        cout<<endl;
    }
}
void VortexFilamentSystem::Decimate_ring(int num_edge_threshold,double min_circum_threshold,double max_circum_threshold)
{
    cout<<"Enter Decimate_ring"<<endl;
    int seg_size=segSize();
    int new_size=seg_size;
    std::vector<bool> mark(seg_size);
    for(int i=0;i<seg_size;i++)
    {
        VortexRing& ring=mSegments[i];
        int num_edge=ring.size();
        /*if(num_edge<num_edge_threshold)
        {
           //new_size--;
           //mark[i]=true; 
        }
        else*/
        {
            Real circum=0;
            for(int j=0;j<num_edge;j++)
                circum+=sqrt(normSquare(mData[ring.idx0(j)].pos-mData[ring.idx1(j)].pos));
            if(circum<min_circum_threshold)
            {
                mark[i]=true;
                new_size--;
                for(int j=0;j<num_edge;j++)
                    mData[ring.idx0(j)].flag|=PDELETE;
            }
            else if(circum>max_circum_threshold)
            {
                mark[i]=true;
                new_size--;
                for(int j=0;j<num_edge;j++)
                    mData[ring.idx0(j)].flag|=PDELETE;
            }
        }
        if(ring.isClosed==false)
        {
            cout<<"Error unclosed ring with "<<num_edge<<" edges"<<endl;
            exit(0);
        }
    }
    for(int i=0,copyFrom=0;i<new_size;i++,copyFrom++)
    {
       while(mark[copyFrom])
           copyFrom++;
       if(i!=copyFrom)
          mSegments[i]=mSegments[copyFrom]; 
    }
    mSegments.resize(new_size);
    compress();
    cout<<"Exit Decimate_ring"<<endl;
}
void VortexFilamentSystem::reconnect_ring(const Real cosine_threshold,const Real dist_threshold)
{
    cout<<"Enter reconnect segSize="<<segSize()<<endl;
    int seg_size=segSize();
    std::vector<Vec3> avg(seg_size);
    std::vector<Real> radius(seg_size);
    for(int i=0;i<seg_size;i++)
    {
        VortexRing& ring=mSegments[i];
        int index_size=ring.size();
        for(int j=0;j<index_size;j++)
        {
           avg[i]=avg[i]+mData[ring.indices[j]].pos; 
        }
        avg[i]=avg[i]/index_size;
        for(int j=0;j<index_size;j++)
           radius[i]=max(radius[i],normSquare(avg[i]-mData[ring.indices[j]].pos)); 
    }
    int new_size=seg_size;
    std::vector<bool> mark(seg_size);
    for(int i=0;i<seg_size;i++)
        mark[i]=false;
    for(int i=0;i<seg_size;i++)
    {
        if(mark[i]==true||dirty.find(i)!=dirty.end())
            continue;
        for(int j=i+1;j<seg_size;j++)
        {
            if(mark[j]==true||dirty.find(j)!=dirty.end())
                continue;
            if(sqrt(normSquare(avg[i]-avg[j]))>(sqrt(radius[i])+sqrt(radius[j])))
                continue;
            bool res=merge_ring(i,j,cosine_threshold,dist_threshold);
            if(res==true)
            {
                cout<<"connect "<<i<<" and "<<j<<endl;
                mark[j]=true;
                mark[i]=true;
                //dirty.insert(i);
                //dirty.insert(j);
                //dirty.insert(segSize());
                new_size--;
                break;
            }
        }
    }
    for(int i=0,copyFrom=0;i<new_size;i++,copyFrom++)
    {
       while(mark[copyFrom])
           copyFrom++;
       if(i!=copyFrom)
          mSegments[i]=mSegments[copyFrom]; 
    }
    mSegments.resize(new_size);
    cout<<"Exit reconnect segSize "<<segSize()<<" seg_size "<<seg_size<<" new_size "<<new_size<<endl;
}
void VortexFilamentSystem::divide_ring(Real avg_edge_length,Real max_edge_length){
    cout<<"Enter divide_ring function "<<endl;
    int seg_size=segSize();
    for(int i=0;i<seg_size;i++)
    {
        VortexRing& r=mSegments[i];
        const int oldLen=r.size();
        std::vector<int> new_index;
        for(int j=0;j<oldLen;j++)
        {
            new_index.push_back(r.idx0(j));
            const Vec3 p0=mData[r.idx0(j)].pos;
            const Vec3 p1=mData[r.idx1(j)].pos;
            const Real dist=sqrt(normSquare(p1-p0));
            if(dirty_edge[i].find(r.idx0(j))==dirty_edge[i].end()&&dist>max_edge_length)
            {
                const Vec3 dif=p1-p0;
                int num=floor(dist/avg_edge_length);
                for(int k=1;k<num;k++)
                {
                    int pos=add(BasicParticleData(p0+k/num*dif));
                    new_index.push_back(pos);
                }            
                cout<<"Divide one ring edge"<<endl;
            }

        }
        r.indices=new_index; 
    }
}
void VortexFilamentSystem::merge_adj_edge(Real cosine_threshold,double max_edge_length){
    cout<<"Enter merge_adj_edge"<<endl;
    int seg_size=segSize();
    dirty_edge.resize(seg_size);
    for(int sub=0;sub<seg_size;sub++)
    {
        VortexRing& r=mSegments[sub];
        const int oldLen=r.size();
        int newLen=oldLen;
        std::vector<bool> mark(oldLen);
        for(int i=0;i<oldLen;i++){
            if(mData[r.idx0(i)].flag&PDELETE||mData[r.idx1(i)].flag&PDELETE||mData[r.idx1(i+1)].flag&PDELETE) 
                continue;
            const Vec3 p0=mData[r.idx0(i)].pos;
            const Vec3 p1=mData[r.idx1(i)].pos;
            const Vec3 p2=mData[r.idx1((i+1)%oldLen)].pos;
            const Real cosine_res=dot(getNormalized(p1-p0),getNormalized(p2-p1));
            const Real distance=sqrt(normSquare(p0-p2));
            if(p0==p1)
            {
                mData[r.idx1(i)].flag|=PDELETE;
                mark[(i+1)%oldLen]=true;
                newLen--;
                dirty_edge[sub].insert(i);
            }
            else if(cosine_res>cosine_threshold||distance<max_edge_length){
                mData[r.idx1(i)].flag|=PDELETE; 
                mark[(i+1)%oldLen]=true;
                newLen--;
                dirty_edge[sub].insert(i);
                }
        }
        for (int j=0, copyFrom=0; j<newLen; j++,copyFrom++) {
            while (mark[copyFrom]) 
                    copyFrom++;
                if (j!=copyFrom)
                    r.indices[j] = r.indices[copyFrom];
        }            
        r.indices.resize(newLen);
    }
    compress();
    cout<<"Exit merge_adj_edge"<<endl;
    return ;
}

void VortexFilamentSystem::remesh(Real maxLen, Real minLen) {
    cout<<"Enter VortexFilamentSystem:remesh"<<endl;
    const Real maxLen2 = maxLen*maxLen, minLen2 = minLen*minLen;
    max_length=maxLen;
    min_length=minLen;
    
    for (int i=0; i < segSize(); i++) {
        VortexRing& r = mSegments[i];
            
        // insert edges
        for(;;) {
            const int oldLen = r.size();
            map<int,int> insert;        
            int offset = 1;
            
            for (int j=0; j<oldLen; j++) {
                const Vec3 p0 = mData[r.idx0(j)].pos;
                const Vec3 p1 = mData[r.idx1(j)].pos;
                const Real l2 = normSquare(p1-p0); //return x*x+y*y+z*z
                
                if (l2 > maxLen2) {
                    // insert midpoint
                    const Vec3 p_1 = mData[r.idx(j-1)].pos;
                    const Vec3 p2 = mData[r.idx(j+2)].pos;
                    const Vec3 mp = hermiteSpline(p0,p1,crTangent(p_1,p0,p1),crTangent(p0,p1,p2), 0.5);
                    insert.insert(pair<int,int>(j+offset, add(mp)));
                    offset++;
                }
            }
            if (insert.empty()) 
                break;
            
            // renumber indices
            const int newLen = oldLen + insert.size();
            int num=oldLen-1;
            r.indices.resize(newLen);
            for (int j=newLen-1; j>=0; j--) {
                map<int,int>::const_iterator f = insert.find(j);
                if (f==insert.end())
                    r.indices[j] = r.indices[num--];
                else
                    r.indices[j] = f->second;
            }
        }
        
        // remove edges
        for(;;) {
            const int oldLen = r.size();
            const int N = r.isClosed ? oldLen : (oldLen-1);
            std::vector<bool> deleted(r.size());
            
            int newLen=oldLen;
            for (int j=0; j<N; j++) {
                if (mData[r.idx0(j)].flag & PDELETE || mData[r.idx1(j)].flag & PDELETE) continue;
                const Vec3 p0 = mData[r.idx0(j)].pos;
                const Vec3 p1 = mData[r.idx1(j)].pos;
                const Real l2 = normSquare(p1-p0);
                
                if (l2 < minLen2) {
                    // kill edge
                    mData[r.idx0(j)].flag |= PDELETE;
                    mData[r.idx1(j)].pos = 0.5*(p0+p1);
                    deleted[j] = true;
                    newLen--;
                    j++;
                }
            }
            if (newLen == oldLen)
                break;
            
            // renumber indices
            for (int j=0, copyFrom=0; j<newLen; j++,copyFrom++) {
                while (deleted[copyFrom]) 
                    copyFrom++;
                if (j!=copyFrom)
                    r.indices[j] = r.indices[copyFrom];
            }            
            r.indices.resize(newLen);
        }
    }
    
    // remove deleted particles
    compress();
    cout<<"Exit VortexFilamentSystem:remesh"<<endl;
}

VortexFilamentSystem::VortexFilamentSystem(FluidSolver* parent) :
    ConnectedParticleSystem<BasicParticleData, VortexRing>(parent)
{     

    to_debug=parent->to_debug;
    cout<<"Exit sheldon's comment: parent->to_debug"<<parent->to_debug<<endl;
    
}

ParticleBase* VortexFilamentSystem::clone() {
    VortexFilamentSystem* nm = new VortexFilamentSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->mSegments = mSegments;
    nm->setName(getName());
    return nm;
}

// ------------------------------------------------------------------------------
// Functions needed for doubly-discrete smoke flow using Darboux transforms
// see [Weissmann,Pinkall 2009]
// doesn't really work yet (can't reverse rotation dir)
// ------------------------------------------------------------------------------

Real evaluateRefU(int N, Real L, Real circ, Real reg) {
    // construct regular n-polygon
    const Real l = L/(Real)N; //Here L is the length of the polygon--sheldon
    const Real r = 0.5*l/sin(M_PI/(Real)N);
    cout << r << " " << l << endl;
    // build vortex ring
    VortexRing ring (circ);
    vector<BasicParticleData> pos(N);
    for(int i=0; i<N; i++) {
        pos[i].pos = Vec3( r*cos(2.0*M_PI*(Real)i/N), r*sin(2.0*M_PI*(Real)i/N), 0);
        pos[i].flag =0;
        ring.indices.push_back(i);
    }
    
    // Build kernel
    vector<VortexRing> rings;
    rings.push_back(ring);
    
    // evaluate impact on pos[0]
    return norm(FilamentKernel(pos[0].pos, rings, pos, reg, 1e10, 1.0));    
    //inline Vec3 FilamentKernel(const Vec3& pos, const vector<VortexRing>& rings, const vector<BasicParticleData>& fp, Real reg, Real cutoff, Real scale) 
}

Vec3 darbouxStep(const Vec3& Si, const Vec3& lTi, Real r) {
    Quaternion rlTS (lTi - Si, -r);
    Quaternion lT (lTi, 0);
    Quaternion lTnext = rlTS * lT * rlTS.inverse();
    return lTnext.imag();
}

Vec3 monodromy(const vector<Vec3>& gamma, const Vec3& lT_1, Real r) {
    const int N = gamma.size();
    Vec3 lT (lT_1);
    
    for (int i=0; i<N; i++) {
        Vec3 Si = gamma[(i+1)%N]-gamma[i];
        lT = darbouxStep(Si, lT, r);
    }
    return lT;
}

bool powerMethod(const vector<Vec3>& gamma, Real l, Real r, Vec3& lT) {
    const int maxIter = 100;
    const Real epsilon = 1e-4;
    
    for (int i=0; i<maxIter; i++) {
        Vec3 lastLT (lT);
        lT = monodromy(gamma, lT, r);
        //if ((i%1) == 0) cout << "iteration " << i << " residual: " << norm(lT-lastLT) << endl;
        if (norm(lT-lastLT) < epsilon) 
            return true;
    }   
    return false;
}

bool darboux(const vector<Vec3>& from, vector<Vec3>& to, Real l, Real r) {
    const int N = from.size();
    Vec3 lT(0,0,l);
    if (!powerMethod(from, l, r, lT)) // r stands for the real part in contrast to the image part
        return false;
    cout << "iniLT " << lT << " norm " << lT/l<< endl;
    
    for (int i=0; i<N; i++) {
        to[i] = from[i] + lT;
        Vec3 Si = from[(i+1)%N] - from[i];
        lT = darbouxStep(Si, lT, r);
    }
    return true;
}
        

void VortexFilamentSystem::doublyDiscreteUpdate(Real reg) {
    const Real dt = getParent()->getDt();
    
    for (int rc=0; rc<segSize(); rc++) {
        if (!isSegActive(rc) || !mSegments[rc].isClosed) continue;
        
         VortexRing& r = mSegments[rc];
         int N = r.size();
        
        // compute arc length
        Real L=0;
        for (int i=0; i<N; i++)
            L += norm(mData[r.idx0(i)].pos - mData[r.idx1(i)].pos);
        
        // build gamma
        vector<Vec3> gamma(N); //here note that gamma stands for the position coordinate
        for (int i=0; i<N; i++) gamma[i] = mData[r.indices[i]].pos;
        
        //N=1000; L=2.0*M_PI; reg=0.1; r.circulation=1;
        
        // compute reference parameters 
        //the correct speed is obtained when 2d=delt_t(U-U~)
        //with a choice of the rotation angle =2pi/n, the parameters turn out to be
        //l=sqrt(power(L/n,2)+power(d,2)), r=dcot(pi/n);   this exceprt is from Real-time interactive simulation of
        //smoke using discrete integrable vortex filaments
        const Real U = 0.5*r.circulation/L * (log(4.0*L/(M_PI*reg)) - 1.0);
        const Real Ur = evaluateRefU(N, L, r.circulation, reg);// it seems what evaluateRefU actually does is to use
        //r.circulation to generate one ring (denoted as R) instead of using the ring r. Then compute the velocity of the first
       //node in the R which is influenced by other nodes in the R, not including other nodes of other rings. 
        const Real d = 0.5*dt*(U-Ur);
        const Real l = sqrt( square(L/N) + square(d) ); //l=sqrt((L/n)^2+d^2)
        const Real ra = d*tan(M_PI * (0.5 - 1.0/N)); // d*cot(pi/n)
        cout << U << " <-< " << Ur << endl;
        
        // fwd darboux transform
        vector<Vec3> eta(N);
        if (!darboux(gamma, eta, l, ra)) {
            cout << "Fwd Darboux correction failed, skipped." << endl; 
            continue;
        }
        
        // bwd darboux transform
        if (!darboux(eta, gamma, l, ra)) {
            cout << "Bwd Darboux correction failed, skipped." << endl; 
            continue;
        }
        
        // copy back
        for (int i=0; i<N; i++) {
            mData[r.indices[i]].pos = gamma[i];
        }
    }
}

void VortexFilamentSystem::addLine(const Vec3& p0, const Vec3& p1, Real circulation) {
    VortexRing ring(circulation, false);
    
    ring.indices.push_back(add(BasicParticleData(p0)));
    ring.indices.push_back(add(BasicParticleData(p1)));
    mSegments.push_back(ring);
}

void VortexFilamentSystem::addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number) {
    if(to_debug)
    cout<<"sheldon VortexFilamentSystem::addRing"<<endl;
    normalize(normal);
    Vec3 worldup (0,1,0);
    if (norm(normal - worldup) < 1e-5) worldup = Vec3(1,0,0);
    
    Vec3 u = cross(normal, worldup); normalize(u);
    Vec3 v = cross(normal, u); normalize(v);
    
    VortexRing ring(circulation*number*2*radius*sin(M_PI/number), true);
    
    ring.total_circulation=circulation*number*2*radius*sin(1*M_PI/number);
    ring.total_length=number*2*radius*sin(M_PI/number);
    for (int i=0; i<number; i++) {
        Real phi = (Real)i/(Real)number * M_PI * 2.0; //it seems that here this divides the ring into a number-edge polygon 
        Vec3 p = position + radius * (u*cos(phi) + v*sin(phi));
        
        int num = add(BasicParticleData(p)); //BasicParticleData generate one particle with the partilce's position and
        //a flag........add() executes the code mData.push_back(data); return mData.size()-1  
        ring.indices.push_back(num);
    }
    mSegments.push_back(ring);
}
void VortexFilamentSystem::revise_circulation()
{
    std::vector<bool> mark;
    int seg_size=segSize();
    for(int i=0;i<seg_size;i++)
    {
        VortexRing& ring=mSegments[i];
        int index_size=ring.indices.size();
        Real circum=0;
        for(int j=0;j<index_size;j++)
        {
            circum+=sqrt(normSquare(mData[ring.idx0(j)].pos-mData[ring.idx1(j)].pos));
        }
        if(ring.total_length/circum>3)
            mark[i]=true;
    }
    return;
}

    
} // namespace
/*void VortexFilamentSystem::split(Real cosine_threshold,Real dist_threshold){
    cout<<"Enter split function"<<endl;
    int seg_size=segSize();
    for(int i=0;i<seg_size;i++){
        VortexRing& r=mSegments[i];
        int size=r.size();
        for(int j=0;j<size;j++){
            bool mark=false;
            for(int k=j+3;k<size;k++){
                const Vec3 fir=mData[r.idx1(j)].pos-mData[r.idx0(j)].pos;
                const Vec3 mid1=(mData[r.idx1(j)].pos+mData[r.idx0(j)].pos)/2;
                const Vec3 sec=mData[r.idx1(k)].pos-mData[r.idx0(k)].pos;
                const Vec3 mid2=(mData[r.idx1(k)].pos+mData[r.idx0(k)].pos)/2;
                const Real cosine_res=dot(getNormalized(fir),getNormalized(sec));
                if(cosine_res>cosine_threshold)
                    continue;
                const Real dist1=normSquare(cross(mid1-mData[r.idx0(k)].pos,sec))/normSquare(sec); 
                const Real dist2=normSquare(cross(mid2-mData[r.idx0(j)].pos,fir))/normSquare(fir);
                Real dist=(sqrt(dist1)+sqrt(dist2))/2;
                if(normSquare(mid1-mid2)<dist_threshold*dist_threshold)
                   dist=sqrt(normSquare(mid1-mid2)); 
                if(dist>dist_threshold)
                    continue;
                if(((k-j)>5&&(size-k+j)>5)==false)
                    continue;
                int pos1=add(mid1);
                int pos2=add(mid2);
                int pos3=add(mid1);
                int pos4=add(mid2);
                vector<int> index1;
                for(int sub=j+1;sub<=k;sub++)
                    index1.push_back(r.indices[sub]);
                VortexRing new_ring;
                for(int sub=k+1;sub<=size+j;sub++)
                    new_ring.indices.push_back(r.indices[sub%size]);
                
                for(int sub=0;sub<=(k-j-1);sub++)
                    r.indices[sub]=index1[sub];
                r.indices[k-j]=pos2;
                r.indices[k-j+1]=pos1;
                r.indices.resize(k-j+2);
                new_ring.indices.push_back(pos3);
                new_ring.indices.push_back(pos4);
                mSegments.push_back(new_ring); 
                dirty.insert(mSegments.size());
                dirty.insert(i);
                
                mark=true;
                cout<<"cosine_res "<<cosine_res<<"  normSquare(sec) "<<normSquare(sec)<<"normSquare(fir)"<<normSquare(fir)<<"  "<<r.indices.size()<<"  "<<new_ring.indices.size()<<endl; 
                cout<<"successful split "<<i<<endl;
                break;

            }
            if(mark==true)
                break;
        }
        }
        cout<<"Exit split function"<<endl;
}*/
