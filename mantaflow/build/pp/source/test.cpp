




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/test.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include <cmath>

using namespace std;

namespace Manta {

//! Kernel: get component (not shifted)

struct GetComponent2 : public KernelBase {  GetComponent2 (const Grid<Vec3 >& _grid, int _dim) : KernelBase(&_grid, 0), parent((_grid).getParent()), m_grid(_grid), m_dim(_dim), ret(parent) { run(); } GetComponent2 (const GetComponent2& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_grid(o.m_grid), m_dim(o.m_dim), ret(parent) {} inline void op(int idx,const Grid<Vec3 >& grid, int dim, Grid<Real >& ret)  {
    ret[idx] = grid[idx][dim];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_grid, m_dim, ret); } operator Grid<Real>  () { return ret; } const Grid<Real >& getRet() const { return ret; } const Grid<Vec3 >& getArg0() { return m_grid; } typedef Grid<Vec3 > type0; int& getArg1() { return m_dim; } typedef int type1; FluidSolver* parent; const Grid<Vec3 >& m_grid; int m_dim; Grid<Real > ret;  }; ;

void testp(Grid<Vec3>& b, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Grid<Real> d(parent);
    b(20,20,20) = Vec3(21,22,23); 
    {
        cout <<"middle" << endl;        
        Grid<Real> a = GetComponent2(b,0);
        cout << a(20,20,20) << endl;        
        cout <<"middle" << endl;        
    }
    cout << "end" << endl;errMsg("f");
}PyObject* _plugin_testp (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "testp"); PyObject *_retval = NULL; { ArgLocker _lock; Grid<Vec3>& b = *_args.get< Grid<Vec3>* > (0,"b", &_lock); _retval = getPyNone();testp(b, parent);_args.check(); } pbFinalizePlugin(parent,"testp"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("testp",e.what()); return 0; } } 




struct ddtest : public KernelBase {  ddtest (const Grid<Real >& _v) : KernelBase(&_v, 0), parent((_v).getParent()), m_v(_v), sum(0) { run(); } ddtest (const ddtest& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_v(o.m_v), sum(0) {} inline void op(int idx,const Grid<Real >& v, double& sum)  {
    sum += v[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_v, sum); } operator double  () { return sum; } const double& getRet() const { return sum; } const Grid<Real >& getArg0() { return m_v; } typedef Grid<Real > type0; FluidSolver* parent; const Grid<Real >& m_v; double sum;  }; 



struct detest : public KernelBase {  detest (const Grid<Real >& _v) : KernelBase(&_v, 0), parent((_v).getParent()), m_v(_v), sum(0) { run(); } detest (const detest& o) : KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), parent(o.parent), m_v(o.m_v), sum(0) {} inline void op(int idx,const Grid<Real >& v, double& sum)  {
    if (sum < v[idx])
        sum = v[idx];
} void run() { const int _maxCells = maxCells; for (int idx=0; idx < _maxCells; idx++) op(idx, m_v, sum); } operator double  () { return sum; } const double& getRet() const { return sum; } const Grid<Real >& getArg0() { return m_v; } typedef Grid<Real > type0; FluidSolver* parent; const Grid<Real >& m_v; double sum;  }; 

void checkGrids(Grid<int>& flags1,Grid<int>& flags2,Grid<Real>& phi1,Grid<Real>& phi2,Grid<Vec3>& vel1,Grid<Vec3>& vel2, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    FOR_IJK(flags1) {
        assertMsg(flags1(i,j,k) == flags2(i,j,k), "flags mismatch");
        assertMsg(norm(vel1(i,j,k)-vel2(i,j,k)) < 1e-1, "vel mismatch");
        assertMsg( fabs(phi1(i,j,k)-phi2(i,j,k)) < 1e-4, "phi mismatch");
    }
}PyObject* _plugin_checkGrids (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "checkGrids"); PyObject *_retval = NULL; { ArgLocker _lock; Grid<int>& flags1 = *_args.get< Grid<int>* > (0,"flags1", &_lock); Grid<int>& flags2 = *_args.get< Grid<int>* > (1,"flags2", &_lock); Grid<Real>& phi1 = *_args.get< Grid<Real>* > (2,"phi1", &_lock); Grid<Real>& phi2 = *_args.get< Grid<Real>* > (3,"phi2", &_lock); Grid<Vec3>& vel1 = *_args.get< Grid<Vec3>* > (4,"vel1", &_lock); Grid<Vec3>& vel2 = *_args.get< Grid<Vec3>* > (5,"vel2", &_lock); _retval = getPyNone();checkGrids(flags1, flags2, phi1, phi2, vel1, vel2, parent);_args.check(); } pbFinalizePlugin(parent,"checkGrids"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("checkGrids",e.what()); return 0; } } 


struct myvec {
    myvec(int n) : x(n) { cout << "constructor" << endl; };
    myvec(const myvec& a) : x(a.x) { cout << "copy constructor" << endl; }
    myvec& operator=(const myvec& a) { x=a.x; cout << "copy operator" << endl; return *this;}
    int& operator[](int idx) { return x[idx]; }
    
    vector<int> x;
};


struct testy : public ParticleKernelBase {  testy (vector<int >& _a) : ParticleKernelBase((_a).size()), m_a(_a), vec(size) { run(); } testy (const testy& o) : ParticleKernelBase(o.size), m_a(o.m_a), vec(size) {} inline void op(int i,vector<int >& a, myvec& vec)  {
    vec[i] = a[i];
} void run() { const int _sz = size; for (int i=0; i < _sz; i++) op(i, m_a, vec); } operator myvec  () { return vec; } const myvec& getRet() const { return vec; } vector<int >& getArg0() { return m_a; } typedef vector<int > type0; vector<int >& m_a; myvec vec;  }; 

void kernelTest( FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    cout << "kernel test" << endl;
    vector<int> a(10);
    for (int i=0;i<10;i++) a[i]=i;
    
    //testy xx(a);
    myvec b = testy(a);
    for (int i=0;i<10;i++) cout << b[i] << endl;
    cout << "kernel end" << endl;
}PyObject* _plugin_kernelTest (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "kernelTest"); PyObject *_retval = NULL; { ArgLocker _lock; _retval = getPyNone();kernelTest(parent);_args.check(); } pbFinalizePlugin(parent,"kernelTest"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("kernelTest",e.what()); return 0; } } 

void getCurl(MACGrid& vel,Grid<Real>& vort,int comp, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    Grid<Vec3> velCenter(parent), curl(parent);
    
    GetCentered(velCenter, vel);
    CurlOp(velCenter, curl);
    GetComponent(curl, vort, comp);
}PyObject* _plugin_getCurl (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "getCurl"); PyObject *_retval = NULL; { ArgLocker _lock; MACGrid& vel = *_args.get< MACGrid* > (0,"vel", &_lock); Grid<Real>& vort = *_args.get< Grid<Real>* > (1,"vort", &_lock); int comp = _args.get< int > (2,"comp", &_lock); _retval = getPyNone();getCurl(vel, vort, comp, parent);_args.check(); } pbFinalizePlugin(parent,"getCurl"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("getCurl",e.what()); return 0; } } 

void setinflow(FlagGrid& flags,MACGrid& vel,LevelsetGrid& phi,Real h, FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY) {
    FOR_IJK(vel) {
        if (i<=2) {
            if (j < h*flags.getSizeY()) {
                vel(i,j,k).x = 1;            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 1;        
                    phi(i,j,k) = -1;
                }                
            } else {
                vel(i,j,k).x = 0;                            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 4;
                    phi(i,j,k) = 1;
                }
            }
        }
        else if (i>=flags.getSizeX()-2) {
            vel(i,j,k).x = 1;            
            /*if (j < 30-12) {
                vel(i,j,k).x = 1;            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 1;        
                    phi(i,j,k) = -1;
                }                
            } else {
                vel(i,j,k).x = 0;                            
                if (!flags.isObstacle(i,j,k)) { 
                    flags(i,j,k) = 4;
                    phi(i,j,k) = 1;
                }
            }*/
        }
    }
}PyObject* _plugin_setinflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); pbPreparePlugin(parent, "setinflow"); PyObject *_retval = NULL; { ArgLocker _lock; FlagGrid& flags = *_args.get< FlagGrid* > (0,"flags", &_lock); MACGrid& vel = *_args.get< MACGrid* > (1,"vel", &_lock); LevelsetGrid& phi = *_args.get< LevelsetGrid* > (2,"phi", &_lock); Real h = _args.get< Real > (3,"h", &_lock); _retval = getPyNone();setinflow(flags, vel, phi, h, parent);_args.check(); } pbFinalizePlugin(parent,"setinflow"); return (_retval) ? _retval : getPyNone(); } catch(std::exception& e) { pbSetError("setinflow",e.what()); return 0; } } 
    

} //namespace

