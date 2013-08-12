




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/shapes.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * shapes classes
 *
 ******************************************************************************/

#ifndef _SHAPES_H
#define _SHAPES_H

#include "pclass.h"
#include "vectorbase.h"
#include "levelset.h"

namespace Manta {

// forward declaration
class Mesh;
    
//! Base class for all shapes
class Shape : public PbClass {
public:
    enum GridType { TypeNone = 0, TypeBox = 1, TypeSphere = 2, TypeCylinder };
    
     Shape(FluidSolver* parent) ;
    
    //! Get the type of grid
    inline GridType getType() const { return mType; }
    
    //! Apply shape to flag grid, set inside cells to <value>
    void applyToGrid(GridBase* grid,FlagGrid* respectFlags=0) ;PyObject* _applyToGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Shape::applyToGrid"); { ArgLocker _lock; GridBase* grid = __args.get< GridBase* > (0,"grid", &_lock); FlagGrid* respectFlags = __args.getOpt< FlagGrid* > (1,"respectFlags", 0, &_lock); this->_args.copy(__args); _retval = getPyNone();applyToGrid(grid, respectFlags);this->_args.check(); } pbFinalizePlugin(this->mParent,"Shape::applyToGrid"); return _retval; } catch(std::exception& e) { pbSetError("Shape::applyToGrid",e.what()); return 0; } } 
    void applyToGridSmooth(GridBase* grid,Real sigma=1.0,Real shift=0,FlagGrid* respectFlags=0) ;PyObject* _applyToGridSmooth (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Shape::applyToGridSmooth"); { ArgLocker _lock; GridBase* grid = __args.get< GridBase* > (0,"grid", &_lock); Real sigma = __args.getOpt< Real > (1,"sigma", 1.0, &_lock); Real shift = __args.getOpt< Real > (2,"shift", 0, &_lock); FlagGrid* respectFlags = __args.getOpt< FlagGrid* > (3,"respectFlags", 0, &_lock); this->_args.copy(__args); _retval = getPyNone();applyToGridSmooth(grid, sigma, shift, respectFlags);this->_args.check(); } pbFinalizePlugin(this->mParent,"Shape::applyToGridSmooth"); return _retval; } catch(std::exception& e) { pbSetError("Shape::applyToGridSmooth",e.what()); return 0; } } 
    LevelsetGrid computeLevelset() ;PyObject* _computeLevelset (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Shape::computeLevelset"); { ArgLocker _lock; this->_args.copy(__args); _retval = d_toPy(computeLevelset() );this->_args.check(); } pbFinalizePlugin(this->mParent,"Shape::computeLevelset"); return _retval; } catch(std::exception& e) { pbSetError("Shape::computeLevelset",e.what()); return 0; } } 
    void collideMesh(Mesh& mesh) ;PyObject* _collideMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Shape::collideMesh"); { ArgLocker _lock; Mesh& mesh = *__args.get< Mesh* > (0,"mesh", &_lock); this->_args.copy(__args); _retval = getPyNone();collideMesh(mesh);this->_args.check(); } pbFinalizePlugin(this->mParent,"Shape::collideMesh"); return _retval; } catch(std::exception& e) { pbSetError("Shape::collideMesh",e.what()); return 0; } } 
    
    //! Inside test of the shape
    virtual bool isInside(const Vec3& pos) const;
    inline bool isInsideGrid(int i, int j, int k) const { return isInside(Vec3(i+0.5,j+0.5,k+0.5)); };
    
    virtual void generateMesh(Mesh* mesh) {} ;    
protected:
    virtual void generateLevelset(Grid<Real>& phi) {};    
    
    GridType mType;
protected:PbArgs _args;};;

//! Cylindrical shape
class NullShape : public Shape {    
public:
     NullShape(FluidSolver* parent) : Shape(parent) {}
    
    virtual bool isInside(const Vec3& pos) const { return false; }
    virtual void generateMesh(Mesh* mesh) {}
       
protected:
    virtual void generateLevelset(Grid<Real>& phi) { phi = 1000.0f; }
protected:PbArgs _args;};;

//! Box shape
class Box : public Shape {    
public:
     Box(FluidSolver* parent,Vec3 center = Vec3::Invalid,Vec3 p0 = Vec3::Invalid,Vec3 p1 = Vec3::Invalid,Vec3 size = Vec3::Invalid) ;
    
    inline Vec3 getCenter() const { return 0.5*(mP0+mP1); }
    inline Vec3 getSize() const { return mP1-mP0; }
    inline Vec3 getP0() const { return mP0; }
    inline Vec3 getP1() const { return mP1; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
        
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mP0, mP1;
protected:PbArgs _args;};;

//! Spherical shape
class Sphere : public Shape {    
public:
     Sphere(FluidSolver* parent,Vec3 center,Real radius,Vec3 scale=Vec3(1,1,1)) ;
    
    inline Vec3 getCenter() const { return mCenter; }
    inline Real getRadius() const { return mRadius; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
       
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mCenter, mScale;
    Real mRadius;
protected:PbArgs _args;};;

//! Cylindrical shape
class Cylinder : public Shape {    
public:
     Cylinder(FluidSolver* parent,Vec3 center,Real radius,Vec3 z) ;
    
    void setCenter(Vec3 center) { mCenter=center; }PyObject* _setCenter (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Cylinder::setCenter"); { ArgLocker _lock; Vec3 center = __args.get< Vec3 > (0,"center", &_lock); this->_args.copy(__args); _retval = getPyNone();setCenter(center);this->_args.check(); } pbFinalizePlugin(this->mParent,"Cylinder::setCenter"); return _retval; } catch(std::exception& e) { pbSetError("Cylinder::setCenter",e.what()); return 0; } } 
    void setRadius(Real r) { mRadius = r; }PyObject* _setRadius (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Cylinder::setRadius"); { ArgLocker _lock; Real r = __args.get< Real > (0,"r", &_lock); this->_args.copy(__args); _retval = getPyNone();setRadius(r);this->_args.check(); } pbFinalizePlugin(this->mParent,"Cylinder::setRadius"); return _retval; } catch(std::exception& e) { pbSetError("Cylinder::setRadius",e.what()); return 0; } } 
    void setZ(Vec3 z) { mZDir=z; mZ=normalize(mZDir); }PyObject* _setZ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Cylinder::setZ"); { ArgLocker _lock; Vec3 z = __args.get< Vec3 > (0,"z", &_lock); this->_args.copy(__args); _retval = getPyNone();setZ(z);this->_args.check(); } pbFinalizePlugin(this->mParent,"Cylinder::setZ"); return _retval; } catch(std::exception& e) { pbSetError("Cylinder::setZ",e.what()); return 0; } } 
    
    inline Vec3 getCenter() const { return mCenter; }
    inline Real getRadius() const { return mRadius; }
    inline Vec3 getZ() const { return mZ*mZDir; }
    virtual bool isInside(const Vec3& pos) const;
    virtual void generateMesh(Mesh* mesh);
        
protected:
    virtual void generateLevelset(Grid<Real>& phi);
    
    Vec3 mCenter, mZDir;
    Real mRadius, mZ;
protected:PbArgs _args;};;

    

} //namespace
#endif

