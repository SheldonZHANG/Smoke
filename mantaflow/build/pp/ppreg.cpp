




// This file is generated using the MantaFlow preprocessor merge stage (prep merge). Do not edit.




#include "pclass.h"
#include "fluidsolver.h"
#include "grid.h"
#include "mesh.h"
#include "particle.h"
#include "levelset.h"
#include "shapes.h"
#include "noisefield.h"
#include "vortexsheet.h"
#include "vortexfilament.h"
#include "flip.h"
#include "vortexpart.h"
#include "customctrl.h"

namespace Manta {
template<> FluidSolver* fromPy<FluidSolver*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("FluidSolver"))) throw Error("can't convert argument to type 'FluidSolver'"); return dynamic_cast<FluidSolver*>(pbo); }
template<> PyObject* toPy< FluidSolver >( FluidSolver& v) { if (v.getPyObject()) return v.getPyObject(); FluidSolver* co = new FluidSolver (v); return co->assignNewPyObject("FluidSolver"); }
template<> GridBase* fromPy<GridBase*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("GridBase"))) throw Error("can't convert argument to type 'GridBase'"); return dynamic_cast<GridBase*>(pbo); }
template<> PyObject* toPy< GridBase >( GridBase& v) { if (v.getPyObject()) return v.getPyObject(); GridBase* co = new GridBase (v); return co->assignNewPyObject("GridBase"); }
template<> Grid<int>* fromPy<Grid<int>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Grid<int>"))) throw Error("can't convert argument to type 'Grid<int>'"); return dynamic_cast<Grid<int>*>(pbo); }
template<> PyObject* toPy< Grid<int> >( Grid<int>& v) { if (v.getPyObject()) return v.getPyObject(); Grid<int>* co = new Grid<int> (v); return co->assignNewPyObject("Grid<int>"); }

template<> Grid<Real>* fromPy<Grid<Real>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Grid<Real>"))) throw Error("can't convert argument to type 'Grid<Real>'"); return dynamic_cast<Grid<Real>*>(pbo); }
template<> PyObject* toPy< Grid<Real> >( Grid<Real>& v) { if (v.getPyObject()) return v.getPyObject(); Grid<Real>* co = new Grid<Real> (v); return co->assignNewPyObject("Grid<Real>"); }

template<> Grid<Vec3>* fromPy<Grid<Vec3>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Grid<Vec3>"))) throw Error("can't convert argument to type 'Grid<Vec3>'"); return dynamic_cast<Grid<Vec3>*>(pbo); }
template<> PyObject* toPy< Grid<Vec3> >( Grid<Vec3>& v) { if (v.getPyObject()) return v.getPyObject(); Grid<Vec3>* co = new Grid<Vec3> (v); return co->assignNewPyObject("Grid<Vec3>"); }

template<> MACGrid* fromPy<MACGrid*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("MACGrid"))) throw Error("can't convert argument to type 'MACGrid'"); return dynamic_cast<MACGrid*>(pbo); }
template<> PyObject* toPy< MACGrid >( MACGrid& v) { if (v.getPyObject()) return v.getPyObject(); MACGrid* co = new MACGrid (v); return co->assignNewPyObject("MACGrid"); }
template<> FlagGrid* fromPy<FlagGrid*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("FlagGrid"))) throw Error("can't convert argument to type 'FlagGrid'"); return dynamic_cast<FlagGrid*>(pbo); }
template<> PyObject* toPy< FlagGrid >( FlagGrid& v) { if (v.getPyObject()) return v.getPyObject(); FlagGrid* co = new FlagGrid (v); return co->assignNewPyObject("FlagGrid"); }
template<> Mesh* fromPy<Mesh*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Mesh"))) throw Error("can't convert argument to type 'Mesh'"); return dynamic_cast<Mesh*>(pbo); }
template<> PyObject* toPy< Mesh >( Mesh& v) { if (v.getPyObject()) return v.getPyObject(); Mesh* co = new Mesh (v); return co->assignNewPyObject("Mesh"); }
template<> ParticleBase* fromPy<ParticleBase*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("ParticleBase"))) throw Error("can't convert argument to type 'ParticleBase'"); return dynamic_cast<ParticleBase*>(pbo); }
template<> PyObject* toPy< ParticleBase >( ParticleBase& v) { if (v.getPyObject()) return v.getPyObject(); ParticleBase* co = new ParticleBase (v); return co->assignNewPyObject("ParticleBase"); }
template<> ParticleSystem<BasicParticleData>* fromPy<ParticleSystem<BasicParticleData>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("ParticleSystem<BasicParticleData>"))) throw Error("can't convert argument to type 'ParticleSystem<BasicParticleData>'"); return dynamic_cast<ParticleSystem<BasicParticleData>*>(pbo); }
template<> PyObject* toPy< ParticleSystem<BasicParticleData> >( ParticleSystem<BasicParticleData>& v) { if (v.getPyObject()) return v.getPyObject(); ParticleSystem<BasicParticleData>* co = new ParticleSystem<BasicParticleData> (v); return co->assignNewPyObject("ParticleSystem<BasicParticleData>"); }

template<> TracerParticleSystem* fromPy<TracerParticleSystem*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("TracerParticleSystem"))) throw Error("can't convert argument to type 'TracerParticleSystem'"); return dynamic_cast<TracerParticleSystem*>(pbo); }
template<> PyObject* toPy< TracerParticleSystem >( TracerParticleSystem& v) { if (v.getPyObject()) return v.getPyObject(); TracerParticleSystem* co = new TracerParticleSystem (v); return co->assignNewPyObject("TracerParticleSystem"); }
template<> LevelsetGrid* fromPy<LevelsetGrid*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("LevelsetGrid"))) throw Error("can't convert argument to type 'LevelsetGrid'"); return dynamic_cast<LevelsetGrid*>(pbo); }
template<> PyObject* toPy< LevelsetGrid >( LevelsetGrid& v) { if (v.getPyObject()) return v.getPyObject(); LevelsetGrid* co = new LevelsetGrid (v); return co->assignNewPyObject("LevelsetGrid"); }
template<> Shape* fromPy<Shape*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Shape"))) throw Error("can't convert argument to type 'Shape'"); return dynamic_cast<Shape*>(pbo); }
template<> PyObject* toPy< Shape >( Shape& v) { if (v.getPyObject()) return v.getPyObject(); Shape* co = new Shape (v); return co->assignNewPyObject("Shape"); }
template<> NullShape* fromPy<NullShape*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("NullShape"))) throw Error("can't convert argument to type 'NullShape'"); return dynamic_cast<NullShape*>(pbo); }
template<> PyObject* toPy< NullShape >( NullShape& v) { if (v.getPyObject()) return v.getPyObject(); NullShape* co = new NullShape (v); return co->assignNewPyObject("NullShape"); }
template<> Box* fromPy<Box*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Box"))) throw Error("can't convert argument to type 'Box'"); return dynamic_cast<Box*>(pbo); }
template<> PyObject* toPy< Box >( Box& v) { if (v.getPyObject()) return v.getPyObject(); Box* co = new Box (v); return co->assignNewPyObject("Box"); }
template<> Sphere* fromPy<Sphere*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Sphere"))) throw Error("can't convert argument to type 'Sphere'"); return dynamic_cast<Sphere*>(pbo); }
template<> PyObject* toPy< Sphere >( Sphere& v) { if (v.getPyObject()) return v.getPyObject(); Sphere* co = new Sphere (v); return co->assignNewPyObject("Sphere"); }
template<> Cylinder* fromPy<Cylinder*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Cylinder"))) throw Error("can't convert argument to type 'Cylinder'"); return dynamic_cast<Cylinder*>(pbo); }
template<> PyObject* toPy< Cylinder >( Cylinder& v) { if (v.getPyObject()) return v.getPyObject(); Cylinder* co = new Cylinder (v); return co->assignNewPyObject("Cylinder"); }
template<> WaveletNoiseField* fromPy<WaveletNoiseField*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("WaveletNoiseField"))) throw Error("can't convert argument to type 'WaveletNoiseField'"); return dynamic_cast<WaveletNoiseField*>(pbo); }
template<> PyObject* toPy< WaveletNoiseField >( WaveletNoiseField& v) { if (v.getPyObject()) return v.getPyObject(); WaveletNoiseField* co = new WaveletNoiseField (v); return co->assignNewPyObject("WaveletNoiseField"); }
template<> VortexSheetMesh* fromPy<VortexSheetMesh*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("VortexSheetMesh"))) throw Error("can't convert argument to type 'VortexSheetMesh'"); return dynamic_cast<VortexSheetMesh*>(pbo); }
template<> PyObject* toPy< VortexSheetMesh >( VortexSheetMesh& v) { if (v.getPyObject()) return v.getPyObject(); VortexSheetMesh* co = new VortexSheetMesh (v); return co->assignNewPyObject("VortexSheetMesh"); }
template<> ConnectedParticleSystem<BasicParticleData,VortexRing>* fromPy<ConnectedParticleSystem<BasicParticleData,VortexRing>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("ConnectedParticleSystem<BasicParticleData,VortexRing>"))) throw Error("can't convert argument to type 'ConnectedParticleSystem<BasicParticleData,VortexRing>'"); return dynamic_cast<ConnectedParticleSystem<BasicParticleData,VortexRing>*>(pbo); }
template<> PyObject* toPy< ConnectedParticleSystem<BasicParticleData,VortexRing> >( ConnectedParticleSystem<BasicParticleData,VortexRing>& v) { if (v.getPyObject()) return v.getPyObject(); ConnectedParticleSystem<BasicParticleData,VortexRing>* co = new ConnectedParticleSystem<BasicParticleData,VortexRing> (v); return co->assignNewPyObject("ConnectedParticleSystem<BasicParticleData,VortexRing>"); }

template<> VortexFilamentSystem* fromPy<VortexFilamentSystem*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("VortexFilamentSystem"))) throw Error("can't convert argument to type 'VortexFilamentSystem'"); return dynamic_cast<VortexFilamentSystem*>(pbo); }
template<> PyObject* toPy< VortexFilamentSystem >( VortexFilamentSystem& v) { if (v.getPyObject()) return v.getPyObject(); VortexFilamentSystem* co = new VortexFilamentSystem (v); return co->assignNewPyObject("VortexFilamentSystem"); }
template<> ParticleSystem<FlipData>* fromPy<ParticleSystem<FlipData>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("ParticleSystem<FlipData>"))) throw Error("can't convert argument to type 'ParticleSystem<FlipData>'"); return dynamic_cast<ParticleSystem<FlipData>*>(pbo); }
template<> PyObject* toPy< ParticleSystem<FlipData> >( ParticleSystem<FlipData>& v) { if (v.getPyObject()) return v.getPyObject(); ParticleSystem<FlipData>* co = new ParticleSystem<FlipData> (v); return co->assignNewPyObject("ParticleSystem<FlipData>"); }

template<> FlipSystem* fromPy<FlipSystem*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("FlipSystem"))) throw Error("can't convert argument to type 'FlipSystem'"); return dynamic_cast<FlipSystem*>(pbo); }
template<> PyObject* toPy< FlipSystem >( FlipSystem& v) { if (v.getPyObject()) return v.getPyObject(); FlipSystem* co = new FlipSystem (v); return co->assignNewPyObject("FlipSystem"); }
template<> ParticleSystem<VortexParticleData>* fromPy<ParticleSystem<VortexParticleData>*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("ParticleSystem<VortexParticleData>"))) throw Error("can't convert argument to type 'ParticleSystem<VortexParticleData>'"); return dynamic_cast<ParticleSystem<VortexParticleData>*>(pbo); }
template<> PyObject* toPy< ParticleSystem<VortexParticleData> >( ParticleSystem<VortexParticleData>& v) { if (v.getPyObject()) return v.getPyObject(); ParticleSystem<VortexParticleData>* co = new ParticleSystem<VortexParticleData> (v); return co->assignNewPyObject("ParticleSystem<VortexParticleData>"); }

template<> VortexParticleSystem* fromPy<VortexParticleSystem*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("VortexParticleSystem"))) throw Error("can't convert argument to type 'VortexParticleSystem'"); return dynamic_cast<VortexParticleSystem*>(pbo); }
template<> PyObject* toPy< VortexParticleSystem >( VortexParticleSystem& v) { if (v.getPyObject()) return v.getPyObject(); VortexParticleSystem* co = new VortexParticleSystem (v); return co->assignNewPyObject("VortexParticleSystem"); }
template<> CustomControl* fromPy<CustomControl*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("CustomControl"))) throw Error("can't convert argument to type 'CustomControl'"); return dynamic_cast<CustomControl*>(pbo); }
template<> PyObject* toPy< CustomControl >( CustomControl& v) { if (v.getPyObject()) return v.getPyObject(); CustomControl* co = new CustomControl (v); return co->assignNewPyObject("CustomControl"); }
template<> CustomSlider* fromPy<CustomSlider*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("CustomSlider"))) throw Error("can't convert argument to type 'CustomSlider'"); return dynamic_cast<CustomSlider*>(pbo); }
template<> PyObject* toPy< CustomSlider >( CustomSlider& v) { if (v.getPyObject()) return v.getPyObject(); CustomSlider* co = new CustomSlider (v); return co->assignNewPyObject("CustomSlider"); }
template<> Gui* fromPy<Gui*>(PyObject* obj) { if (PbClass::isNullRef(obj)) return 0; PbClass* pbo = PbClass::fromPyObject(obj); if (!pbo || !(pbo->canConvertTo("Gui"))) throw Error("can't convert argument to type 'Gui'"); return dynamic_cast<Gui*>(pbo); }
template<> PyObject* toPy< Gui >( Gui& v) { if (v.getPyObject()) return v.getPyObject(); Gui* co = new Gui (v); return co->assignNewPyObject("Gui"); }
extern PyObject* _plugin_solvePressure (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_densityInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_advectSemiLagrange (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_addGravity (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_addBuoyancy (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_setWallBcs (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_setLiquidBcs (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_vorticityConfinement (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonComputeProduction (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonSources (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonSourcesMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonInit (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonInitMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_KEpsilonGradientDiffusion (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_smoothMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_subdivideMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_killSmallComponents (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_markAsFixed (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_texcoordInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_meshSmokeInflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_vorticitySource (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_smoothVorticity (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_VPseedK41 (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_VICintegration (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_densityFromLevelset (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_testp (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_checkGrids (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_kernelTest (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_getCurl (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
extern PyObject* _plugin_setinflow (PyObject* _self, PyObject* _linargs, PyObject* _kwds);
int _FluidSolver_FluidSolver (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "FluidSolver::FluidSolver"); { ArgLocker _lock; Vec3i gridSize = _args.get< Vec3i > (0,"gridSize", &_lock); int dim = _args.getOpt< int > (1,"dim", 3, &_lock); int is_debug = _args.getOpt< int > (2,"is_debug", 1, &_lock); obj = new FluidSolver(gridSize, dim, is_debug);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"FluidSolver::FluidSolver"); return 0; } catch(std::exception& e) { pbSetError("FluidSolver::FluidSolver",e.what()); return -1; } } 
PyObject* _FluidSolver_getGridSize (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FluidSolver*>(_self)->_getGridSize(_self, _linargs, _kwds); }
PyObject* _FluidSolver_printTimings (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FluidSolver*>(_self)->_printTimings(_self, _linargs, _kwds); }
PyObject* _FluidSolver_saveMeanTimings (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FluidSolver*>(_self)->_saveMeanTimings(_self, _linargs, _kwds); }
PyObject* _FluidSolver_step (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FluidSolver*>(_self)->_step(_self, _linargs, _kwds); }
PyObject* _FluidSolver_create (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FluidSolver*>(_self)->_create(_self, _linargs, _kwds); }
PyObject* _get_FluidSolver_mDt(PyObject* self, void* cl) { return d_toPy(fromPy<FluidSolver*>(self)->mDt); }
int _set_FluidSolver_mDt(PyObject* self, PyObject* val, void* cl) { fromPy<FluidSolver*>(self)->mDt=fromPy<Real >(val); return 0;}
int _GridBase_GridBase (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "GridBase::GridBase"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new GridBase(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"GridBase::GridBase"); return 0; } catch(std::exception& e) { pbSetError("GridBase::GridBase",e.what()); return -1; } } 
int _GridIntGrid_GridIntGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Grid::Grid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); bool show = _args.getOpt< bool > (1,"show", true, &_lock); obj = new Grid<int> (parent, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Grid::Grid"); return 0; } catch(std::exception& e) { pbSetError("Grid::Grid",e.what()); return -1; } } 
PyObject* _GridIntGrid_save (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<int>*>(PbClass::fromPyObject(_self))->_save(_self, _linargs, _kwds); }
PyObject* _GridIntGrid_load (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<int>*>(PbClass::fromPyObject(_self))->_load(_self, _linargs, _kwds); }
PyObject* _GridIntGrid_add (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<int>*>(PbClass::fromPyObject(_self))->_add(_self, _linargs, _kwds); }
PyObject* _GridIntGrid_scale (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<int>*>(PbClass::fromPyObject(_self))->_scale(_self, _linargs, _kwds); }
PyObject* _GridIntGrid_copyFrom (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<int>*>(PbClass::fromPyObject(_self))->_copyFrom(_self, _linargs, _kwds); }
int _GridRealGrid_GridRealGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Grid::Grid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); bool show = _args.getOpt< bool > (1,"show", true, &_lock); obj = new Grid<Real> (parent, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Grid::Grid"); return 0; } catch(std::exception& e) { pbSetError("Grid::Grid",e.what()); return -1; } } 
PyObject* _GridRealGrid_save (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Real>*>(PbClass::fromPyObject(_self))->_save(_self, _linargs, _kwds); }
PyObject* _GridRealGrid_load (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Real>*>(PbClass::fromPyObject(_self))->_load(_self, _linargs, _kwds); }
PyObject* _GridRealGrid_add (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Real>*>(PbClass::fromPyObject(_self))->_add(_self, _linargs, _kwds); }
PyObject* _GridRealGrid_scale (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Real>*>(PbClass::fromPyObject(_self))->_scale(_self, _linargs, _kwds); }
PyObject* _GridRealGrid_copyFrom (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Real>*>(PbClass::fromPyObject(_self))->_copyFrom(_self, _linargs, _kwds); }
int _GridVecGrid_GridVecGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Grid::Grid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); bool show = _args.getOpt< bool > (1,"show", true, &_lock); obj = new Grid<Vec3> (parent, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Grid::Grid"); return 0; } catch(std::exception& e) { pbSetError("Grid::Grid",e.what()); return -1; } } 
PyObject* _GridVecGrid_save (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Vec3>*>(PbClass::fromPyObject(_self))->_save(_self, _linargs, _kwds); }
PyObject* _GridVecGrid_load (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Vec3>*>(PbClass::fromPyObject(_self))->_load(_self, _linargs, _kwds); }
PyObject* _GridVecGrid_add (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Vec3>*>(PbClass::fromPyObject(_self))->_add(_self, _linargs, _kwds); }
PyObject* _GridVecGrid_scale (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Vec3>*>(PbClass::fromPyObject(_self))->_scale(_self, _linargs, _kwds); }
PyObject* _GridVecGrid_copyFrom (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<Grid<Vec3>*>(PbClass::fromPyObject(_self))->_copyFrom(_self, _linargs, _kwds); }
int _MACGrid_MACGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "MACGrid::MACGrid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); bool show = _args.getOpt< bool > (1,"show", true, &_lock); obj = new MACGrid(parent, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"MACGrid::MACGrid"); return 0; } catch(std::exception& e) { pbSetError("MACGrid::MACGrid",e.what()); return -1; } } 
int _FlagGrid_FlagGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "FlagGrid::FlagGrid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); int dim = _args.getOpt< int > (1,"dim", 3, &_lock); bool show = _args.getOpt< bool > (2,"show", true, &_lock); obj = new FlagGrid(parent, dim, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"FlagGrid::FlagGrid"); return 0; } catch(std::exception& e) { pbSetError("FlagGrid::FlagGrid",e.what()); return -1; } } 
PyObject* _FlagGrid_initDomain (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlagGrid*>(_self)->_initDomain(_self, _linargs, _kwds); }
PyObject* _FlagGrid_initBoundaries (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlagGrid*>(_self)->_initBoundaries(_self, _linargs, _kwds); }
PyObject* _FlagGrid_updateFromLevelset (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlagGrid*>(_self)->_updateFromLevelset(_self, _linargs, _kwds); }
PyObject* _FlagGrid_fillGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlagGrid*>(_self)->_fillGrid(_self, _linargs, _kwds); }
int _Mesh_Mesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Mesh::Mesh"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new Mesh(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Mesh::Mesh"); return 0; } catch(std::exception& e) { pbSetError("Mesh::Mesh",e.what()); return -1; } } 
PyObject* _Mesh_load (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_load(_self, _linargs, _kwds); }
PyObject* _Mesh_fromShape (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_fromShape(_self, _linargs, _kwds); }
PyObject* _Mesh_save (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_save(_self, _linargs, _kwds); }
PyObject* _Mesh_advectInGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_advectInGrid(_self, _linargs, _kwds); }
PyObject* _Mesh_scale (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_scale(_self, _linargs, _kwds); }
PyObject* _Mesh_offset (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Mesh*>(_self)->_offset(_self, _linargs, _kwds); }
int _ParticleBase_ParticleBase (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleBase::ParticleBase"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new ParticleBase(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleBase::ParticleBase"); return 0; } catch(std::exception& e) { pbSetError("ParticleBase::ParticleBase",e.what()); return -1; } } 
int _ParticleSystem_ParticleSystem_BasicParticleData_ParticleSystem_ParticleSystem_BasicParticleData (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleSystem::ParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new ParticleSystem<BasicParticleData> (parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleSystem::ParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("ParticleSystem::ParticleSystem",e.what()); return -1; } } 
PyObject* _ParticleSystem_ParticleSystem_BasicParticleData_size (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<BasicParticleData>*>(PbClass::fromPyObject(_self))->_size(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_BasicParticleData_setPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<BasicParticleData>*>(PbClass::fromPyObject(_self))->_setPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_BasicParticleData_getPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<BasicParticleData>*>(PbClass::fromPyObject(_self))->_getPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_BasicParticleData_advectInGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<BasicParticleData>*>(PbClass::fromPyObject(_self))->_advectInGrid(_self, _linargs, _kwds); }
int _TracerParticleSystem_TracerParticleSystem (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "TracerParticleSystem::TracerParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new TracerParticleSystem(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"TracerParticleSystem::TracerParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("TracerParticleSystem::TracerParticleSystem",e.what()); return -1; } } 
PyObject* _TracerParticleSystem_addParticle (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<TracerParticleSystem*>(_self)->_addParticle(_self, _linargs, _kwds); }
PyObject* _TracerParticleSystem_add_Particle (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<TracerParticleSystem*>(_self)->_add_Particle(_self, _linargs, _kwds); }
int _LevelsetGrid_LevelsetGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "LevelsetGrid::LevelsetGrid"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); bool show = _args.getOpt< bool > (1,"show", true, &_lock); obj = new LevelsetGrid(parent, show);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"LevelsetGrid::LevelsetGrid"); return 0; } catch(std::exception& e) { pbSetError("LevelsetGrid::LevelsetGrid",e.what()); return -1; } } 
PyObject* _LevelsetGrid_reinitMarching (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<LevelsetGrid*>(_self)->_reinitMarching(_self, _linargs, _kwds); }
PyObject* _LevelsetGrid_createMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<LevelsetGrid*>(_self)->_createMesh(_self, _linargs, _kwds); }
PyObject* _LevelsetGrid_join (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<LevelsetGrid*>(_self)->_join(_self, _linargs, _kwds); }
PyObject* _LevelsetGrid_initFromFlags (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<LevelsetGrid*>(_self)->_initFromFlags(_self, _linargs, _kwds); }
int _Shape_Shape (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Shape::Shape"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new Shape(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Shape::Shape"); return 0; } catch(std::exception& e) { pbSetError("Shape::Shape",e.what()); return -1; } } 
PyObject* _Shape_applyToGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Shape*>(_self)->_applyToGrid(_self, _linargs, _kwds); }
PyObject* _Shape_applyToGridSmooth (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Shape*>(_self)->_applyToGridSmooth(_self, _linargs, _kwds); }
PyObject* _Shape_computeLevelset (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Shape*>(_self)->_computeLevelset(_self, _linargs, _kwds); }
PyObject* _Shape_collideMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Shape*>(_self)->_collideMesh(_self, _linargs, _kwds); }
int _NullShape_NullShape (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "NullShape::NullShape"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new NullShape(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"NullShape::NullShape"); return 0; } catch(std::exception& e) { pbSetError("NullShape::NullShape",e.what()); return -1; } } 
int _Box_Box (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Box::Box"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); Vec3 center = _args.getOpt< Vec3 > (1,"center", Vec3::Invalid, &_lock); Vec3 p0 = _args.getOpt< Vec3 > (2,"p0", Vec3::Invalid, &_lock); Vec3 p1 = _args.getOpt< Vec3 > (3,"p1", Vec3::Invalid, &_lock); Vec3 size = _args.getOpt< Vec3 > (4,"size", Vec3::Invalid, &_lock); obj = new Box(parent, center, p0, p1, size);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Box::Box"); return 0; } catch(std::exception& e) { pbSetError("Box::Box",e.what()); return -1; } } 
int _Sphere_Sphere (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Sphere::Sphere"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); Vec3 center = _args.get< Vec3 > (1,"center", &_lock); Real radius = _args.get< Real > (2,"radius", &_lock); Vec3 scale = _args.getOpt< Vec3 > (3,"scale", Vec3(1,1,1), &_lock); obj = new Sphere(parent, center, radius, scale);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Sphere::Sphere"); return 0; } catch(std::exception& e) { pbSetError("Sphere::Sphere",e.what()); return -1; } } 
int _Cylinder_Cylinder (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Cylinder::Cylinder"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); Vec3 center = _args.get< Vec3 > (1,"center", &_lock); Real radius = _args.get< Real > (2,"radius", &_lock); Vec3 z = _args.get< Vec3 > (3,"z", &_lock); obj = new Cylinder(parent, center, radius, z);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Cylinder::Cylinder"); return 0; } catch(std::exception& e) { pbSetError("Cylinder::Cylinder",e.what()); return -1; } } 
PyObject* _Cylinder_setCenter (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Cylinder*>(_self)->_setCenter(_self, _linargs, _kwds); }
PyObject* _Cylinder_setRadius (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Cylinder*>(_self)->_setRadius(_self, _linargs, _kwds); }
PyObject* _Cylinder_setZ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Cylinder*>(_self)->_setZ(_self, _linargs, _kwds); }
int _WaveletNoiseField_WaveletNoiseField (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "WaveletNoiseField::WaveletNoiseField"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new WaveletNoiseField(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"WaveletNoiseField::WaveletNoiseField"); return 0; } catch(std::exception& e) { pbSetError("WaveletNoiseField::WaveletNoiseField",e.what()); return -1; } } 
PyObject* _get_WaveletNoiseField_mPosOffset(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mPosOffset); }
int _set_WaveletNoiseField_mPosOffset(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mPosOffset=fromPy<Vec3 >(val); return 0;}
PyObject* _get_WaveletNoiseField_mPosScale(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mPosScale); }
int _set_WaveletNoiseField_mPosScale(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mPosScale=fromPy<Vec3 >(val); return 0;}
PyObject* _get_WaveletNoiseField_mValOffset(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mValOffset); }
int _set_WaveletNoiseField_mValOffset(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mValOffset=fromPy<Real >(val); return 0;}
PyObject* _get_WaveletNoiseField_mValScale(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mValScale); }
int _set_WaveletNoiseField_mValScale(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mValScale=fromPy<Real >(val); return 0;}
PyObject* _get_WaveletNoiseField_mClamp(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mClamp); }
int _set_WaveletNoiseField_mClamp(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mClamp=fromPy<bool >(val); return 0;}
PyObject* _get_WaveletNoiseField_mClampNeg(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mClampNeg); }
int _set_WaveletNoiseField_mClampNeg(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mClampNeg=fromPy<Real >(val); return 0;}
PyObject* _get_WaveletNoiseField_mClampPos(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mClampPos); }
int _set_WaveletNoiseField_mClampPos(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mClampPos=fromPy<Real >(val); return 0;}
PyObject* _get_WaveletNoiseField_mTimeAnim(PyObject* self, void* cl) { return d_toPy(fromPy<WaveletNoiseField*>(self)->mTimeAnim); }
int _set_WaveletNoiseField_mTimeAnim(PyObject* self, PyObject* val, void* cl) { fromPy<WaveletNoiseField*>(self)->mTimeAnim=fromPy<Real >(val); return 0;}
int _VortexSheetMesh_VortexSheetMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "VortexSheetMesh::VortexSheetMesh"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new VortexSheetMesh(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"VortexSheetMesh::VortexSheetMesh"); return 0; } catch(std::exception& e) { pbSetError("VortexSheetMesh::VortexSheetMesh",e.what()); return -1; } } 
PyObject* _VortexSheetMesh_calcCirculation (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexSheetMesh*>(_self)->_calcCirculation(_self, _linargs, _kwds); }
PyObject* _VortexSheetMesh_calcVorticity (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexSheetMesh*>(_self)->_calcVorticity(_self, _linargs, _kwds); }
PyObject* _VortexSheetMesh_reinitTexCoords (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexSheetMesh*>(_self)->_reinitTexCoords(_self, _linargs, _kwds); }
int _ConnectedParticleSystem_ConnectedParticleSystem_BasicParticleData_VortexRing_ConnectedParticleSystem_ConnectedParticleSystem_BasicParticleData_VortexRing (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ConnectedParticleSystem::ConnectedParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new ConnectedParticleSystem<BasicParticleData,VortexRing> (parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"ConnectedParticleSystem::ConnectedParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("ConnectedParticleSystem::ConnectedParticleSystem",e.what()); return -1; } } 
int _VortexFilamentSystem_VortexFilamentSystem (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "VortexFilamentSystem::VortexFilamentSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new VortexFilamentSystem(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"VortexFilamentSystem::VortexFilamentSystem"); return 0; } catch(std::exception& e) { pbSetError("VortexFilamentSystem::VortexFilamentSystem",e.what()); return -1; } } 
PyObject* _VortexFilamentSystem_advectSelf (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_advectSelf(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_advectParticles (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_advectParticles(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_advectMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_advectMesh(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_doublyDiscreteUpdate (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_doublyDiscreteUpdate(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_remesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_remesh(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_split_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_split_ring(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_reconnect_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_reconnect_ring(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_merge_adj_edge (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_merge_adj_edge(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_divide_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_divide_ring(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_reset_dirty (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_reset_dirty(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_Debug_fun (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_Debug_fun(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_revise_circulation (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_revise_circulation(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_Decimate_ring (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_Decimate_ring(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_Direct_motion (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_Direct_motion(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_addRing (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_addRing(_self, _linargs, _kwds); }
PyObject* _VortexFilamentSystem_addLine (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexFilamentSystem*>(_self)->_addLine(_self, _linargs, _kwds); }
int _ParticleSystem_ParticleSystem_FlipData_ParticleSystem_ParticleSystem_FlipData (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleSystem::ParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new ParticleSystem<FlipData> (parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleSystem::ParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("ParticleSystem::ParticleSystem",e.what()); return -1; } } 
PyObject* _ParticleSystem_ParticleSystem_FlipData_size (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<FlipData>*>(PbClass::fromPyObject(_self))->_size(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_FlipData_setPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<FlipData>*>(PbClass::fromPyObject(_self))->_setPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_FlipData_getPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<FlipData>*>(PbClass::fromPyObject(_self))->_getPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_FlipData_advectInGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<FlipData>*>(PbClass::fromPyObject(_self))->_advectInGrid(_self, _linargs, _kwds); }
int _FlipSystem_FlipSystem (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "FlipSystem::FlipSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new FlipSystem(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"FlipSystem::FlipSystem"); return 0; } catch(std::exception& e) { pbSetError("FlipSystem::FlipSystem",e.what()); return -1; } } 
PyObject* _FlipSystem_velocitiesFromGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlipSystem*>(_self)->_velocitiesFromGrid(_self, _linargs, _kwds); }
PyObject* _FlipSystem_velocitiesToGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlipSystem*>(_self)->_velocitiesToGrid(_self, _linargs, _kwds); }
PyObject* _FlipSystem_adjustNumber (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlipSystem*>(_self)->_adjustNumber(_self, _linargs, _kwds); }
PyObject* _FlipSystem_markFluidCells (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<FlipSystem*>(_self)->_markFluidCells(_self, _linargs, _kwds); }
int _ParticleSystem_ParticleSystem_VortexParticleData_ParticleSystem_ParticleSystem_VortexParticleData (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "ParticleSystem::ParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new ParticleSystem<VortexParticleData> (parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"ParticleSystem::ParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("ParticleSystem::ParticleSystem",e.what()); return -1; } } 
PyObject* _ParticleSystem_ParticleSystem_VortexParticleData_size (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<VortexParticleData>*>(PbClass::fromPyObject(_self))->_size(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_VortexParticleData_setPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<VortexParticleData>*>(PbClass::fromPyObject(_self))->_setPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_VortexParticleData_getPos (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<VortexParticleData>*>(PbClass::fromPyObject(_self))->_getPos(_self, _linargs, _kwds); }
PyObject* _ParticleSystem_ParticleSystem_VortexParticleData_advectInGrid (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return dynamic_cast<ParticleSystem<VortexParticleData>*>(PbClass::fromPyObject(_self))->_advectInGrid(_self, _linargs, _kwds); }
int _VortexParticleSystem_VortexParticleSystem (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "VortexParticleSystem::VortexParticleSystem"); { ArgLocker _lock; FluidSolver* parent = _args.get< FluidSolver* > (0,"parent", &_lock); obj = new VortexParticleSystem(parent);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"VortexParticleSystem::VortexParticleSystem"); return 0; } catch(std::exception& e) { pbSetError("VortexParticleSystem::VortexParticleSystem",e.what()); return -1; } } 
PyObject* _VortexParticleSystem_advectSelf (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexParticleSystem*>(_self)->_advectSelf(_self, _linargs, _kwds); }
PyObject* _VortexParticleSystem_applyToMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<VortexParticleSystem*>(_self)->_applyToMesh(_self, _linargs, _kwds); }
int _CustomControl_CustomControl (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "CustomControl::CustomControl"); { ArgLocker _lock; obj = new CustomControl();std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"CustomControl::CustomControl"); return 0; } catch(std::exception& e) { pbSetError("CustomControl::CustomControl",e.what()); return -1; } } 
int _CustomSlider_CustomSlider (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "CustomSlider::CustomSlider"); { ArgLocker _lock; std::string text = _args.get< std::string > (0,"text", &_lock); float val = _args.get< float > (1,"val", &_lock); float min = _args.get< float > (2,"min", &_lock); float max = _args.get< float > (3,"max", &_lock); obj = new CustomSlider(text, val, min, max);std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"CustomSlider::CustomSlider"); return 0; } catch(std::exception& e) { pbSetError("CustomSlider::CustomSlider",e.what()); return -1; } } 
PyObject* _CustomSlider_get (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<CustomSlider*>(_self)->_get(_self, _linargs, _kwds); }
PyObject* _CustomSlider_set (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<CustomSlider*>(_self)->_set(_self, _linargs, _kwds); }
int _Gui_Gui (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { PbClass* obj = PbClass::fromPyObject(_self); if (obj) delete obj; try { PbArgs _args(_linargs, _kwds); pbPreparePlugin(0, "Gui::Gui"); { ArgLocker _lock; obj = new Gui();std::string _name = _args.getOpt<std::string>("name",""); obj->setPyObject(_self); if (!_name.empty()) obj->setName(_name); _args.check(); } pbFinalizePlugin(obj->getParent(),"Gui::Gui"); return 0; } catch(std::exception& e) { pbSetError("Gui::Gui",e.what()); return -1; } } 
PyObject* _Gui_setBackgroundMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_setBackgroundMesh(_self, _linargs, _kwds); }
PyObject* _Gui_show (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_show(_self, _linargs, _kwds); }
PyObject* _Gui_update (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_update(_self, _linargs, _kwds); }
PyObject* _Gui_pause (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_pause(_self, _linargs, _kwds); }
PyObject* _Gui_addControl (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_addControl(_self, _linargs, _kwds); }
PyObject* _Gui_screenshot (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { return fromPy<Gui*>(_self)->_screenshot(_self, _linargs, _kwds); }

struct _pbRegister {
	_pbRegister(int dummy) {
		PbWrapperRegistry::instance().addMethod("", "solvePressure", _plugin_solvePressure);
		PbWrapperRegistry::instance().addMethod("", "densityInflow", _plugin_densityInflow);
		PbWrapperRegistry::instance().addMethod("", "advectSemiLagrange", _plugin_advectSemiLagrange);
		PbWrapperRegistry::instance().addMethod("", "addGravity", _plugin_addGravity);
		PbWrapperRegistry::instance().addMethod("", "addBuoyancy", _plugin_addBuoyancy);
		PbWrapperRegistry::instance().addMethod("", "setWallBcs", _plugin_setWallBcs);
		PbWrapperRegistry::instance().addMethod("", "setLiquidBcs", _plugin_setLiquidBcs);
		PbWrapperRegistry::instance().addMethod("", "vorticityConfinement", _plugin_vorticityConfinement);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonComputeProduction", _plugin_KEpsilonComputeProduction);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonSources", _plugin_KEpsilonSources);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonSourcesMesh", _plugin_KEpsilonSourcesMesh);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonInit", _plugin_KEpsilonInit);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonInitMesh", _plugin_KEpsilonInitMesh);
		PbWrapperRegistry::instance().addMethod("", "KEpsilonGradientDiffusion", _plugin_KEpsilonGradientDiffusion);
		PbWrapperRegistry::instance().addMethod("", "smoothMesh", _plugin_smoothMesh);
		PbWrapperRegistry::instance().addMethod("", "subdivideMesh", _plugin_subdivideMesh);
		PbWrapperRegistry::instance().addMethod("", "killSmallComponents", _plugin_killSmallComponents);
		PbWrapperRegistry::instance().addMethod("", "markAsFixed", _plugin_markAsFixed);
		PbWrapperRegistry::instance().addMethod("", "texcoordInflow", _plugin_texcoordInflow);
		PbWrapperRegistry::instance().addMethod("", "meshSmokeInflow", _plugin_meshSmokeInflow);
		PbWrapperRegistry::instance().addMethod("", "vorticitySource", _plugin_vorticitySource);
		PbWrapperRegistry::instance().addMethod("", "smoothVorticity", _plugin_smoothVorticity);
		PbWrapperRegistry::instance().addMethod("", "VPseedK41", _plugin_VPseedK41);
		PbWrapperRegistry::instance().addMethod("", "VICintegration", _plugin_VICintegration);
		PbWrapperRegistry::instance().addMethod("", "densityFromLevelset", _plugin_densityFromLevelset);
		PbWrapperRegistry::instance().addPythonCode("/home/mzhang/mantaflow/source/python/defines.py", "################################################################################\n#\n# MantaFlow fluid solver framework\n# Copyright 2011 Tobias Pfaff, Nils Thuerey \n#\n# This program is free software, distributed under the terms of the\n# GNU General Public License (GPL) \n# http://www.gnu.org/licenses\n#\n# Defines some constants for use in python subprograms\n#\n#################################################################################\n\n# grid flags\nFlagFluid = 1\nFlagObstacle = 2\nFlagEmpty = 4\nFlagStick = 128\n\n# integration mode\nIntEuler = 0\nIntRK2 = 1\nIntRK4 = 2\n\n\n\n");
		PbWrapperRegistry::instance().addMethod("", "testp", _plugin_testp);
		PbWrapperRegistry::instance().addMethod("", "checkGrids", _plugin_checkGrids);
		PbWrapperRegistry::instance().addMethod("", "kernelTest", _plugin_kernelTest);
		PbWrapperRegistry::instance().addMethod("", "getCurl", _plugin_getCurl);
		PbWrapperRegistry::instance().addMethod("", "setinflow", _plugin_setinflow);
		PbWrapperRegistry::instance().addClass("Solver", "FluidSolver", "PbClass");
		PbWrapperRegistry::instance().addConstructor("FluidSolver", _FluidSolver_FluidSolver);
		PbWrapperRegistry::instance().addMethod("FluidSolver", "getGridSize", _FluidSolver_getGridSize);
		PbWrapperRegistry::instance().addMethod("FluidSolver", "printTimings", _FluidSolver_printTimings);
		PbWrapperRegistry::instance().addMethod("FluidSolver", "saveMeanTimings", _FluidSolver_saveMeanTimings);
		PbWrapperRegistry::instance().addMethod("FluidSolver", "step", _FluidSolver_step);
		PbWrapperRegistry::instance().addMethod("FluidSolver", "create", _FluidSolver_create);
		PbWrapperRegistry::instance().addGetSet("FluidSolver","timestep",_get_FluidSolver_mDt,_set_FluidSolver_mDt);
		PbWrapperRegistry::instance().addClass("GridBase", "GridBase", "PbClass");
		PbWrapperRegistry::instance().addConstructor("GridBase", _GridBase_GridBase);
		PbWrapperRegistry::instance().addClass("IntGrid", "Grid<int>", "GridBase");
		PbWrapperRegistry::instance().addConstructor("Grid<int>", _GridIntGrid_GridIntGrid);
		PbWrapperRegistry::instance().addMethod("Grid<int>", "save", _GridIntGrid_save);
		PbWrapperRegistry::instance().addMethod("Grid<int>", "load", _GridIntGrid_load);
		PbWrapperRegistry::instance().addMethod("Grid<int>", "add", _GridIntGrid_add);
		PbWrapperRegistry::instance().addMethod("Grid<int>", "scale", _GridIntGrid_scale);
		PbWrapperRegistry::instance().addMethod("Grid<int>", "copyFrom", _GridIntGrid_copyFrom);
		PbWrapperRegistry::instance().addClass("RealGrid", "Grid<Real>", "GridBase");
		PbWrapperRegistry::instance().addConstructor("Grid<Real>", _GridRealGrid_GridRealGrid);
		PbWrapperRegistry::instance().addMethod("Grid<Real>", "save", _GridRealGrid_save);
		PbWrapperRegistry::instance().addMethod("Grid<Real>", "load", _GridRealGrid_load);
		PbWrapperRegistry::instance().addMethod("Grid<Real>", "add", _GridRealGrid_add);
		PbWrapperRegistry::instance().addMethod("Grid<Real>", "scale", _GridRealGrid_scale);
		PbWrapperRegistry::instance().addMethod("Grid<Real>", "copyFrom", _GridRealGrid_copyFrom);
		PbWrapperRegistry::instance().addClass("VecGrid", "Grid<Vec3>", "GridBase");
		PbWrapperRegistry::instance().addConstructor("Grid<Vec3>", _GridVecGrid_GridVecGrid);
		PbWrapperRegistry::instance().addMethod("Grid<Vec3>", "save", _GridVecGrid_save);
		PbWrapperRegistry::instance().addMethod("Grid<Vec3>", "load", _GridVecGrid_load);
		PbWrapperRegistry::instance().addMethod("Grid<Vec3>", "add", _GridVecGrid_add);
		PbWrapperRegistry::instance().addMethod("Grid<Vec3>", "scale", _GridVecGrid_scale);
		PbWrapperRegistry::instance().addMethod("Grid<Vec3>", "copyFrom", _GridVecGrid_copyFrom);
		PbWrapperRegistry::instance().addClass("MACGrid", "MACGrid", "Grid<Vec3>");
		PbWrapperRegistry::instance().addConstructor("MACGrid", _MACGrid_MACGrid);
		PbWrapperRegistry::instance().addClass("FlagGrid", "FlagGrid", "Grid<int>");
		PbWrapperRegistry::instance().addConstructor("FlagGrid", _FlagGrid_FlagGrid);
		PbWrapperRegistry::instance().addMethod("FlagGrid", "initDomain", _FlagGrid_initDomain);
		PbWrapperRegistry::instance().addMethod("FlagGrid", "initBoundaries", _FlagGrid_initBoundaries);
		PbWrapperRegistry::instance().addMethod("FlagGrid", "updateFromLevelset", _FlagGrid_updateFromLevelset);
		PbWrapperRegistry::instance().addMethod("FlagGrid", "fillGrid", _FlagGrid_fillGrid);
		PbWrapperRegistry::instance().addClass("Mesh", "Mesh", "PbClass");
		PbWrapperRegistry::instance().addConstructor("Mesh", _Mesh_Mesh);
		PbWrapperRegistry::instance().addMethod("Mesh", "load", _Mesh_load);
		PbWrapperRegistry::instance().addMethod("Mesh", "fromShape", _Mesh_fromShape);
		PbWrapperRegistry::instance().addMethod("Mesh", "save", _Mesh_save);
		PbWrapperRegistry::instance().addMethod("Mesh", "advectInGrid", _Mesh_advectInGrid);
		PbWrapperRegistry::instance().addMethod("Mesh", "scale", _Mesh_scale);
		PbWrapperRegistry::instance().addMethod("Mesh", "offset", _Mesh_offset);
		PbWrapperRegistry::instance().addClass("ParticleBase", "ParticleBase", "PbClass");
		PbWrapperRegistry::instance().addConstructor("ParticleBase", _ParticleBase_ParticleBase);
		PbWrapperRegistry::instance().addClass("_ParticleSystem_BasicParticleData", "ParticleSystem<BasicParticleData>", "ParticleBase");
		PbWrapperRegistry::instance().addConstructor("ParticleSystem<BasicParticleData>", _ParticleSystem_ParticleSystem_BasicParticleData_ParticleSystem_ParticleSystem_BasicParticleData);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<BasicParticleData>", "size", _ParticleSystem_ParticleSystem_BasicParticleData_size);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<BasicParticleData>", "setPos", _ParticleSystem_ParticleSystem_BasicParticleData_setPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<BasicParticleData>", "getPos", _ParticleSystem_ParticleSystem_BasicParticleData_getPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<BasicParticleData>", "advectInGrid", _ParticleSystem_ParticleSystem_BasicParticleData_advectInGrid);
		PbWrapperRegistry::instance().addClass("TracerParticleSystem", "TracerParticleSystem", "ParticleSystem<BasicParticleData>");
		PbWrapperRegistry::instance().addConstructor("TracerParticleSystem", _TracerParticleSystem_TracerParticleSystem);
		PbWrapperRegistry::instance().addMethod("TracerParticleSystem", "addParticle", _TracerParticleSystem_addParticle);
		PbWrapperRegistry::instance().addMethod("TracerParticleSystem", "add_Particle", _TracerParticleSystem_add_Particle);
		PbWrapperRegistry::instance().addClass("LevelsetGrid", "LevelsetGrid", "Grid<Real>");
		PbWrapperRegistry::instance().addConstructor("LevelsetGrid", _LevelsetGrid_LevelsetGrid);
		PbWrapperRegistry::instance().addMethod("LevelsetGrid", "reinitMarching", _LevelsetGrid_reinitMarching);
		PbWrapperRegistry::instance().addMethod("LevelsetGrid", "createMesh", _LevelsetGrid_createMesh);
		PbWrapperRegistry::instance().addMethod("LevelsetGrid", "join", _LevelsetGrid_join);
		PbWrapperRegistry::instance().addMethod("LevelsetGrid", "initFromFlags", _LevelsetGrid_initFromFlags);
		PbWrapperRegistry::instance().addClass("Shape", "Shape", "PbClass");
		PbWrapperRegistry::instance().addConstructor("Shape", _Shape_Shape);
		PbWrapperRegistry::instance().addMethod("Shape", "applyToGrid", _Shape_applyToGrid);
		PbWrapperRegistry::instance().addMethod("Shape", "applyToGridSmooth", _Shape_applyToGridSmooth);
		PbWrapperRegistry::instance().addMethod("Shape", "computeLevelset", _Shape_computeLevelset);
		PbWrapperRegistry::instance().addMethod("Shape", "collideMesh", _Shape_collideMesh);
		PbWrapperRegistry::instance().addClass("NullShape", "NullShape", "Shape");
		PbWrapperRegistry::instance().addConstructor("NullShape", _NullShape_NullShape);
		PbWrapperRegistry::instance().addClass("Box", "Box", "Shape");
		PbWrapperRegistry::instance().addConstructor("Box", _Box_Box);
		PbWrapperRegistry::instance().addClass("Sphere", "Sphere", "Shape");
		PbWrapperRegistry::instance().addConstructor("Sphere", _Sphere_Sphere);
		PbWrapperRegistry::instance().addClass("Cylinder", "Cylinder", "Shape");
		PbWrapperRegistry::instance().addConstructor("Cylinder", _Cylinder_Cylinder);
		PbWrapperRegistry::instance().addMethod("Cylinder", "setCenter", _Cylinder_setCenter);
		PbWrapperRegistry::instance().addMethod("Cylinder", "setRadius", _Cylinder_setRadius);
		PbWrapperRegistry::instance().addMethod("Cylinder", "setZ", _Cylinder_setZ);
		PbWrapperRegistry::instance().addClass("NoiseField", "WaveletNoiseField", "PbClass");
		PbWrapperRegistry::instance().addConstructor("WaveletNoiseField", _WaveletNoiseField_WaveletNoiseField);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","posOffset",_get_WaveletNoiseField_mPosOffset,_set_WaveletNoiseField_mPosOffset);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","posScale",_get_WaveletNoiseField_mPosScale,_set_WaveletNoiseField_mPosScale);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","valOffset",_get_WaveletNoiseField_mValOffset,_set_WaveletNoiseField_mValOffset);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","valScale",_get_WaveletNoiseField_mValScale,_set_WaveletNoiseField_mValScale);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","clamp",_get_WaveletNoiseField_mClamp,_set_WaveletNoiseField_mClamp);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","clampNeg",_get_WaveletNoiseField_mClampNeg,_set_WaveletNoiseField_mClampNeg);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","clampPos",_get_WaveletNoiseField_mClampPos,_set_WaveletNoiseField_mClampPos);
		PbWrapperRegistry::instance().addGetSet("WaveletNoiseField","timeAnim",_get_WaveletNoiseField_mTimeAnim,_set_WaveletNoiseField_mTimeAnim);
		PbWrapperRegistry::instance().addClass("VortexSheetMesh", "VortexSheetMesh", "Mesh");
		PbWrapperRegistry::instance().addConstructor("VortexSheetMesh", _VortexSheetMesh_VortexSheetMesh);
		PbWrapperRegistry::instance().addMethod("VortexSheetMesh", "calcCirculation", _VortexSheetMesh_calcCirculation);
		PbWrapperRegistry::instance().addMethod("VortexSheetMesh", "calcVorticity", _VortexSheetMesh_calcVorticity);
		PbWrapperRegistry::instance().addMethod("VortexSheetMesh", "reinitTexCoords", _VortexSheetMesh_reinitTexCoords);
		PbWrapperRegistry::instance().addClass("_ConnectedParticleSystem_BasicParticleData_VortexRing", "ConnectedParticleSystem<BasicParticleData,VortexRing>", "ParticleSystem<BasicParticleData>");
		PbWrapperRegistry::instance().addConstructor("ConnectedParticleSystem<BasicParticleData,VortexRing>", _ConnectedParticleSystem_ConnectedParticleSystem_BasicParticleData_VortexRing_ConnectedParticleSystem_ConnectedParticleSystem_BasicParticleData_VortexRing);
		PbWrapperRegistry::instance().addClass("VortexFilamentSystem", "VortexFilamentSystem", "ConnectedParticleSystem<BasicParticleData,VortexRing>");
		PbWrapperRegistry::instance().addConstructor("VortexFilamentSystem", _VortexFilamentSystem_VortexFilamentSystem);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "advectSelf", _VortexFilamentSystem_advectSelf);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "advectParticles", _VortexFilamentSystem_advectParticles);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "advectMesh", _VortexFilamentSystem_advectMesh);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "doublyDiscreteUpdate", _VortexFilamentSystem_doublyDiscreteUpdate);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "remesh", _VortexFilamentSystem_remesh);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "split_ring", _VortexFilamentSystem_split_ring);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "reconnect_ring", _VortexFilamentSystem_reconnect_ring);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "merge_adj_edge", _VortexFilamentSystem_merge_adj_edge);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "divide_ring", _VortexFilamentSystem_divide_ring);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "reset_dirty", _VortexFilamentSystem_reset_dirty);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "Debug_fun", _VortexFilamentSystem_Debug_fun);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "revise_circulation", _VortexFilamentSystem_revise_circulation);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "Decimate_ring", _VortexFilamentSystem_Decimate_ring);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "Direct_motion", _VortexFilamentSystem_Direct_motion);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "addRing", _VortexFilamentSystem_addRing);
		PbWrapperRegistry::instance().addMethod("VortexFilamentSystem", "addLine", _VortexFilamentSystem_addLine);
		PbWrapperRegistry::instance().addClass("_ParticleSystem_FlipData", "ParticleSystem<FlipData>", "ParticleBase");
		PbWrapperRegistry::instance().addConstructor("ParticleSystem<FlipData>", _ParticleSystem_ParticleSystem_FlipData_ParticleSystem_ParticleSystem_FlipData);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<FlipData>", "size", _ParticleSystem_ParticleSystem_FlipData_size);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<FlipData>", "setPos", _ParticleSystem_ParticleSystem_FlipData_setPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<FlipData>", "getPos", _ParticleSystem_ParticleSystem_FlipData_getPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<FlipData>", "advectInGrid", _ParticleSystem_ParticleSystem_FlipData_advectInGrid);
		PbWrapperRegistry::instance().addClass("FlipSystem", "FlipSystem", "ParticleSystem<FlipData>");
		PbWrapperRegistry::instance().addConstructor("FlipSystem", _FlipSystem_FlipSystem);
		PbWrapperRegistry::instance().addMethod("FlipSystem", "velocitiesFromGrid", _FlipSystem_velocitiesFromGrid);
		PbWrapperRegistry::instance().addMethod("FlipSystem", "velocitiesToGrid", _FlipSystem_velocitiesToGrid);
		PbWrapperRegistry::instance().addMethod("FlipSystem", "adjustNumber", _FlipSystem_adjustNumber);
		PbWrapperRegistry::instance().addMethod("FlipSystem", "markFluidCells", _FlipSystem_markFluidCells);
		PbWrapperRegistry::instance().addClass("_ParticleSystem_VortexParticleData", "ParticleSystem<VortexParticleData>", "ParticleBase");
		PbWrapperRegistry::instance().addConstructor("ParticleSystem<VortexParticleData>", _ParticleSystem_ParticleSystem_VortexParticleData_ParticleSystem_ParticleSystem_VortexParticleData);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<VortexParticleData>", "size", _ParticleSystem_ParticleSystem_VortexParticleData_size);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<VortexParticleData>", "setPos", _ParticleSystem_ParticleSystem_VortexParticleData_setPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<VortexParticleData>", "getPos", _ParticleSystem_ParticleSystem_VortexParticleData_getPos);
		PbWrapperRegistry::instance().addMethod("ParticleSystem<VortexParticleData>", "advectInGrid", _ParticleSystem_ParticleSystem_VortexParticleData_advectInGrid);
		PbWrapperRegistry::instance().addClass("VortexParticleSystem", "VortexParticleSystem", "ParticleSystem<VortexParticleData>");
		PbWrapperRegistry::instance().addConstructor("VortexParticleSystem", _VortexParticleSystem_VortexParticleSystem);
		PbWrapperRegistry::instance().addMethod("VortexParticleSystem", "advectSelf", _VortexParticleSystem_advectSelf);
		PbWrapperRegistry::instance().addMethod("VortexParticleSystem", "applyToMesh", _VortexParticleSystem_applyToMesh);
		PbWrapperRegistry::instance().addClass("CustomControl", "CustomControl", "PbClass");
		PbWrapperRegistry::instance().addConstructor("CustomControl", _CustomControl_CustomControl);
		PbWrapperRegistry::instance().addClass("Slider", "CustomSlider", "CustomControl");
		PbWrapperRegistry::instance().addConstructor("CustomSlider", _CustomSlider_CustomSlider);
		PbWrapperRegistry::instance().addMethod("CustomSlider", "get", _CustomSlider_get);
		PbWrapperRegistry::instance().addMethod("CustomSlider", "set", _CustomSlider_set);
		PbWrapperRegistry::instance().addClass("Gui", "Gui", "PbClass");
		PbWrapperRegistry::instance().addConstructor("Gui", _Gui_Gui);
		PbWrapperRegistry::instance().addMethod("Gui", "setBackgroundMesh", _Gui_setBackgroundMesh);
		PbWrapperRegistry::instance().addMethod("Gui", "show", _Gui_show);
		PbWrapperRegistry::instance().addMethod("Gui", "update", _Gui_update);
		PbWrapperRegistry::instance().addMethod("Gui", "pause", _Gui_pause);
		PbWrapperRegistry::instance().addMethod("Gui", "addControl", _Gui_addControl);
		PbWrapperRegistry::instance().addMethod("Gui", "screenshot", _Gui_screenshot);
	}
};
static const Manta::_pbRegister _doRegister(12);
} //namespace
