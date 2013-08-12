




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/mzhang/mantaflow/source/gui/customctrl.h"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * GUI extension from python
 *
 ******************************************************************************/

#ifndef _CUSTOMCTRL_H__
#define _CUSTOMCTRL_H__

#include <QSlider>
#include <QLabel>
#include <QLayout>
#include "pclass.h"

namespace Manta {

// fwd decl.
class Mesh;
class GuiThread;
class MainThread;
    
//! Interface for python declared controls
class CustomControl : public PbClass {
public:
     CustomControl() ;
    
    virtual void init(QLayout* layout) {};

protected:
protected:PbArgs _args;};;

//! Slider with attached text display
class TextSlider : public QSlider {
    Q_OBJECT
public:
    TextSlider(const std::string& name, float val, float min, float max);
    void attach(QLayout* layout);
    void set(float v);
    float get();
    
public slots:
    void update(int v);
        
protected:
    float mMin, mMax, mScale;
    QLabel* mLabel;    
    QString mSName;    
};
    
//! Links a slider control

class CustomSlider : public CustomControl {
public:
     CustomSlider(std::string text,float val,float min,float max) ;
    virtual void init(QLayout* layout);
    
    float get() ;PyObject* _get (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "CustomSlider::get"); { ArgLocker _lock; this->_args.copy(__args); _retval = d_toPy(get() );this->_args.check(); } pbFinalizePlugin(this->mParent,"CustomSlider::get"); return _retval; } catch(std::exception& e) { pbSetError("CustomSlider::get",e.what()); return 0; } } 
    void set(float v) ;PyObject* _set (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "CustomSlider::set"); { ArgLocker _lock; float v = __args.get< float > (0,"v", &_lock); this->_args.copy(__args); _retval = getPyNone();set(v);this->_args.check(); } pbFinalizePlugin(this->mParent,"CustomSlider::set"); return _retval; } catch(std::exception& e) { pbSetError("CustomSlider::set",e.what()); return 0; } } 
    
protected:
    float mMin, mMax, mVal;
    std::string mSName;
    TextSlider* mSlider;
protected:PbArgs _args;};;
    

//! GUI adapter class to call from Python
class Gui : public PbClass {
public:
     Gui() ;
    
    void setBackgroundMesh(Mesh* m) ;PyObject* _setBackgroundMesh (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::setBackgroundMesh"); { ArgLocker _lock; Mesh* m = __args.get< Mesh* > (0,"m", &_lock); this->_args.copy(__args); _retval = getPyNone();setBackgroundMesh(m);this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::setBackgroundMesh"); return _retval; } catch(std::exception& e) { pbSetError("Gui::setBackgroundMesh",e.what()); return 0; } } 
    void show() ;PyObject* _show (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::show"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();show();this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::show"); return _retval; } catch(std::exception& e) { pbSetError("Gui::show",e.what()); return 0; } } 
    void update() ;PyObject* _update (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::update"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();update();this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::update"); return _retval; } catch(std::exception& e) { pbSetError("Gui::update",e.what()); return 0; } } 
    void pause() ;PyObject* _pause (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::pause"); { ArgLocker _lock; this->_args.copy(__args); _retval = getPyNone();pause();this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::pause"); return _retval; } catch(std::exception& e) { pbSetError("Gui::pause",e.what()); return 0; } } 
    PbClass* addControl(PbType t) ;PyObject* _addControl (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::addControl"); { ArgLocker _lock; PbType t = __args.get< PbType > (0,"t", &_lock); this->_args.copy(__args); _retval = d_toPy(addControl(t) );this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::addControl"); return _retval; } catch(std::exception& e) { pbSetError("Gui::addControl",e.what()); return 0; } } 
    void screenshot(std::string filename) ;PyObject* _screenshot (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs __args(_linargs, _kwds); PyObject *_retval = NULL; pbPreparePlugin(this->mParent, "Gui::screenshot"); { ArgLocker _lock; std::string filename = __args.get< std::string > (0,"filename", &_lock); this->_args.copy(__args); _retval = getPyNone();screenshot(filename);this->_args.check(); } pbFinalizePlugin(this->mParent,"Gui::screenshot"); return _retval; } catch(std::exception& e) { pbSetError("Gui::screenshot",e.what()); return 0; } } 
    
protected:
    GuiThread* mGuiPtr;
    MainThread* mMainPtr;
protected:PbArgs _args;};;
    
} // namespace

#endif



