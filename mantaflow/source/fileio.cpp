/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Loading and writing grids and meshes to disk
 *
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <zlib.h>
#include "fileio.h"
#include "grid.h"
#include "mesh.h"
#include "vortexsheet.h"
#include  <cstring>

using namespace std;

namespace Manta {
/*
void trim(string& s) {
    size_t pl = s.find_first_not_of("\t\r ");
    if (pl != string::npos) s=s.substr(pl);
    size_t pr = s.find_last_not_of("\t\r ");
    if (pr != string::npos) s=s.substr(0,pr);
}*/

void writeObjFile(const string& name, Mesh* mesh) {
    errMsg("obj exporter not yet implemented");
}

void writeBobjFile(const string& name, Mesh* mesh) {
    cout << "writing mesh file " << name << endl;
    const Real dx = mesh->getParent()->getDx();
    const Vec3i gs = mesh->getParent()->getGridSize();
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf)
        errMsg("writeBobj: unable to open file");
    
    // write vertices
    int numVerts = mesh->numNodes();
    gzwrite(gzf, &numVerts, sizeof(int));
    for (int i=0; i<numVerts; i++) {
        Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
        // normalize
        pos -= toVec3f(gs)*0.5;
        pos *= dx;
        gzwrite(gzf, &pos.value[0], sizeof(float)*3);
    }
    
    // normals
    mesh->computeVertexNormals();
    gzwrite(gzf, &numVerts, sizeof(int));
    for (int i=0; i<numVerts; i++) {
        Vector3D<float> pos = toVec3f(mesh->nodes(i).normal);
        gzwrite(gzf, &pos.value[0], sizeof(float)*3);
    }
    
    // write tris
    int numTris = mesh->numTris();
    gzwrite(gzf, &numTris, sizeof(int));
    for(int t=0; t<numTris; t++) {
        for(int j=0; j<3; j++) { 
            int trip = mesh->tris(t).c[j];
            gzwrite(gzf, &trip, sizeof(int)); 
        }
    }
    
    // per vertex smoke densities
    if (mesh->getType() == Mesh::TypeVortexSheet) {
        VortexSheetMesh* vmesh = (VortexSheetMesh*) mesh;
        int densId[4] = {0, 'v','d','e'};
        gzwrite(gzf, &densId[0], sizeof(int) * 4); 

        // compute densities
        vector<float> triDensity(numTris);
        for (int tri=0; tri < numTris; tri++) {
            Real area = vmesh->getFaceArea(tri);
            if (area>0)
                triDensity[tri] = vmesh->sheet(tri).smokeAmount;
        }
        
        // project triangle data to vertex
        vector<int> triPerVertex(numVerts);
        vector<float> density(numVerts);
        for (int tri=0; tri < numTris; tri++) {
            for (int c=0; c<3; c++) {
                int vertex = mesh->tris(tri).c[c];
                density[vertex] += triDensity[tri];
                triPerVertex[vertex]++;
            }
        }
        
        // averaged smoke densities
        for(int point=0; point<numVerts; point++) {
            float dens = 0;
            if (triPerVertex[point]>0)
                dens = density[point] / triPerVertex[point];
            gzwrite(gzf, &dens, sizeof(float));             
        }
    }
    
    // vertex flags
    if (mesh->getType() == Mesh::TypeVortexSheet) {
        int Id[4] = {0, 'v','x','f'};
        gzwrite(gzf, &Id[0], sizeof(int) * 4); 

        // averaged smoke densities
        for(int point=0; point<numVerts; point++) {
            float alpha = (mesh->nodes(point).flags & Mesh::NfMarked) ? 1: 0;
            gzwrite(gzf, &alpha, sizeof(float));             
        }
    }

    gzclose( gzf );    
}

void readObjFile(const std::string& name, Mesh* mesh, bool append) {
    ifstream ifs (name.c_str());
    
    if (!ifs.good())
        errMsg("can't open file '" + name + "'");
    
    if (!append)
        mesh->clear();
    int nodebase = mesh->numNodes();
    while(ifs.good() && !ifs.eof()) {
        string id;
        ifs >> id;
        
        if (id[0] == '#') {
            // comment
            getline(ifs, id);
            continue;
        }
        if (id == "vt") {
            // tex coord, ignore            
        } else if (id == "vn") {
            // normals, ignore            
        } else if (id == "v") {
            // vertex
            Node n;
            ifs >> n.pos.x >> n.pos.y >> n.pos.z;
            mesh->addNode(n);
        } else if (id == "g") {
            // group
            string group;
            ifs >> group;
        } else if (id == "f") {
            // face
            string face;
            Triangle t;
            for (int i=0; i<3; i++) {
                ifs >> face;
                if (face.find('/') != string::npos)
                    face = face.substr(0, face.find('/')); // ignore other indices
                int idx = atoi(face.c_str()) - 1;
                if (idx < 0)
                    errMsg("invalid face encountered");
                idx += nodebase;
                t.c[i] = idx;
            }
            mesh->addTri(t);
        } else {
            // whatever, ignore
        }
        // kill rest of line
        getline(ifs, id);   
    }
    ifs.close();    
}

template<class T>
void writeGridRaw(const string& name, Grid<T>* grid) {
    cout << "writing grid " << grid->getName() << " to raw file " << name << endl;
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) errMsg("can't open file");
    gzwrite(gzf, &((*grid)[0]), sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ());
    gzclose(gzf);
}

template<class T>
void readGridRaw(const string& name, Grid<T>* grid) {
    cout << "reading grid " << grid->getName() << " from raw file " << name << endl;
    
    gzFile gzf = gzopen(name.c_str(), "rb");
    if (!gzf) errMsg("can't open file");
    
    int bytes = sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ();
    int readBytes = gzread(gzf, &((*grid)[0]), bytes);
    assertMsg(bytes==readBytes, "can't read raw file, stream length does not match");
    gzclose(gzf);
}


typedef struct {
    int dimX, dimY, dimZ;
    int frames, elements, elementType, bytesPerElement, bytesPerFrame;
} UniLegacyHeader;

typedef struct {
    int dimX, dimY, dimZ;
    int gridType, elementType, bytesPerElement;
} UniHeader;

template <class T>
void writeGridUni(const string& name, Grid<T>* grid) {
    cout << "writing grid " << grid->getName() << " to uni file " << name << endl;
    
    char ID[5] = "MNT1";
    UniHeader head;
	head.dimX = grid->getSizeX();
    head.dimY = grid->getSizeY();
    head.dimZ = grid->getSizeZ();
    head.gridType = grid->getType();
    head.bytesPerElement = sizeof(T);
    
    if (grid->getType() & GridBase::TypeInt)
        head.elementType = 0;
    else if (grid->getType() & GridBase::TypeReal)
        head.elementType = 1;
    else if (grid->getType() & GridBase::TypeVec3)
        head.elementType = 2;
    else 
        errMsg("unknown element type");
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) errMsg("can't open file");
    
    gzwrite(gzf, ID, 4);
    gzwrite(gzf, &head, sizeof(UniHeader));
    gzwrite(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
    gzclose(gzf);
};

template <class T>
void readGridUni(const string& name, Grid<T>* grid) {
    cout << "reading grid " << grid->getName() << " from uni file " << name << endl;
    
    gzFile gzf = gzopen(name.c_str(), "rb");
    if (!gzf) errMsg("can't open file");
    
    char ID[5]={0,0,0,0,0};
    gzread(gzf, ID, 4);
    
    if (!strcmp(ID, "DDF2")) {
        // legacy file format
        UniLegacyHeader head;
        assertMsg (gzread(gzf, &head, sizeof(UniLegacyHeader)) == sizeof(UniLegacyHeader), "can't read file, no header present");
        assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match");
        assertMsg (head.bytesPerElement * head.elements == sizeof(T), "grid type doesn't match");
        // skip flags
        int numEl = head.dimX*head.dimY*head.dimZ;
        gzseek(gzf, numEl, SEEK_CUR);
        // actual grid read
        gzread(gzf, &((*grid)[0]), sizeof(T)*numEl);
    } 
    else if (!strcmp(ID, "MNT1")) {
        // current file format
        UniHeader head;
        assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present");
        assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match");
        assertMsg (head.gridType == grid->getType() && head.bytesPerElement == sizeof(T), "grid type doesn't match");
        gzread(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
    }
    gzclose(gzf);
};

// explicit instantiation
template void writeGridRaw<int>(const string& name, Grid<int>* grid);
template void writeGridRaw<Real>(const string& name, Grid<Real>* grid);
template void writeGridRaw<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridUni<int>(const string& name, Grid<int>* grid);
template void writeGridUni<Real>(const string& name, Grid<Real>* grid);
template void writeGridUni<Vec3>(const string& name, Grid<Vec3>* grid);
template void readGridRaw<int>(const string& name, Grid<int>* grid);
template void readGridRaw<Real>(const string& name, Grid<Real>* grid);
template void readGridRaw<Vec3>(const string& name, Grid<Vec3>* grid);
template void readGridUni<int>(const string& name, Grid<int>* grid);
template void readGridUni<Real>(const string& name, Grid<Real>* grid);
template void readGridUni<Vec3>(const string& name, Grid<Vec3>* grid);


} //namespace
