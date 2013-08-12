/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Wavelet noise field
 *
 ******************************************************************************/
 
#include "vectorbase.h"
#include "pclass.h"

namespace Manta {

#define NOISE_TILE_SIZE 128

// wrapper for a parametrized field of wavelet noise
PYTHON(name=NoiseField) 
class WaveletNoiseField : public PbClass {
    public:     
        PYTHON WaveletNoiseField(FluidSolver* parent);
        ~WaveletNoiseField() {};

        //! evaluate noise
        inline Real evaluate(Vec3 pos);
        //! evaluate noise as a vector
        inline Vec3 evaluateVec(Vec3 pos);

        //! direct data access
        float* data() { return mNoiseTile; }
        
        // helper
        std::string toString();
        
        // texcoord position and scale
        PYTHON(name=posOffset) Vec3 mPosOffset;
        PYTHON(name=posScale) Vec3 mPosScale;
        // value offset & scale
        PYTHON(name=valOffset) Real mValOffset;
        PYTHON(name=valScale) Real mValScale;
        // clamp? (default 0-1)
        PYTHON(name=clamp) bool mClamp;
        PYTHON(name=clampNeg) Real mClampNeg;
        PYTHON(name=clampPos) Real mClampPos;
        
        PYTHON(name=timeAnim) Real mTimeAnim;        
    protected:
        static inline float WNoiseDx(Vec3 p, float *data);
        static inline Vec3 WNoiseVec(const Vec3& p, float *data);
        static inline float WNoise(const Vec3& p, float *data);
        static void Downsample(float *from, float *to, int n, int stride);
        static void Upsample(float *from, float *to, int n, int stride);
        static inline int modSlow(int x, int n) { int m = x % n; return (m<0) ? m+n : m; }
        // warning - noiseTileSize has to be 128^3!
        #define modFast128(x)  ((x) & 127)

        inline Real getTime() { return mParent->getTime() * mParent->getDx() * mTimeAnim; }
        void generateTile();
                
        // animation over time
        // grid size normalization (inverse size)
        Real mGsInvX, mGsInvY, mGsInvZ;
        
        float* mNoiseTile;
        static int seed;
};



// **************************************************************************
// Implementation

#define ADD_WEIGHTED(x,y,z)\
  weight = 1.0f;\
  xC = modFast128(midX + (x));\
  weight *= w[0][(x) + 1];\
  yC = modFast128(midY + (y));\
  weight *= w[1][(y) + 1];\
  zC = modFast128(midZ + (z));\
  weight *= w[2][(z) + 1];\
  result += weight * data[(zC * NOISE_TILE_SIZE + yC) * NOISE_TILE_SIZE + xC];

//////////////////////////////////////////////////////////////////////////////////////////
// derivatives of 3D noise - unrolled for performance
//////////////////////////////////////////////////////////////////////////////////////////
inline float WaveletNoiseField::WNoiseDx(Vec3 p, float *data) {
    float w[3][3], t, result = 0;

    // Evaluate quadratic B-spline basis functions
    int midX = (int)ceil(p[0] - 0.5f); 
    t        =   midX - (p[0] - 0.5f);
    w[0][0] = -t;
    w[0][2] = (1.f - t);
    w[0][1] = 2.0f * t - 1.0f;

    int midY = (int)ceil(p[1] - 0.5f); 
    t        =   midY - (p[1] - 0.5f);
    w[1][0] = t * t * 0.5f; 
    w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
    w[1][1] = 1.f - w[1][0] - w[1][2];

    int midZ = (int)ceil(p[2] - 0.5f); 
    t        =   midZ - (p[2] - 0.5f);
    w[2][0] = t * t * 0.5f; 
    w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
    w[2][1] = 1.f - w[2][0] - w[2][2];

    // Evaluate noise by weighting noise coefficients by basis function values
    int xC, yC, zC;
    float weight = 1;

    ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
    ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
    ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

    ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
    ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
    ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

    ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
    ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
    ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

    return result;
}

inline float WaveletNoiseField::WNoise(const Vec3& p, float *data) {
    float w[3][3], t, result = 0;

    // Evaluate quadratic B-spline basis functions
    int midX = (int)ceilf(p[0] - 0.5f); 
    t        =   midX - (p[0] - 0.5f);
    w[0][0] = t * t * 0.5f; 
    w[0][2] = (1.f - t) * (1.f - t) *0.5f; 
    w[0][1] = 1.f - w[0][0] - w[0][2];

    int midY = (int)ceilf(p[1] - 0.5f); 
    t        =   midY - (p[1] - 0.5f);
    w[1][0] = t * t * 0.5f; 
    w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
    w[1][1] = 1.f - w[1][0] - w[1][2];

    int midZ = (int)ceilf(p[2] - 0.5f); 
    t        =   midZ - (p[2] - 0.5f);
    w[2][0] = t * t * 0.5f; 
    w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
    w[2][1] = 1.f - w[2][0] - w[2][2];

    // Evaluate noise by weighting noise coefficients by basis function values
    int xC, yC, zC;
    float weight = 1;

    ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
    ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
    ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

    ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
    ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
    ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

    ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
    ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
    ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

    return result;
}



#define ADD_WEIGHTEDX(x,y,z)\
  weight = dw[0][(x) + 1] * w[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDY(x,y,z)\
  weight = w[0][(x) + 1] * dw[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDZ(x,y,z)\
  weight = w[0][(x) + 1] * w[1][(y) + 1] * dw[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

//////////////////////////////////////////////////////////////////////////////////////////
// compute all derivatives in at once
//////////////////////////////////////////////////////////////////////////////////////////
inline Vec3 WaveletNoiseField::WNoiseVec(const Vec3& p, float *data)
{
    Vec3 final(0.);
    float w[3][3];
    float dw[3][3];
    float result = 0;
    int xC, yC, zC;
    float weight;

    int midX = (int)ceil(p[0] - 0.5f); 
    int midY = (int)ceil(p[1] - 0.5f); 
    int midZ = (int)ceil(p[2] - 0.5f);

    float t0 =   midX - (p[0] - 0.5f);
    float t1 =   midY - (p[1] - 0.5f);
    float t2 =   midZ - (p[2] - 0.5f);

    // precache all the neighbors for fast access
    float neighbors[3][3][3];
    for (int z = -1; z <=1; z++)
        for (int y = -1; y <= 1; y++)
            for (int x = -1; x <= 1; x++)
            {
                xC = modFast128(midX + (x));
                yC = modFast128(midY + (y));
                zC = modFast128(midZ + (z));
                neighbors[x + 1][y + 1][z + 1] = data[zC * NOISE_TILE_SIZE * NOISE_TILE_SIZE + yC * NOISE_TILE_SIZE + xC];
            }

    ///////////////////////////////////////////////////////////////////////////////////////
    // evaluate splines
    ///////////////////////////////////////////////////////////////////////////////////////
    dw[0][0] = -t0;
    dw[0][2] = (1.f - t0);
    dw[0][1] = 2.0f * t0 - 1.0f;

    dw[1][0] = -t1;
    dw[1][2] = (1.0f - t1);
    dw[1][1] = 2.0f * t1 - 1.0f;

    dw[2][0] = -t2;
    dw[2][2] = (1.0f - t2);
    dw[2][1] = 2.0f * t2 - 1.0f;

    w[0][0] = t0 * t0 * 0.5f; 
    w[0][2] = (1.f - t0) * (1.f - t0) *0.5f; 
    w[0][1] = 1.f - w[0][0] - w[0][2];

    w[1][0] = t1 * t1 * 0.5f; 
    w[1][2] = (1.f - t1) * (1.f - t1) *0.5f; 
    w[1][1] = 1.f - w[1][0] - w[1][2];

    w[2][0] = t2 * t2 * 0.5f; 
    w[2][2] = (1.f - t2) * (1.f - t2) *0.5f;
    w[2][1] = 1.f - w[2][0] - w[2][2];

    ///////////////////////////////////////////////////////////////////////////////////////
    // x derivative
    ///////////////////////////////////////////////////////////////////////////////////////
    result = 0.0f;
    ADD_WEIGHTEDX(-1,-1, -1); ADD_WEIGHTEDX( 0,-1, -1); ADD_WEIGHTEDX( 1,-1, -1);
    ADD_WEIGHTEDX(-1, 0, -1); ADD_WEIGHTEDX( 0, 0, -1); ADD_WEIGHTEDX( 1, 0, -1);
    ADD_WEIGHTEDX(-1, 1, -1); ADD_WEIGHTEDX( 0, 1, -1); ADD_WEIGHTEDX( 1, 1, -1);

    ADD_WEIGHTEDX(-1,-1, 0);  ADD_WEIGHTEDX( 0,-1, 0);  ADD_WEIGHTEDX( 1,-1, 0);
    ADD_WEIGHTEDX(-1, 0, 0);  ADD_WEIGHTEDX( 0, 0, 0);  ADD_WEIGHTEDX( 1, 0, 0);
    ADD_WEIGHTEDX(-1, 1, 0);  ADD_WEIGHTEDX( 0, 1, 0);  ADD_WEIGHTEDX( 1, 1, 0);

    ADD_WEIGHTEDX(-1,-1, 1);  ADD_WEIGHTEDX( 0,-1, 1);  ADD_WEIGHTEDX( 1,-1, 1);
    ADD_WEIGHTEDX(-1, 0, 1);  ADD_WEIGHTEDX( 0, 0, 1);  ADD_WEIGHTEDX( 1, 0, 1);
    ADD_WEIGHTEDX(-1, 1, 1);  ADD_WEIGHTEDX( 0, 1, 1);  ADD_WEIGHTEDX( 1, 1, 1);
    final[0] = result;

    ///////////////////////////////////////////////////////////////////////////////////////
    // y derivative
    ///////////////////////////////////////////////////////////////////////////////////////
    result = 0.0f;
    ADD_WEIGHTEDY(-1,-1, -1); ADD_WEIGHTEDY( 0,-1, -1); ADD_WEIGHTEDY( 1,-1, -1);
    ADD_WEIGHTEDY(-1, 0, -1); ADD_WEIGHTEDY( 0, 0, -1); ADD_WEIGHTEDY( 1, 0, -1);
    ADD_WEIGHTEDY(-1, 1, -1); ADD_WEIGHTEDY( 0, 1, -1); ADD_WEIGHTEDY( 1, 1, -1);

    ADD_WEIGHTEDY(-1,-1, 0);  ADD_WEIGHTEDY( 0,-1, 0);  ADD_WEIGHTEDY( 1,-1, 0);
    ADD_WEIGHTEDY(-1, 0, 0);  ADD_WEIGHTEDY( 0, 0, 0);  ADD_WEIGHTEDY( 1, 0, 0);
    ADD_WEIGHTEDY(-1, 1, 0);  ADD_WEIGHTEDY( 0, 1, 0);  ADD_WEIGHTEDY( 1, 1, 0);

    ADD_WEIGHTEDY(-1,-1, 1);  ADD_WEIGHTEDY( 0,-1, 1);  ADD_WEIGHTEDY( 1,-1, 1);
    ADD_WEIGHTEDY(-1, 0, 1);  ADD_WEIGHTEDY( 0, 0, 1);  ADD_WEIGHTEDY( 1, 0, 1);
    ADD_WEIGHTEDY(-1, 1, 1);  ADD_WEIGHTEDY( 0, 1, 1);  ADD_WEIGHTEDY( 1, 1, 1);
    final[1] = result;

    ///////////////////////////////////////////////////////////////////////////////////////
    // z derivative
    ///////////////////////////////////////////////////////////////////////////////////////
    result = 0.0f;
    ADD_WEIGHTEDZ(-1,-1, -1); ADD_WEIGHTEDZ( 0,-1, -1); ADD_WEIGHTEDZ( 1,-1, -1);
    ADD_WEIGHTEDZ(-1, 0, -1); ADD_WEIGHTEDZ( 0, 0, -1); ADD_WEIGHTEDZ( 1, 0, -1);
    ADD_WEIGHTEDZ(-1, 1, -1); ADD_WEIGHTEDZ( 0, 1, -1); ADD_WEIGHTEDZ( 1, 1, -1);

    ADD_WEIGHTEDZ(-1,-1, 0);  ADD_WEIGHTEDZ( 0,-1, 0);  ADD_WEIGHTEDZ( 1,-1, 0);
    ADD_WEIGHTEDZ(-1, 0, 0);  ADD_WEIGHTEDZ( 0, 0, 0);  ADD_WEIGHTEDZ( 1, 0, 0);
    ADD_WEIGHTEDZ(-1, 1, 0);  ADD_WEIGHTEDZ( 0, 1, 0);  ADD_WEIGHTEDZ( 1, 1, 0);

    ADD_WEIGHTEDZ(-1,-1, 1);  ADD_WEIGHTEDZ( 0,-1, 1);  ADD_WEIGHTEDZ( 1,-1, 1);
    ADD_WEIGHTEDZ(-1, 0, 1);  ADD_WEIGHTEDZ( 0, 0, 1);  ADD_WEIGHTEDZ( 1, 0, 1);
    ADD_WEIGHTEDZ(-1, 1, 1);  ADD_WEIGHTEDZ( 0, 1, 1);  ADD_WEIGHTEDZ( 1, 1, 1);
    final[2] = result;

    //debMsg("FINAL","at "<<p<<" = "<<final); // DEBUG
    return final;
}
#undef ADD_WEIGHTEDX
#undef ADD_WEIGHTEDY
#undef ADD_WEIGHTEDZ

inline Real WaveletNoiseField::evaluate(Vec3 pos) { 
    pos[0] *= mGsInvX;
    pos[1] *= mGsInvY;
    pos[2] *= mGsInvZ;

    // time anim
    pos += Vec3(getTime());

    pos[0] *= mPosScale[0];
    pos[1] *= mPosScale[1];
    pos[2] *= mPosScale[2];
    pos += mPosOffset;

    Real v = WNoise(pos, mNoiseTile);

    v += mValOffset;
    v *= mValScale;
    if (mClamp) {
        if (v< mClampNeg) v = mClampNeg;
        if (v> mClampPos) v = mClampPos;
    }
    return v;
}

inline Vec3 WaveletNoiseField::evaluateVec(Vec3 pos) { 
    pos[0] *= mGsInvX;
    pos[1] *= mGsInvY;
    pos[2] *= mGsInvZ;

    // time anim
    pos += Vec3(getTime());

    pos[0] *= mPosScale[0];
    pos[1] *= mPosScale[1];
    pos[2] *= mPosScale[2];
    pos += mPosOffset;

    Vec3 v = WNoiseVec(pos, mNoiseTile);

    v += Vec3(mValOffset);
    v *= mValScale;
    if (mClamp) {
        for(int i=0; i<3; i++) {
            if (v[i]< mClampNeg) v[i] = mClampNeg;
            if (v[i]> mClampPos) v[i] = mClampPos;
        }
    }
    return v;
}



} // namespace WAVELETNOISE 

