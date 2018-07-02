/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef MONTE_CARLO_RENDERER_CL
#define MONTE_CARLO_RENDERER_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef COMMON_CL
#define COMMON_CL

#define PI 3.14159265358979323846f
#define KERNEL __kernel
#define GLOBAL __global

#ifndef APPLE
#define INLINE __attribute__((always_inline))
#endif

#define HIT_MARKER 1
#define MISS_MARKER -1
#define INVALID_IDX -1

#define CRAZY_LOW_THROUGHPUT 0.0f
#define CRAZY_HIGH_RADIANCE 10.f
#define CRAZY_HIGH_DISTANCE 1000000.f
#define CRAZY_LOW_DISTANCE 0.001f
#define CRAZY_HIGH_DISTANCE_IN_VOLUME 1000000.f
#define REASONABLE_RADIANCE(x) (clamp((x), 0.f, CRAZY_HIGH_RADIANCE))
#define NON_BLACK(x) (length(x) > 0.f)

#define MULTISCATTER

#define RANDOM 1
#define SOBOL 2
#define CMJ 3

#define SAMPLER CMJ

#define CMJ_DIM 16

#define BDPT_MAX_SUBPATH_LEN 3

#ifdef BAIKAL_ATOMIC_RESOLVE
#define ADD_FLOAT3(x,y) atomic_add_float3((x),(y))
#define ADD_FLOAT4(x,y) atomic_add_float4((x),(y))
#else
#define ADD_FLOAT3(x,y) add_float3((x),(y))
#define ADD_FLOAT4(x,y) add_float4((x),(y))
#endif

#define VISIBILITY_MASK_PRIMARY (0x1)
#define VISIBILITY_MASK_SHADOW (0x1 << 15)
#define VISIBILITY_MASK_ALL (0xffffffffu)
#define VISIBILITY_MASK_NONE (0x0u)
#define VISIBILITY_MASK_BOUNCE(i) (VISIBILITY_MASK_PRIMARY << (i))
#define VISIBILITY_MASK_BOUNCE_SHADOW(i) (VISIBILITY_MASK_SHADOW << (i))

#endif // COMMON_CL
/**********************************************************************
 Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 ********************************************************************/
#ifndef RAY_CL
#define RAY_CL


// Ray descriptor
typedef struct
{
    // xyz - origin, w - max range
    float4 o;
    // xyz - direction, w - time
    float4 d;
    // x - ray mask, y - activity flag
    int2 extra;
    // Padding
    float2 padding;
} ray;

// Set ray activity flag
INLINE void Ray_SetInactive(GLOBAL ray* r)
{
    r->extra.y = 0;
}

INLINE bool Ray_IsActive(GLOBAL ray* r)
{
    return r->extra.y != 0;
}

// Set extra data for ray
INLINE void Ray_SetExtra(GLOBAL ray* r, float2 extra)
{
    r->padding = extra;
}

// Set mask
INLINE void Ray_SetMask(GLOBAL ray* r, int mask)
{
    r->extra.x = mask;
}

INLINE int Ray_GetMask(GLOBAL ray* r)
{
    return r->extra.x;
}

// Get extra data for ray
INLINE float2 Ray_GetExtra(GLOBAL ray const* r)
{
    return r->padding;
}

// Initialize ray structure
INLINE void Ray_Init(GLOBAL ray* r, float3 o, float3 d, float maxt, float time, int mask)
{
    r->o.xyz = o;
    r->d.xyz = d;
    r->o.w = maxt;
    r->d.w = time;
    r->extra.x = mask;
    r->extra.y = 0xFFFFFFFF;
}

#endif
/**********************************************************************
 Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 ********************************************************************/
#ifndef ISECT_CL
#define ISECT_CL

/// Intersection data returned by RadeonRays
typedef struct _Intersection
{
    // id of a shape
    int shapeid;
    // Primitive index
    int primid;
    // Padding elements
    int padding0;
    int padding1;
        
    // uv - hit barycentrics, w - ray distance
    float4 uvwt;
} Intersection;

float Intersection_GetDistance(__global Intersection const* isect)
{
    return isect->uvwt.w;
}

float2 Intersection_GetBarycentrics(__global Intersection const* isect)
{
    return isect->uvwt.xy;
}

#endif
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef UTILS_CL
#define UTILS_CL

#define PI 3.14159265358979323846f
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef PAYLOAD_CL
#define PAYLOAD_CL

#define TEXTURED_INPUT(x) union { struct { float4 value; } float_value; struct { int value[4]; } int_value; } x
#define TEXTURED_INPUT_HAS_TEXTURE(x) ((x).int_value.value[3] != -1)
#define TEXTURED_INPUT_GET_COLOR(x) ((x).float_value.value.xyz)

// Matrix
typedef struct
{
    float4 m0;
    float4 m1;
    float4 m2;
    float4 m3;
} matrix4x4;

// Camera
typedef struct
{
    // Coordinate frame
    float3 forward;
    float3 right;
    float3 up;
    // Position
    float3 p;

    // Image plane width & height in current units
    float2 dim;
    // Near and far Z
    float2 zcap;
    // Focal lenght
    float focal_length;
    // Camera aspect_ratio ratio
    float aspect_ratio;
    float focus_distance;
    float aperture;
} Camera;

enum UberMaterialLayers
{
    kEmissionLayer = 0x1,
    kTransparencyLayer = 0x2,
    kCoatingLayer = 0x4,
    kReflectionLayer = 0x8,
    kDiffuseLayer = 0x10,
    kRefractionLayer = 0x20,
    kSSSLayer = 0x40,
    kShadingNormalLayer = 0x80
};

typedef struct
{
    int offset;
    int layers;
    int flags;
    int padding;
} Material;

// Shape description
typedef struct
{
    // Shape starting index
    int startidx;
    // Start vertex
    int startvtx;
    // Number of primitives in the shape
    int volume_idx;
    // unique shape id
    int id;
    // Linear motion vector
    float3 linearvelocity;
    // Angular velocity
    float4 angularvelocity;
    // Transform in row major format
    matrix4x4 transform;
    Material material;
} Shape;

typedef struct
{
    int group_id;
    int padding[3];
} ShapeAdditionalData;

typedef enum
{
    kFloat3 = 0,
    kFloat = 1,
    kInt = 2
} InputMapDataType;

// Input data for input maps
typedef struct _InputMapData
{
    union
    {
        struct
        {
            float3 value;
        } float_value;
        struct
        {
            int idx;
            int placeholder[2];
            int type; //We can use it since float3 is actually float4
        } int_values;
    };
} InputMapData;

enum Bxdf
{
    kZero,
    kUberV2
};

enum LightType
{
    kPoint = 0x1,
    kDirectional,
    kSpot,
    kArea,
    kIbl
};

typedef struct
{
    union
    {
        // Area light
        struct
        {
            int id;
            int shapeidx;
            int primidx;
            int padding0;
        };

        // IBL
        struct
        {
            int tex;
            int tex_reflection;
            int tex_refraction;
            int tex_transparency;
        };

        // Spot
        struct
        {
            float ia;
            float oa;
            float f;
            int padding1;
        };
    };

    float3 p;
    float3 d;
    float3 intensity;
    int type;
    float multiplier;
    int tex_background;
    bool ibl_mirror_x;
} Light;

typedef enum
    {
        kEmpty,
        kHomogeneous,
        kHeterogeneous
    } VolumeType;


typedef struct _Volume
{
    VolumeType type;
    float g;

    // Id of volume data if present
    int data;
    int extra;

    // Absorbtion
    TEXTURED_INPUT(sigma_a);
    // Scattering
    TEXTURED_INPUT(sigma_s);
    // Emission
    TEXTURED_INPUT(sigma_e);
} Volume;

/// Supported formats
enum TextureFormat
{
    UNKNOWN,
    RGBA8,
    RGBA16,
    RGBA32
};

/// Texture description
typedef
struct _Texture
{
    // Width, height and depth
    int w;
    int h;
    int d;
    // Offset in texture data array
    int dataoffset;
    // Format
    int fmt;
    int extra;
} Texture;

// Hit data
typedef struct _DifferentialGeometry
{
    // World space position
    float3 p;
    // Shading normal
    float3 n;
    // Geo normal
    float3 ng;
    // UVs
    float2 uv;
    // Derivatives
    float3 dpdu;
    float3 dpdv;

    matrix4x4 world_to_tangent;
    matrix4x4 tangent_to_world;

    // Material
    Material mat;
    float  area;
    int transfer_mode;
    int padding[2];
} DifferentialGeometry;










#endif // PAYLOAD_CL


#ifndef APPLE
/// These functions are defined on OSX already
float4 make_float4(float x, float y, float z, float w)
{
    float4 res;
    res.x = x;
    res.y = y;
    res.z = z;
    res.w = w;
    return res;
}

float3 make_float3(float x, float y, float z)
{
    float3 res;
    res.x = x;
    res.y = y;
    res.z = z;
    return res;
}

float2 make_float2(float x, float y)
{
    float2 res;
    res.x = x;
    res.y = y;
    return res;
}

int2 make_int2(int x, int y)
{
    int2 res;
    res.x = x;
    res.y = y;
    return res;
}
#endif

matrix4x4 matrix_from_cols(float4 c0, float4 c1, float4 c2, float4 c3)
{
    matrix4x4 m;
    m.m0 = make_float4(c0.x, c1.x, c2.x, c3.x);
    m.m1 = make_float4(c0.y, c1.y, c2.y, c3.y);
    m.m2 = make_float4(c0.z, c1.z, c2.z, c3.z);
    m.m3 = make_float4(c0.w, c1.w, c2.w, c3.w);
    return m;
}

matrix4x4 matrix_from_rows(float4 c0, float4 c1, float4 c2, float4 c3)
{
    matrix4x4 m;
    m.m0 = c0;
    m.m1 = c1;
    m.m2 = c2;
    m.m3 = c3;
    return m;
}

matrix4x4 matrix_from_rows3(float3 c0, float3 c1, float3 c2)
{
    matrix4x4 m;
    m.m0.xyz = c0; m.m0.w = 0;
    m.m1.xyz = c1; m.m1.w = 0;
    m.m2.xyz = c2; m.m2.w = 0;
    m.m3 = make_float4(0.f, 0.f, 0.f, 1.f);
    return m;
}

matrix4x4 matrix_from_cols3(float3 c0, float3 c1, float3 c2)
{
    matrix4x4 m;
    m.m0 = make_float4(c0.x, c1.x, c2.x, 0.f);
    m.m1 = make_float4(c0.y, c1.y, c2.y, 0.f);
    m.m2 = make_float4(c0.z, c1.z, c2.z, 0.f);
    m.m3 = make_float4(0.f, 0.f, 0.f, 1.f);
    return m;
}

matrix4x4 matrix_transpose(matrix4x4 m)
{
    return matrix_from_cols(m.m0, m.m1, m.m2, m.m3);
}

float4 matrix_mul_vector4(matrix4x4 m, float4 v)
{
    float4 res;
    res.x = dot(m.m0, v);
    res.y = dot(m.m1, v);
    res.z = dot(m.m2, v);
    res.w = dot(m.m3, v);
    return res;
}

float3 matrix_mul_vector3(matrix4x4 m, float3 v)
{
    float3 res;
    res.x = dot(m.m0.xyz, v);
    res.y = dot(m.m1.xyz, v);
    res.z = dot(m.m2.xyz, v);
    return res;
}

float3 matrix_mul_point3(matrix4x4 m, float3 v)
{
    float3 res;
    res.x = dot(m.m0.xyz, v) + m.m0.w;
    res.y = dot(m.m1.xyz, v) + m.m1.w;
    res.z = dot(m.m2.xyz, v) + m.m2.w;
    return res;
}

/// Linearly interpolate between two values
float4 lerp(float4 a, float4 b, float w)
{
    return a + w*(b-a);
}

/// Linearly interpolate between two values
float3 lerp3(float3 a, float3 b, float w)
{
	return a + w*(b - a);
}

/// Translate cartesian coordinates to spherical system
void CartesianToSpherical ( float3 cart, float* r, float* phi, float* theta )
{
    float temp = atan2(cart.x, cart.z);
    *r = sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
    // Account for discontinuity
    *phi = (float)((temp >= 0)?temp:(temp + 2*PI));
    *theta = acos(cart.y/ *r);
}

/// Get vector orthogonal to a given one
float3 GetOrthoVector(float3 n)
{
    float3 p;

    if (fabs(n.z) > 0.f) {
        float k = sqrt(n.y*n.y + n.z*n.z);
        p.x = 0; p.y = -n.z/k; p.z = n.y/k;
    }
    else {
        float k = sqrt(n.x*n.x + n.y*n.y);
        p.x = n.y/k; p.y = -n.x/k; p.z = 0;
    }

    return normalize(p);
}

float luminance(float3 v)
{
    // Luminance
    return 0.2126f * v.x + 0.7152f * v.y + 0.0722f * v.z;
}

uint upper_power_of_two(uint v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

INLINE
void atomic_add_float(volatile __global float* addr, float value)
{
    union {
        unsigned int u32;
        float        f32;
    } next, expected, current;
    current.f32 = *addr;
    do {
        expected.f32 = current.f32;
        next.f32 = expected.f32 + value;
        current.u32 = atomic_cmpxchg((volatile __global unsigned int *)addr,
            expected.u32, next.u32);
    } while (current.u32 != expected.u32);
}

void atomic_add_float3(volatile __global float3* ptr, float3 value)
{
    volatile __global float* p = (volatile __global float*)ptr;
    atomic_add_float(p, value.x);
    atomic_add_float(p + 1, value.y);
    atomic_add_float(p + 2, value.z);
}

void atomic_add_float4(volatile __global float4* ptr, float4 value)
{
    volatile __global float* p = (volatile __global float*)ptr;
    atomic_add_float(p, value.x);
    atomic_add_float(p + 1, value.y);
    atomic_add_float(p + 2, value.z);
    atomic_add_float(p + 3, value.w);
}

void add_float3(__global float3* ptr, float3 value)
{
    *ptr += value;
}

void add_float4(__global float4* ptr, float4 value)
{
    *ptr += value;
}


#endif // UTILS_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef TEXTURE_CL
#define TEXTURE_CL




/// To simplify a bit
#define TEXTURE_ARG_LIST __global Texture const* textures, __global char const* texturedata
#define TEXTURE_ARG_LIST_IDX(x) int x, __global Texture const* textures, __global char const* texturedata
#define TEXTURE_ARGS textures, texturedata
#define TEXTURE_ARGS_IDX(x) x, textures, texturedata

/// Sample 2D texture
inline
float4 Texture_Sample2D(float2 uv, TEXTURE_ARG_LIST_IDX(texidx))
{
    // Get width and height
    int width = textures[texidx].w;
    int height = textures[texidx].h;

    // Find the origin of the data in the pool
    __global char const* mydata = texturedata + textures[texidx].dataoffset;

    // Handle UV wrap
    // TODO: need UV mode support
    uv -= floor(uv);

    // Reverse Y:
    // it is needed as textures are loaded with Y axis going top to down
    // and our axis goes from down to top
    uv.y = 1.f - uv.y;

    // Calculate integer coordinates
    int x0 = clamp((int)floor(uv.x * width), 0, width - 1);
    int y0 = clamp((int)floor(uv.y * height), 0, height - 1);

    // Calculate samples for linear filtering
    int x1 = clamp(x0 + 1, 0,  width - 1);
    int y1 = clamp(y0 + 1, 0, height - 1);

    // Calculate weights for linear filtering
    float wx = uv.x * width - floor(uv.x * width);
    float wy = uv.y * height - floor(uv.y * height);

    switch (textures[texidx].fmt)
    {
        case RGBA32:
        {
            __global float4 const* mydataf = (__global float4 const*)mydata;

            // Get 4 values for linear filtering
            float4 val00 = *(mydataf + width * y0 + x0);
            float4 val01 = *(mydataf + width * y0 + x1);
            float4 val10 = *(mydataf + width * y1 + x0);
            float4 val11 = *(mydataf + width * y1 + x1);

            // Filter and return the result
            return lerp(lerp(val00, val01, wx), lerp(val10, val11, wx), wy);
        }

        case RGBA16:
        {
            __global half const* mydatah = (__global half const*)mydata;

            // Get 4 values
            float4 val00 = vload_half4(width * y0 + x0, mydatah);
            float4 val01 = vload_half4(width * y0 + x1, mydatah);
            float4 val10 = vload_half4(width * y1 + x0, mydatah);
            float4 val11 = vload_half4(width * y1 + x1, mydatah);

            // Filter and return the result
            return lerp(lerp(val00, val01, wx), lerp(val10, val11, wx), wy);
        }

        case RGBA8:
        {
            __global uchar4 const* mydatac = (__global uchar4 const*)mydata;

            // Get 4 values and convert to float
            uchar4 valu00 = *(mydatac + width * y0 + x0);
            uchar4 valu01 = *(mydatac + width * y0 + x1);
            uchar4 valu10 = *(mydatac + width * y1 + x0);
            uchar4 valu11 = *(mydatac + width * y1 + x1);

            float4 val00 = make_float4((float)valu00.x / 255.f, (float)valu00.y / 255.f, (float)valu00.z / 255.f, (float)valu00.w / 255.f);
            float4 val01 = make_float4((float)valu01.x / 255.f, (float)valu01.y / 255.f, (float)valu01.z / 255.f, (float)valu01.w / 255.f);
            float4 val10 = make_float4((float)valu10.x / 255.f, (float)valu10.y / 255.f, (float)valu10.z / 255.f, (float)valu10.w / 255.f);
            float4 val11 = make_float4((float)valu11.x / 255.f, (float)valu11.y / 255.f, (float)valu11.z / 255.f, (float)valu11.w / 255.f);

            // Filter and return the result
            return lerp(lerp(val00, val01, wx), lerp(val10, val11, wx), wy);
        }

        default:
        {
            return make_float4(0.f, 0.f, 0.f, 0.f);
        }
    }
}

/// Sample lattitue-longitude environment map using 3d vector
inline
float3 Texture_SampleEnvMap(float3 d, TEXTURE_ARG_LIST_IDX(texidx), bool mirror_x)
{
    // Transform to spherical coords
    float r, phi, theta;
    CartesianToSpherical(d, &r, &phi, &theta);

    // Map to [0,1]x[0,1] range and reverse Y axis
    float2 uv;
    uv.x = (mirror_x) ? (1.f - phi / (2 * PI)) : phi / (2 * PI);
    uv.y = 1.f - theta / PI;

    // Sample the texture
    return Texture_Sample2D(uv, TEXTURE_ARGS_IDX(texidx)).xyz;
}

/// Get data from parameter value or texture
inline
float3 Texture_GetValue3f(
                // Value
                float3 v,
                // Texture coordinate
                float2 uv,
                // Texture args
                TEXTURE_ARG_LIST_IDX(texidx)
                )
{
    // If texture present sample from texture
    if (texidx != -1)
    {
        // Sample texture
        return native_powr(Texture_Sample2D(uv, TEXTURE_ARGS_IDX(texidx)).xyz, 2.2f);
    }

    // Return fixed color otherwise
    return v;
}

/// Get data from parameter value or texture
inline
float4 Texture_GetValue4f(
                // Value
                float4 v,
                // Texture coordinate
                float2 uv,
                // Texture args
                TEXTURE_ARG_LIST_IDX(texidx)
                )
{
    // If texture present sample from texture
    if (texidx != -1)
    {
        // Sample texture
        return native_powr(Texture_Sample2D(uv, TEXTURE_ARGS_IDX(texidx)), 2.2f);
    }

    // Return fixed color otherwise
    return v;
}

/// Get data from parameter value or texture
inline
float Texture_GetValue1f(
                        // Value
                        float v,
                        // Texture coordinate
                        float2 uv,
                        // Texture args
                        TEXTURE_ARG_LIST_IDX(texidx)
                        )
{
    // If texture present sample from texture
    if (texidx != -1)
    {
        // Sample texture
        return Texture_Sample2D(uv, TEXTURE_ARGS_IDX(texidx)).x;
    }

    // Return fixed color otherwise
    return v;
}

inline float3 TextureData_SampleNormalFromBump_uchar4(__global uchar4 const* mydatac, int width, int height, int t0, int s0)
{
	int t0minus = clamp(t0 - 1, 0, height - 1);
	int t0plus = clamp(t0 + 1, 0, height - 1);
	int s0minus = clamp(s0 - 1, 0, width - 1);
	int s0plus = clamp(s0 + 1, 0, width - 1);

	const uchar utex00 = (*(mydatac + width * t0minus + s0minus)).x;
	const uchar utex10 = (*(mydatac + width * t0minus + (s0))).x;
	const uchar utex20 = (*(mydatac + width * t0minus + s0plus)).x;

	const uchar utex01 = (*(mydatac + width * (t0)+s0minus)).x;
	const uchar utex21 = (*(mydatac + width * (t0)+(s0 + 1))).x;

	const uchar utex02 = (*(mydatac + width * t0plus + s0minus)).x;
	const uchar utex12 = (*(mydatac + width * t0plus + (s0))).x;
	const uchar utex22 = (*(mydatac + width * t0plus + s0plus)).x;

	const float tex00 = (float)utex00 / 255.f;
	const float tex10 = (float)utex10 / 255.f;
	const float tex20 = (float)utex20 / 255.f;

	const float tex01 = (float)utex01 / 255.f;
	const float tex21 = (float)utex21 / 255.f;

	const float tex02 = (float)utex02 / 255.f;
	const float tex12 = (float)utex12 / 255.f;
	const float tex22 = (float)utex22 / 255.f;

	const float Gx = tex00 - tex20 + 2.0f * tex01 - 2.0f * tex21 + tex02 - tex22;
	const float Gy = tex00 + 2.0f * tex10 + tex20 - tex02 - 2.0f * tex12 - tex22;
	const float3 n = make_float3(Gx, Gy, 1.f);

	return n;
}

inline float3 TextureData_SampleNormalFromBump_half4(__global half const* mydatah, int width, int height, int t0, int s0)
{
	int t0minus = clamp(t0 - 1, 0, height - 1);
	int t0plus = clamp(t0 + 1, 0, height - 1);
	int s0minus = clamp(s0 - 1, 0, width - 1);
	int s0plus = clamp(s0 + 1, 0, width - 1);

	const float tex00 = vload_half4(width * t0minus + s0minus, mydatah).x;
	const float tex10 = vload_half4(width * t0minus + (s0), mydatah).x;
	const float tex20 = vload_half4(width * t0minus + s0plus, mydatah).x;

	const float tex01 = vload_half4(width * (t0)+s0minus, mydatah).x;
	const float tex21 = vload_half4(width * (t0)+s0plus, mydatah).x;

	const float tex02 = vload_half4(width * t0plus + s0minus, mydatah).x;
	const float tex12 = vload_half4(width * t0plus + (s0), mydatah).x;
	const float tex22 = vload_half4(width * t0plus + s0plus, mydatah).x;

	const float Gx = tex00 - tex20 + 2.0f * tex01 - 2.0f * tex21 + tex02 - tex22;
	const float Gy = tex00 + 2.0f * tex10 + tex20 - tex02 - 2.0f * tex12 - tex22;
	const float3 n = make_float3(Gx, Gy, 1.f);

	return n;
}

inline float3 TextureData_SampleNormalFromBump_float4(__global float4 const* mydataf, int width, int height, int t0, int s0)
{
	int t0minus = clamp(t0 - 1, 0, height - 1);
	int t0plus = clamp(t0 + 1, 0, height - 1);
	int s0minus = clamp(s0 - 1, 0, width - 1);
	int s0plus = clamp(s0 + 1, 0, width - 1);

	const float tex00 = (*(mydataf + width * t0minus + s0minus)).x;
	const float tex10 = (*(mydataf + width * t0minus + (s0))).x;
	const float tex20 = (*(mydataf + width * t0minus + s0plus)).x;

	const float tex01 = (*(mydataf + width * (t0)+s0minus)).x;
	const float tex21 = (*(mydataf + width * (t0)+s0plus)).x;

	const float tex02 = (*(mydataf + width * t0plus + s0minus)).x;
	const float tex12 = (*(mydataf + width * t0plus + (s0))).x;
	const float tex22 = (*(mydataf + width * t0plus + s0plus)).x;

	const float Gx = tex00 - tex20 + 2.0f * tex01 - 2.0f * tex21 + tex02 - tex22;
	const float Gy = tex00 + 2.0f * tex10 + tex20 - tex02 - 2.0f * tex12 - tex22;
	const float3 n = make_float3(Gx, Gy, 1.f);

	return n;
}

/// Sample 2D texture
inline
float3 Texture_SampleBump(float2 uv, TEXTURE_ARG_LIST_IDX(texidx))
{
    // Get width and height
    int width = textures[texidx].w;
    int height = textures[texidx].h;

    // Find the origin of the data in the pool
    __global char const* mydata = texturedata + textures[texidx].dataoffset;

    // Handle UV wrap
    // TODO: need UV mode support
    uv -= floor(uv);

    // Reverse Y:
    // it is needed as textures are loaded with Y axis going top to down
    // and our axis goes from down to top
    uv.y = 1.f - uv.y;

    // Calculate integer coordinates
    int s0 = clamp((int)floor(uv.x * width), 0, width - 1);
    int t0 = clamp((int)floor(uv.y * height), 0, height - 1);

	int s1 = clamp(s0 + 1, 0, width - 1);
	int t1 = clamp(t0 + 1, 0, height - 1);

	// Calculate weights for linear filtering
	float wx = uv.x * width - floor(uv.x * width);
	float wy = uv.y * height - floor(uv.y * height);

    switch (textures[texidx].fmt)
    {
    case RGBA32:
    {
        __global float3 const* mydataf = (__global float3 const*)mydata;

		float3 n00 = TextureData_SampleNormalFromBump_float4(mydataf, width, height, t0, s0);
		float3 n01 = TextureData_SampleNormalFromBump_float4(mydataf, width, height, t0, s1);
		float3 n10 = TextureData_SampleNormalFromBump_float4(mydataf, width, height, t1, s0);
		float3 n11 = TextureData_SampleNormalFromBump_float4(mydataf, width, height, t1, s1);

		float3 n = lerp3(lerp3(n00, n01, wx), lerp3(n10, n11, wx), wy);

		return 0.5f * normalize(n) + make_float3(0.5f, 0.5f, 0.5f);
    }

    case RGBA16:
    {
        __global half const* mydatah = (__global half const*)mydata;

		float3 n00 = TextureData_SampleNormalFromBump_half4(mydatah, width, height, t0, s0);
		float3 n01 = TextureData_SampleNormalFromBump_half4(mydatah, width, height, t0, s1);
		float3 n10 = TextureData_SampleNormalFromBump_half4(mydatah, width, height, t1, s0);
		float3 n11 = TextureData_SampleNormalFromBump_half4(mydatah, width, height, t1, s1);

		float3 n = lerp3(lerp3(n00, n01, wx), lerp3(n10, n11, wx), wy);

		return 0.5f * normalize(n) + make_float3(0.5f, 0.5f, 0.5f);
    }

    case RGBA8:
    {
        __global uchar4 const* mydatac = (__global uchar4 const*)mydata;

		float3 n00 = TextureData_SampleNormalFromBump_uchar4(mydatac, width, height, t0, s0);
		float3 n01 = TextureData_SampleNormalFromBump_uchar4(mydatac, width, height, t0, s1);
		float3 n10 = TextureData_SampleNormalFromBump_uchar4(mydatac, width, height, t1, s0);
		float3 n11 = TextureData_SampleNormalFromBump_uchar4(mydatac, width, height, t1, s1);

		float3 n = lerp3(lerp3(n00, n01, wx), lerp3(n10, n11, wx), wy);

		return 0.5f * normalize(n) + make_float3(0.5f, 0.5f, 0.5f);
    }

    default:
    {
        return make_float3(0.f, 0.f, 0.f);
    }
    }
}



#endif // TEXTURE_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef SAMPLING_CL
#define SAMPLING_CL


#define SAMPLE_DIMS_PER_BOUNCE 300
#define SAMPLE_DIM_CAMERA_OFFSET 1
#define SAMPLE_DIM_SURFACE_OFFSET 5
#define SAMPLE_DIM_VOLUME_APPLY_OFFSET 101
#define SAMPLE_DIM_VOLUME_EVALUATE_OFFSET 201
#define SAMPLE_DIM_IMG_PLANE_EVALUATE_OFFSET 401

typedef struct
{
    uint seq;
    uint s0;
    uint s1;
    uint s2;
} SobolSampler;

typedef struct _Sampler
{
    uint index;
    uint dimension;
    uint scramble;
    uint padding;
} Sampler;

#if SAMPLER == SOBOL
#define SAMPLER_ARG_LIST __global uint const* sobol_mat
#define SAMPLER_ARGS sobol_mat
#elif SAMPLER == RANDOM
#define SAMPLER_ARG_LIST int unused
#define SAMPLER_ARGS 0
#elif SAMPLER == CMJ
#define SAMPLER_ARG_LIST int unused
#define SAMPLER_ARGS 0
#endif

/**
    Sobol sampler
**/
#define MATSIZE 52

// The code is taken from: http://gruenschloss.org/sobol/kuo-2d-proj-single-precision.zip
// 
float SobolSampler_Sample1D(Sampler* sampler, __global uint const* mat)
{
    uint result = sampler->scramble;
    uint index = sampler->index;
    for (uint i = sampler->dimension * MATSIZE; index;  index >>= 1, ++i)
    {
        if (index & 1)
            result ^= mat[i];
    }

    return result * (1.f / (1UL << 32));
}

/**
    Random sampler
**/

/// Hash function
uint WangHash(uint seed)
{
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
}

/// Return random unsigned
uint UniformSampler_SampleUint(Sampler* sampler)
{
    sampler->index = WangHash(1664525U * sampler->index + 1013904223U);
    return sampler->index;
}

/// Return random float
float UniformSampler_Sample1D(Sampler* sampler)
{
    return ((float)UniformSampler_SampleUint(sampler)) / 0xffffffffU;
}


/**
    Correllated multi-jittered 
**/

uint permute(uint i, uint l, uint p)
{
    unsigned w = l - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;

    do
    {
        i ^= p;
        i *= 0xe170893d;
        i ^= p >> 16;
        i ^= (i & w) >> 4;
        i ^= p >> 8;
        i *= 0x0929eb3f;
        i ^= p >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | p >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;
    } while (i >= l);
    return (i + p) % l;
}

float randfloat(uint i, uint p)
{
    i ^= p;
    i ^= i >> 17;
    i ^= i >> 10;
    i *= 0xb36534e5;
    i ^= i >> 12;
    i ^= i >> 21;
    i *= 0x93fc4795;
    i ^= 0xdf6e307f;
    i ^= i >> 17;
    i *= 1 | p >> 18;
    return i * (1.0f / 4294967808.0f);
}

float2 cmj(int s, int n, int p)
{
    int sx = permute(s % n, n, p * 0xa511e9b3);
    int sy = permute(s / n, n, p * 0x63d83595);
    float jx = randfloat(s, p * 0xa399d265);
    float jy = randfloat(s, p * 0x711ad6a5);

    return make_float2((s % n + (sy + jx) / n) / n,
        (s / n + (sx + jy) / n) / n);
}

float2 CmjSampler_Sample2D(Sampler* sampler)
{
    int idx = permute(sampler->index, CMJ_DIM * CMJ_DIM, 0xa399d265 * sampler->dimension * sampler->scramble);
    return cmj(idx, CMJ_DIM, sampler->dimension * sampler->scramble);
}

#if SAMPLER == SOBOL
void Sampler_Init(Sampler* sampler, uint index, uint start_dimension, uint scramble)
{
    sampler->index = index;
    sampler->scramble = scramble;
    sampler->dimension = start_dimension;
}
#elif SAMPLER == RANDOM
void Sampler_Init(Sampler* sampler, uint seed)
{
    sampler->index = seed;
    sampler->scramble = 0;
    sampler->dimension = 0;
}
#elif SAMPLER == CMJ
void Sampler_Init(Sampler* sampler, uint index, uint dimension, uint scramble)
{
    sampler->index = index;
    sampler->scramble = scramble;
    sampler->dimension = dimension;
}
#endif


float2 Sampler_Sample2D(Sampler* sampler, SAMPLER_ARG_LIST)
{
#if SAMPLER == SOBOL
    float2 sample;
    sample.x = SobolSampler_Sample1D(sampler, SAMPLER_ARGS);
    ++(sampler->dimension);
    sample.y = SobolSampler_Sample1D(sampler, SAMPLER_ARGS);
    ++(sampler->dimension);
    return sample;
#elif SAMPLER == RANDOM
    float2 sample;
    sample.x = UniformSampler_Sample1D(sampler);
    sample.y = UniformSampler_Sample1D(sampler);
    return sample;
#elif SAMPLER == CMJ
    float2 sample;
    sample = CmjSampler_Sample2D(sampler);
    ++(sampler->dimension);
    return sample;
#endif
}

float Sampler_Sample1D(Sampler* sampler, SAMPLER_ARG_LIST)
{
#if SAMPLER == SOBOL
    float sample = SobolSampler_Sample1D(sampler, SAMPLER_ARGS);
    ++(sampler->dimension);
    return sample;
#elif SAMPLER == RANDOM
    return UniformSampler_Sample1D(sampler);
#elif SAMPLER == CMJ
    float2 sample;
    sample = CmjSampler_Sample2D(sampler);
    ++(sampler->dimension);
    return sample.x;
#endif
}

/// Sample hemisphere with cos weight
float3 Sample_MapToHemisphere(
                        // Sample
                        float2 sample,
                        // Hemisphere normal
                        float3 n,
                        // Cos power
                        float e
                        )
{
    // Construct basis
    float3 u = GetOrthoVector(n);
    float3 v = cross(u, n);
    u = cross(n, v);
    
    // Calculate 2D sample
    float r1 = sample.x;
    float r2 = sample.y;
    
    // Transform to spherical coordinates
    float sinpsi = sin(2*PI*r1);
    float cospsi = cos(2*PI*r1);
    float costheta = pow(1.f - r2, 1.f/(e + 1.f));
    float sintheta = sqrt(1.f - costheta * costheta);
    
    // Return the result
    return normalize(u * sintheta * cospsi + v * sintheta * sinpsi + n * costheta);
}

float2 Sample_MapToDisk(
    // Sample
    float2 sample
    )
{
    float r = native_sqrt(sample.x); 
    float theta = 2 * PI * sample.y;
    return make_float2(r * native_cos(theta), r * native_sin(theta));
}

float2 Sample_MapToDiskConcentric(
    // Sample
    float2 sample
    )
{
    float2 offset = 2.f * sample - make_float2(1.f, 1.f);

    if (offset.x == 0 && offset.y == 0) return 0.f;

    float theta, r;

    if (fabs(offset.x) > fabs(offset.y)) 
    {
        r = offset.x;
        theta = PI / 4.f * (offset.y / offset.x);
    }
    else 
    {
        r = offset.y;
        theta = PI / 2.f * ( 1.f - 0.5f * (offset.x / offset.y));
    }
    
    return make_float2(r * native_cos(theta), r * native_sin(theta));
}

/// Sample hemisphere with cos weight
float3 Sample_MapToSphere(
                        // Sample
                        float2 sample
                        )
{
    float z = 1.f - 2.f * sample.x;
    float r = native_sqrt(max(0.f, 1.f - z*z));
    float phi = 2.f * PI * sample.y;
    float x = cos(phi);
    float y = sin(phi);
    
    // Return the result
    return make_float3(x,y,z);
}

float2 Sample_MapToPolygon(int n, float2 sample, float sample1)
{
    float theta = 2.f * PI / n;
    int edge = clamp((int)(sample1 * n), 0, n - 1);
    float t = native_sqrt(sample.x);
    float u = 1.f - t;
    float v = t * sample.y;
    float2 v1 = make_float2(native_cos(theta * edge), native_sin(theta * edge));
    float2 v2 = make_float2(native_cos(theta * (edge + 1)), native_sin(theta * (edge + 1)));
    return u*v1 + v*v2;;
}

/// Power heuristic for multiple importance sampling
float PowerHeuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f*f) / (f*f + g*g);
}

/// Balance heuristic for multiple importance sampling
float BalanceHeuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f) / (f + g);
}

int lower_bound(GLOBAL float const* values, int n, float value)
{
    int count = n;
    int b = 0;
    int it = 0;
    int step = 0;

    while (count > 0)
    {
        it = b;
        step = count / 2;
        it += step;
        if (values[it] < value)
        {
            b = ++it;
            count -= step + 1;
        }
        else
        {
            count = step;
        }
    }

    return b;
}
/// Sample 1D distribution
float Distribution1D_Sample(float s, GLOBAL int const* data, float* pdf)
{
    int num_segments = data[0];

    GLOBAL float const* cdf_data = (GLOBAL float const*)&data[1];
    GLOBAL float const* pdf_data = cdf_data + num_segments + 1;

    int segment_idx = max(lower_bound(cdf_data, num_segments + 1, s), 1);

    // Find lerp coefficient
    float du = (s - cdf_data[segment_idx - 1]) / (cdf_data[segment_idx] - cdf_data[segment_idx - 1]);

    // Calc pdf
    *pdf = pdf_data[segment_idx - 1];

    return (segment_idx - 1 + du) / num_segments;;
}

/// Sample 1D distribution
int Distribution1D_SampleDiscrete(float s, GLOBAL int const* data, float* pdf)
{
    int num_segments = data[0];

    GLOBAL float const* cdf_data = (GLOBAL float const*)&data[1];
    GLOBAL float const* pdf_data = cdf_data + num_segments + 1;

    int segment_idx = max(lower_bound(cdf_data, num_segments + 1, s), 1);

    // Find lerp coefficient
    float du = (s - cdf_data[segment_idx - 1]) / (cdf_data[segment_idx] - cdf_data[segment_idx - 1]);

    // Calc pdf
    *pdf = pdf_data[segment_idx - 1] / num_segments;

    return segment_idx - 1;
}

/// PDF of  1D distribution
float Distribution1D_GetPdf(float s, GLOBAL int const* data)
{
    int num_segments = data[0];
    GLOBAL float const* cdf_data = (GLOBAL float const*)&data[1];
    GLOBAL float const* pdf_data = cdf_data + num_segments + 1;

    int segment_idx = max(lower_bound(cdf_data, num_segments + 1, s), 1);

    // Calc pdf
    return pdf_data[segment_idx - 1];
}

/// PDF of  1D distribution
float Distribution1D_GetPdfDiscreet(int d, GLOBAL int const* data)
{
    int num_segments = data[0];
    GLOBAL float const* cdf_data = (GLOBAL float const*)&data[1];
    GLOBAL float const* pdf_data = cdf_data + num_segments + 1;

    // Calc pdf
    return pdf_data[d] / num_segments;
}



#endif // SAMPLING_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef BXDF_CL
#define BXDF_CL
/**********************************************************************
 * Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 ********************************************************************/
#ifndef BXDF_FLAGS_CL
#define BXDF_FLAGS_CL

#define DENOM_EPS 1e-8f
#define ROUGHNESS_EPS 0.0001f

enum BxdfFlags
{
    kBxdfFlagsSingular = (1 << 0),
    kBxdfFlagsBrdf = (1 << 1),
    kBxdfFlagsEmissive = (1 << 2),
    kBxdfFlagsTransparency = (1 << 3),
    kBxdfFlagsDiffuse = (1 << 4),

    //Used to mask value from bxdf_flags
    kBxdfFlagsAll = (kBxdfFlagsSingular | kBxdfFlagsBrdf | kBxdfFlagsEmissive | kBxdfFlagsTransparency | kBxdfFlagsDiffuse)
};

enum BxdfUberV2SampledComponent
{
    kBxdfUberV2SampleTransparency = 0,
    kBxdfUberV2SampleCoating = 1,
    kBxdfUberV2SampleReflection = 2,
    kBxdfUberV2SampleRefraction = 3,
    kBxdfUberV2SampleDiffuse = 4
};

/// Returns BxDF flags. Flags stored in first byte of bxdf_flags
int Bxdf_GetFlags(DifferentialGeometry const* dg)
{
    return (dg->mat.flags & kBxdfFlagsAll);
}

/// Sets BxDF flags. Flags stored in first byte of bxdf_flags
void Bxdf_SetFlags(DifferentialGeometry *dg, int flags)
{
    dg->mat.flags &= 0xffffff00; //Reset flags
    dg->mat.flags |= flags; //Set new flags
}

/// Return BxDF sampled component. Sampled component stored in second byte of bxdf_flags
int Bxdf_UberV2_GetSampledComponent(DifferentialGeometry const* dg)
{
    return (dg->mat.flags >> 8) & 0xff;
}

/// Sets BxDF sampled component. Sampled component stored in second byte of bxdf_flags
void Bxdf_UberV2_SetSampledComponent(DifferentialGeometry *dg, int sampledComponent)
{
    dg->mat.flags &= 0xffff00ff; //Reset sampled component
    dg->mat.flags |= (sampledComponent << 8); //Set new component
}

#endif


/// Schlick's approximation of Fresnel equtions
float SchlickFresnel(float eta, float ndotw)
{
    const float f = ((1.f - eta) / (1.f + eta)) * ((1.f - eta) / (1.f + eta));
    const float m = 1.f - fabs(ndotw);
    const float m2 = m*m;
    return f + (1.f - f) * m2 * m2 * m;
}

/// Full Fresnel equations
float FresnelDielectric(float etai, float etat, float ndotwi, float ndotwt)
{
    // Parallel and perpendicular polarization
    float rparl = ((etat * ndotwi) - (etai * ndotwt)) / ((etat * ndotwi) + (etai * ndotwt));
    float rperp = ((etai * ndotwi) - (etat * ndotwt)) / ((etai * ndotwi) + (etat * ndotwt));
    return (rparl*rparl + rperp*rperp) * 0.5f;
}
#ifndef BXDF_UBERV2_CL
#define BXDF_UBERV2_CL
#ifndef INPUTMAPS_CL
#define INPUTMAPS_CL

float4 ReadInputMap635(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[0].float_value.value, 0.0f))
	);
}
float4 ReadInputMap636(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[1].float_value.value, 0.0f))
	);
}
float4 ReadInputMap637(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[2].float_value.value, 0.0f))
	);
}
float4 ReadInputMap638(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[3].float_value.value, 0.0f))
	);
}
float4 ReadInputMap640(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[4].int_values.idx % 10) / 10.0f, (input_map_values[4].int_values.idx % 20) / 20.0f, (input_map_values[4].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[4].int_values.idx))
	);
}
float4 ReadInputMap642(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[5].int_values.idx % 10) / 10.0f, (input_map_values[5].int_values.idx % 20) / 20.0f, (input_map_values[5].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[5].int_values.idx))
	);
}
float4 ReadInputMap643(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[6].float_value.value, 0.0f))
	);
}
float4 ReadInputMap645(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[7].int_values.idx % 10) / 10.0f, (input_map_values[7].int_values.idx % 20) / 20.0f, (input_map_values[7].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[7].int_values.idx))
	);
}
float4 ReadInputMap647(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[8].int_values.idx % 10) / 10.0f, (input_map_values[8].int_values.idx % 20) / 20.0f, (input_map_values[8].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[8].int_values.idx))
	);
}
float4 ReadInputMap649(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[9].int_values.idx % 10) / 10.0f, (input_map_values[9].int_values.idx % 20) / 20.0f, (input_map_values[9].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[9].int_values.idx))
	);
}
float4 ReadInputMap650(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[10].float_value.value, 0.0f))
	);
}
float4 ReadInputMap651(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[11].float_value.value, 0.0f))
	);
}
float4 ReadInputMap653(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[12].int_values.idx % 10) / 10.0f, (input_map_values[12].int_values.idx % 20) / 20.0f, (input_map_values[12].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[12].int_values.idx))
	);
}
float4 ReadInputMap654(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[13].float_value.value, 0.0f))
	);
}
float4 ReadInputMap655(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[14].float_value.value, 0.0f))
	);
}
float4 ReadInputMap656(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[15].float_value.value, 0.0f))
	);
}
float4 ReadInputMap657(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[16].float_value.value, 0.0f))
	);
}
float4 ReadInputMap658(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[17].float_value.value, 0.0f))
	);
}
float4 ReadInputMap659(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[18].float_value.value, 0.0f))
	);
}
float4 ReadInputMap660(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[19].float_value.value, 0.0f))
	);
}
float4 ReadInputMap661(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[20].float_value.value, 0.0f))
	);
}
float4 ReadInputMap662(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[21].float_value.value, 0.0f))
	);
}
float4 ReadInputMap663(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[22].float_value.value, 0.0f))
	);
}
float4 ReadInputMap664(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[23].float_value.value, 0.0f))
	);
}
float4 ReadInputMap665(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[24].float_value.value, 0.0f))
	);
}
float4 ReadInputMap667(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[25].int_values.idx % 10) / 10.0f, (input_map_values[25].int_values.idx % 20) / 20.0f, (input_map_values[25].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[25].int_values.idx))
	);
}
float4 ReadInputMap668(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[26].float_value.value, 0.0f))
	);
}
float4 ReadInputMap669(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[27].float_value.value, 0.0f))
	);
}
float4 ReadInputMap670(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[28].float_value.value, 0.0f))
	);
}
float4 ReadInputMap671(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[29].float_value.value, 0.0f))
	);
}
float4 ReadInputMap673(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[30].int_values.idx % 10) / 10.0f, (input_map_values[30].int_values.idx % 20) / 20.0f, (input_map_values[30].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[30].int_values.idx))
	);
}
float4 ReadInputMap674(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[31].float_value.value, 0.0f))
	);
}
float4 ReadInputMap675(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[32].float_value.value, 0.0f))
	);
}
float4 ReadInputMap676(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[33].float_value.value, 0.0f))
	);
}
float4 ReadInputMap677(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[34].float_value.value, 0.0f))
	);
}
float4 ReadInputMap678(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[35].float_value.value, 0.0f))
	);
}
float4 ReadInputMap680(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[36].int_values.idx % 10) / 10.0f, (input_map_values[36].int_values.idx % 20) / 20.0f, (input_map_values[36].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[36].int_values.idx))
	);
}
float4 ReadInputMap682(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[37].int_values.idx % 10) / 10.0f, (input_map_values[37].int_values.idx % 20) / 20.0f, (input_map_values[37].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[37].int_values.idx))
	);
}
float4 ReadInputMap683(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[38].float_value.value, 0.0f))
	);
}
float4 ReadInputMap684(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[39].float_value.value, 0.0f))
	);
}
float4 ReadInputMap686(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[40].int_values.idx % 10) / 10.0f, (input_map_values[40].int_values.idx % 20) / 20.0f, (input_map_values[40].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[40].int_values.idx))
	);
}
float4 ReadInputMap687(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[41].float_value.value, 0.0f))
	);
}
float4 ReadInputMap688(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[42].float_value.value, 0.0f))
	);
}
float4 ReadInputMap689(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[43].float_value.value, 0.0f))
	);
}
float4 ReadInputMap690(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[44].float_value.value, 0.0f))
	);
}
float4 ReadInputMap691(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[45].float_value.value, 0.0f))
	);
}
float4 ReadInputMap692(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[46].float_value.value, 0.0f))
	);
}
float4 ReadInputMap693(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[47].float_value.value, 0.0f))
	);
}
float4 ReadInputMap695(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[48].int_values.idx % 10) / 10.0f, (input_map_values[48].int_values.idx % 20) / 20.0f, (input_map_values[48].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[48].int_values.idx))
	);
}
float4 ReadInputMap696(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[49].float_value.value, 0.0f))
	);
}
float4 ReadInputMap697(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[50].float_value.value, 0.0f))
	);
}
float4 ReadInputMap698(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[51].float_value.value, 0.0f))
	);
}
float4 ReadInputMap699(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[52].float_value.value, 0.0f))
	);
}
float4 ReadInputMap700(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[53].float_value.value, 0.0f))
	);
}
float4 ReadInputMap702(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	(float4)((input_map_values[54].int_values.idx % 10) / 10.0f, (input_map_values[54].int_values.idx % 20) / 20.0f, (input_map_values[54].int_values.idx % 30), 0.0f) * Texture_Sample2D(dg->uv, TEXTURE_ARGS_IDX(input_map_values[54].int_values.idx))
	);
}
float4 ReadInputMap703(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[55].float_value.value, 0.0f))
	);
}
float4 ReadInputMap704(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[56].float_value.value, 0.0f))
	);
}
float4 ReadInputMap705(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[57].float_value.value, 0.0f))
	);
}
float4 ReadInputMap706(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[58].float_value.value, 0.0f))
	);
}
float4 ReadInputMap707(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[59].float_value.value, 0.0f))
	);
}
float4 ReadInputMap708(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[60].float_value.value, 0.0f))
	);
}
float4 ReadInputMap709(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[61].float_value.value, 0.0f))
	);
}
float4 ReadInputMap710(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return (float4)(
	((float4)(input_map_values[62].float_value.value, 0.0f))
	);
}
float4 GetInputMapFloat4(uint input_id, DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	switch(input_id)
	{
		case 635: return ReadInputMap635(dg, input_map_values, TEXTURE_ARGS);
		case 636: return ReadInputMap636(dg, input_map_values, TEXTURE_ARGS);
		case 637: return ReadInputMap637(dg, input_map_values, TEXTURE_ARGS);
		case 638: return ReadInputMap638(dg, input_map_values, TEXTURE_ARGS);
		case 640: return ReadInputMap640(dg, input_map_values, TEXTURE_ARGS);
		case 642: return ReadInputMap642(dg, input_map_values, TEXTURE_ARGS);
		case 643: return ReadInputMap643(dg, input_map_values, TEXTURE_ARGS);
		case 645: return ReadInputMap645(dg, input_map_values, TEXTURE_ARGS);
		case 647: return ReadInputMap647(dg, input_map_values, TEXTURE_ARGS);
		case 649: return ReadInputMap649(dg, input_map_values, TEXTURE_ARGS);
		case 650: return ReadInputMap650(dg, input_map_values, TEXTURE_ARGS);
		case 651: return ReadInputMap651(dg, input_map_values, TEXTURE_ARGS);
		case 653: return ReadInputMap653(dg, input_map_values, TEXTURE_ARGS);
		case 654: return ReadInputMap654(dg, input_map_values, TEXTURE_ARGS);
		case 655: return ReadInputMap655(dg, input_map_values, TEXTURE_ARGS);
		case 656: return ReadInputMap656(dg, input_map_values, TEXTURE_ARGS);
		case 657: return ReadInputMap657(dg, input_map_values, TEXTURE_ARGS);
		case 658: return ReadInputMap658(dg, input_map_values, TEXTURE_ARGS);
		case 659: return ReadInputMap659(dg, input_map_values, TEXTURE_ARGS);
		case 660: return ReadInputMap660(dg, input_map_values, TEXTURE_ARGS);
		case 661: return ReadInputMap661(dg, input_map_values, TEXTURE_ARGS);
		case 662: return ReadInputMap662(dg, input_map_values, TEXTURE_ARGS);
		case 663: return ReadInputMap663(dg, input_map_values, TEXTURE_ARGS);
		case 664: return ReadInputMap664(dg, input_map_values, TEXTURE_ARGS);
		case 665: return ReadInputMap665(dg, input_map_values, TEXTURE_ARGS);
		case 667: return ReadInputMap667(dg, input_map_values, TEXTURE_ARGS);
		case 668: return ReadInputMap668(dg, input_map_values, TEXTURE_ARGS);
		case 669: return ReadInputMap669(dg, input_map_values, TEXTURE_ARGS);
		case 670: return ReadInputMap670(dg, input_map_values, TEXTURE_ARGS);
		case 671: return ReadInputMap671(dg, input_map_values, TEXTURE_ARGS);
		case 673: return ReadInputMap673(dg, input_map_values, TEXTURE_ARGS);
		case 674: return ReadInputMap674(dg, input_map_values, TEXTURE_ARGS);
		case 675: return ReadInputMap675(dg, input_map_values, TEXTURE_ARGS);
		case 676: return ReadInputMap676(dg, input_map_values, TEXTURE_ARGS);
		case 677: return ReadInputMap677(dg, input_map_values, TEXTURE_ARGS);
		case 678: return ReadInputMap678(dg, input_map_values, TEXTURE_ARGS);
		case 680: return ReadInputMap680(dg, input_map_values, TEXTURE_ARGS);
		case 682: return ReadInputMap682(dg, input_map_values, TEXTURE_ARGS);
		case 683: return ReadInputMap683(dg, input_map_values, TEXTURE_ARGS);
		case 684: return ReadInputMap684(dg, input_map_values, TEXTURE_ARGS);
		case 686: return ReadInputMap686(dg, input_map_values, TEXTURE_ARGS);
		case 687: return ReadInputMap687(dg, input_map_values, TEXTURE_ARGS);
		case 688: return ReadInputMap688(dg, input_map_values, TEXTURE_ARGS);
		case 689: return ReadInputMap689(dg, input_map_values, TEXTURE_ARGS);
		case 690: return ReadInputMap690(dg, input_map_values, TEXTURE_ARGS);
		case 691: return ReadInputMap691(dg, input_map_values, TEXTURE_ARGS);
		case 692: return ReadInputMap692(dg, input_map_values, TEXTURE_ARGS);
		case 693: return ReadInputMap693(dg, input_map_values, TEXTURE_ARGS);
		case 695: return ReadInputMap695(dg, input_map_values, TEXTURE_ARGS);
		case 696: return ReadInputMap696(dg, input_map_values, TEXTURE_ARGS);
		case 697: return ReadInputMap697(dg, input_map_values, TEXTURE_ARGS);
		case 698: return ReadInputMap698(dg, input_map_values, TEXTURE_ARGS);
		case 699: return ReadInputMap699(dg, input_map_values, TEXTURE_ARGS);
		case 700: return ReadInputMap700(dg, input_map_values, TEXTURE_ARGS);
		case 702: return ReadInputMap702(dg, input_map_values, TEXTURE_ARGS);
		case 703: return ReadInputMap703(dg, input_map_values, TEXTURE_ARGS);
		case 704: return ReadInputMap704(dg, input_map_values, TEXTURE_ARGS);
		case 705: return ReadInputMap705(dg, input_map_values, TEXTURE_ARGS);
		case 706: return ReadInputMap706(dg, input_map_values, TEXTURE_ARGS);
		case 707: return ReadInputMap707(dg, input_map_values, TEXTURE_ARGS);
		case 708: return ReadInputMap708(dg, input_map_values, TEXTURE_ARGS);
		case 709: return ReadInputMap709(dg, input_map_values, TEXTURE_ARGS);
		case 710: return ReadInputMap710(dg, input_map_values, TEXTURE_ARGS);
	}
	return 0.0f;
}
float GetInputMapFloat(uint input_id, DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values, TEXTURE_ARG_LIST)
{
	return GetInputMapFloat4(input_id, dg, input_map_values, TEXTURE_ARGS).x;
}
#endif



typedef struct _UberV2ShaderData
{
    float4 diffuse_color;
    float4 reflection_color;
    float4 coating_color;
    float4 refraction_color;
    float4 emission_color;
    float4 sss_absorption_color;
    float4 sss_scatter_color;
    float4 sss_subsurface_color;
    float4 shading_normal;

    float reflection_roughness;
    float reflection_anisotropy;
    float reflection_anisotropy_rotation;
    float reflection_ior;

    float reflection_metalness;
    float coating_ior;
    float refraction_roughness;
    float refraction_ior;

    float transparency;
    float sss_absorption_distance;
    float sss_scatter_distance;
    float sss_scatter_direction;

} UberV2ShaderData;

bool UberV2IsTransmissive(
    // Layers
    int layers
    )
{
    return (layers & (kTransparencyLayer | kRefractionLayer)) != 0;
}

float4 GetUberV2EmissionColor(
    // Material offset
    int offset,
    // Geometry
    DifferentialGeometry const* dg,
    // Values for input maps
    GLOBAL InputMapData const* restrict input_map_values,
    // Material attributes
    GLOBAL int const* restrict material_attributes,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return GetInputMapFloat4(material_attributes[offset+1], dg, input_map_values, TEXTURE_ARGS);
}

#ifndef BXDF_UBERV2_BRICKS
#define BXDF_UBERV2_BRICKS

// Utility functions for Uberv2
/// Calculates Fresnel for provided parameters. Swaps IORs if needed
float CalculateFresnel(
    // IORs
    float top_ior,
    float bottom_ior,
    // Angle between normal and incoming ray
    float ndotwi
)
{
    float etai = top_ior;
    float etat = bottom_ior;
    float cosi = ndotwi;

    // Revert normal and eta if needed
    if (cosi < 0.f)
    {
        float tmp = etai;
        etai = etat;
        etat = tmp;
        cosi = -cosi;
    }

    float eta = etai / etat;
    float sini2 = 1.f - cosi * cosi;
    float sint2 = eta * eta * sini2;
    float fresnel = 1.f;

    if (sint2 < 1.f)
    {
        float cost = native_sqrt(max(0.f, 1.f - sint2));
        fresnel = FresnelDielectric(etai, etat, cosi, cost);
    }

    return fresnel;
}

// Calucates Fresnel blend for two float3 values.
// F(top_ior, bottom_ior) * top_value + (1 - F(top_ior, bottom_ior) * bottom_value)
float3 Fresnel_Blend(
    // IORs
    float top_ior,
    float bottom_ior,
    // Values to blend
    float3 top_value,
    float3 bottom_value,
    // Incoming direction
    float3 wi
)
{
    float fresnel = CalculateFresnel(top_ior, bottom_ior, wi.y);
    return fresnel * top_value + (1.f - fresnel) * bottom_value;
}

// Calucates Fresnel blend for two float values.
// F(top_ior, bottom_ior) * top_value + (1 - F(top_ior, bottom_ior) * bottom_value)
float Fresnel_Blend_F(
    // IORs
    float top_ior,
    float bottom_ior,
    // Values to blend
    float top_value,
    float bottom_value,
    // Incoming direction
    float3 wi
)
{
    float fresnel = CalculateFresnel(top_ior, bottom_ior, wi.y);
    return fresnel * top_value + (1.f - fresnel) * bottom_value;
}

// Diffuse layer
float3 UberV2_Lambert_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return shader_data->diffuse_color.xyz / PI;
}

float UberV2_Lambert_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return fabs(wo.y) / PI;
}

/// Lambert BRDF sampling
float3 UberV2_Lambert_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const float3 kd = UberV2_Lambert_Evaluate(shader_data, wi, *wo, TEXTURE_ARGS);

    *wo = Sample_MapToHemisphere(sample, make_float3(0.f, 1.f, 0.f), 1.f);

    *pdf = fabs((*wo).y) / PI;

    return kd;
}

// Reflection/Coating
/*
Microfacet GGX
*/
// Distribution fucntion
float UberV2_MicrofacetDistribution_GGX_D(float roughness, float3 m)
{
    float ndotm = fabs(m.y);
    float ndotm2 = ndotm * ndotm;
    float sinmn = native_sqrt(1.f - clamp(ndotm * ndotm, 0.f, 1.f));
    float tanmn = ndotm > DENOM_EPS ? sinmn / ndotm : 0.f;
    float a2 = roughness * roughness;
    float denom = (PI * ndotm2 * ndotm2 * (a2 + tanmn * tanmn) * (a2 + tanmn * tanmn));
    return denom > DENOM_EPS ? (a2 / denom) : 1.f;
}

// PDF of the given direction
float UberV2_MicrofacetDistribution_GGX_GetPdf(
    // Halfway vector
    float3 m,
    // Rougness
    float roughness,
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    float mpdf = UberV2_MicrofacetDistribution_GGX_D(roughness, m) * fabs(m.y);
    // See Humphreys and Pharr for derivation
    float denom = (4.f * fabs(dot(wo, m)));

    return denom > DENOM_EPS ? mpdf / denom : 0.f;
}

// Sample the distribution
void UberV2_MicrofacetDistribution_GGX_SampleNormal(
    // Roughness
    float roughness,
    // Differential geometry
    UberV2ShaderData const* shader_data,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wh
)
{
    float r1 = sample.x;
    float r2 = sample.y;

    // Sample halfway vector first, then reflect wi around that
    float theta = atan2(roughness * native_sqrt(r1), native_sqrt(1.f - r1));
    float costheta = native_cos(theta);
    float sintheta = native_sin(theta);

    // phi = 2*PI*ksi2
    float phi = 2.f * PI * r2;
    float cosphi = native_cos(phi);
    float sinphi = native_sin(phi);

    // Calculate wh
    *wh = make_float3(sintheta * cosphi, costheta, sintheta * sinphi);
}

//
float UberV2_MicrofacetDistribution_GGX_G1(float roughness, float3 v, float3 m)
{
    float ndotv = fabs(v.y);
    float mdotv = fabs(dot(m, v));

    float sinnv = native_sqrt(1.f - clamp(ndotv * ndotv, 0.f, 1.f));
    float tannv = ndotv > DENOM_EPS ? sinnv / ndotv : 0.f;
    float a2 = roughness * roughness;
    return 2.f / (1.f + native_sqrt(1.f + a2 * tannv * tannv));
}

// Shadowing function also depends on microfacet distribution
float UberV2_MicrofacetDistribution_GGX_G(float roughness, float3 wi, float3 wo, float3 wh)
{
    return UberV2_MicrofacetDistribution_GGX_G1(roughness, wi, wh) * UberV2_MicrofacetDistribution_GGX_G1(roughness, wo, wh);
}

float3 UberV2_MicrofacetGGX_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST,
    float3 ks
)
{
    // Incident and reflected zenith angles
    float costhetao = fabs(wo.y);
    float costhetai = fabs(wi.y);

    // Calc halfway vector
    float3 wh = normalize(wi + wo);

    float denom = (4.f * costhetao * costhetai);

    return denom > DENOM_EPS ? ks * UberV2_MicrofacetDistribution_GGX_G(shader_data->reflection_roughness, wi, wo, wh) * UberV2_MicrofacetDistribution_GGX_D(shader_data->reflection_roughness, wh) / denom : 0.f;
}


float UberV2_MicrofacetGGX_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    float3 wh = normalize(wo + wi);

    return UberV2_MicrofacetDistribution_GGX_GetPdf(wh, shader_data->reflection_roughness, shader_data, wi, wo, TEXTURE_ARGS);
}

float3 UberV2_MicrofacetGGX_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf,
    float3 ks
)
{
    float3 wh;
    UberV2_MicrofacetDistribution_GGX_SampleNormal(shader_data->reflection_roughness, shader_data, TEXTURE_ARGS, sample, &wh);

    *wo = -wi + 2.f*fabs(dot(wi, wh)) * wh;

    *pdf = UberV2_MicrofacetDistribution_GGX_GetPdf(wh, shader_data->reflection_roughness, shader_data, wi, *wo, TEXTURE_ARGS);

    return UberV2_MicrofacetGGX_Evaluate(shader_data, wi, *wo, TEXTURE_ARGS, ks);
}

/*
Ideal reflection BRDF
*/
float3 UberV2_IdealReflect_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

float UberV2_IdealReflect_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

float3 UberV2_IdealReflect_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf,
    float3 ks)
{
    // Mirror reflect wi
    *wo = normalize(make_float3(-wi.x, wi.y, -wi.z));

    // PDF is infinite at that point, but deltas are going to cancel out while evaluating
    // so set it to 1.f
    *pdf = 1.f;

    float coswo = fabs((*wo).y);

    // Return reflectance value
    return coswo > DENOM_EPS ? (ks * (1.f / coswo)) : 0.f;
}

// Refraction
/*
Ideal refraction BTDF
*/

float3 UberV2_IdealRefract_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

float UberV2_IdealRefract_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

float3 UberV2_IdealRefract_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const float3 ks = shader_data->refraction_color.xyz;

    float etai = 1.f;
    float etat = shader_data->refraction_ior;
    float cosi = wi.y;

    bool entering = cosi > 0.f;

    // Revert normal and eta if needed
    if (!entering)
    {
        float tmp = etai;
        etai = etat;
        etat = tmp;
    }

    float eta = etai / etat;
    float sini2 = 1.f - cosi * cosi;
    float sint2 = eta * eta * sini2;

    if (sint2 >= 1.f)
    {
        *pdf = 0.f;
        return 0.f;
    }

    float cost = native_sqrt(max(0.f, 1.f - sint2));

    // Transmitted ray
    *wo = normalize(make_float3(eta * -wi.x, entering ? -cost : cost, eta * -wi.z));

    *pdf = 1.f;

    return cost > DENOM_EPS ? (eta * eta * ks / cost) : 0.f;
}


float3 UberV2_MicrofacetRefractionGGX_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const float3 ks = shader_data->refraction_color.xyz;
    const float roughness = max(shader_data->refraction_roughness, ROUGHNESS_EPS);

    float ndotwi = wi.y;
    float ndotwo = wo.y;

    if (ndotwi * ndotwo >= 0.f)
    {
        return 0.f;
    }

    float etai = 1.f;
    float etat = shader_data->refraction_ior;

    // Revert normal and eta if needed
    if (ndotwi < 0.f)
    {
        float tmp = etai;
        etai = etat;
        etat = tmp;
    }

    // Calc halfway vector
    float3 ht = -(etai * wi + etat * wo);
    float3 wh = normalize(ht);

    float widotwh = fabs(dot(wh, wi));
    float wodotwh = fabs(dot(wh, wo));

    float denom = dot(ht, ht);
    denom *= (fabs(ndotwi) * fabs(ndotwo));

    return denom > DENOM_EPS ? (ks * (widotwh * wodotwh)  * (etat)* (etat)*
        UberV2_MicrofacetDistribution_GGX_G(roughness, wi, wo, wh) * UberV2_MicrofacetDistribution_GGX_D(roughness, wh) / denom) : 0.f;
}

float UberV2_MicrofacetRefractionGGX_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const float roughness = max(shader_data->refraction_roughness, ROUGHNESS_EPS);

    float ndotwi = wi.y;
    float ndotwo = wo.y;

    if (ndotwi * ndotwo >= 0.f)
    {
        return 0.f;
    }

    float etai = 1.f;
    float etat = shader_data->refraction_ior;

    // Revert normal and eta if needed
    if (ndotwi < 0.f)
    {
        float tmp = etai;
        etai = etat;
        etat = tmp;
    }

    // Calc halfway vector
    float3 ht = -(etai * wi + etat * wo);

    float3 wh = normalize(ht);

    float wodotwh = fabs(dot(wo, wh));

    float whpdf = UberV2_MicrofacetDistribution_GGX_D(roughness, wh) * fabs(wh.y);

    float whwo = wodotwh * etat * etat;

    float denom = dot(ht, ht);

    return denom > DENOM_EPS ? whpdf * whwo / denom : 0.f;
}

float3 UberV2_MicrofacetRefractionGGX_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const float3 ks = shader_data->refraction_color.xyz;
    const float roughness = max(shader_data->refraction_roughness, ROUGHNESS_EPS);

    float ndotwi = wi.y;

    if (ndotwi == 0.f)
    {
        *pdf = 0.f;
        return 0.f;
    }

    float etai = 1.f;
    float etat = shader_data->refraction_ior;
    float s = 1.f;

    // Revert normal and eta if needed
    if (ndotwi < 0.f)
    {
        float tmp = etai;
        etai = etat;
        etat = tmp;
        s = -s;
    }

    float3 wh;
    UberV2_MicrofacetDistribution_GGX_SampleNormal(roughness, shader_data, TEXTURE_ARGS, sample, &wh);

    float c = dot(wi, wh);
    float eta = etai / etat;

    float d = 1 + eta * (c * c - 1);

    if (d <= 0.f)
    {
        *pdf = 0.f;
        return 0.f;
    }

    *wo = normalize((eta * c - s * native_sqrt(d)) * wh - eta * wi);

    *pdf = UberV2_MicrofacetRefractionGGX_GetPdf(shader_data, wi, *wo, TEXTURE_ARGS);

    return UberV2_MicrofacetRefractionGGX_Evaluate(shader_data, wi, *wo, TEXTURE_ARGS);
}

float3 UberV2_Passthrough_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    *wo = -wi;
    float coswo = fabs((*wo).y);

    // PDF is infinite at that point, but deltas are going to cancel out while evaluating
    // so set it to 1.f
    *pdf = 1.f;

    return coswo > 1e-5f ? (1.f / coswo) : 0.f;
}

float3 UberV2_Reflection_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const bool is_singular = (shader_data->reflection_roughness < ROUGHNESS_EPS);
    const float metalness = shader_data->reflection_metalness;

    const float3 ks = shader_data->reflection_color.xyz;

    float3 color = mix((float3)(1.0f, 1.0f, 1.0f), ks, metalness);

    return is_singular ?
        UberV2_IdealReflect_Evaluate(shader_data, wi, wo, TEXTURE_ARGS) :
        UberV2_MicrofacetGGX_Evaluate(shader_data, wi, wo, TEXTURE_ARGS, color);
}

float3 UberV2_Refraction_Evaluate(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const bool is_singular = (shader_data->refraction_roughness < ROUGHNESS_EPS);

    return is_singular ?
        UberV2_IdealRefract_Evaluate(shader_data, wi, wo, TEXTURE_ARGS) :
        UberV2_MicrofacetRefractionGGX_Evaluate(shader_data, wi, wo, TEXTURE_ARGS);
}

float UberV2_Reflection_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const bool is_singular = (shader_data->reflection_roughness < ROUGHNESS_EPS);

    return is_singular ?
        UberV2_IdealReflect_GetPdf(shader_data, wi, wo, TEXTURE_ARGS) :
        UberV2_MicrofacetGGX_GetPdf(shader_data, wi, wo, TEXTURE_ARGS);
}

float UberV2_Refraction_GetPdf(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Outgoing direction
    float3 wo,
    // Texture args
    TEXTURE_ARG_LIST
)
{
    const bool is_singular = (shader_data->refraction_roughness < ROUGHNESS_EPS);

    return is_singular ?
        UberV2_IdealRefract_GetPdf(shader_data, wi, wo, TEXTURE_ARGS) :
        UberV2_MicrofacetRefractionGGX_GetPdf(shader_data, wi, wo, TEXTURE_ARGS);
}

float3 UberV2_Reflection_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const float3 ks = shader_data->reflection_color.xyz;
    const bool is_singular = (shader_data->reflection_roughness < ROUGHNESS_EPS);
    const float metalness = shader_data->reflection_metalness;

    float3 color = mix((float3)(1.0f, 1.0f, 1.0f), ks, metalness);

    return is_singular ?
        UberV2_IdealReflect_Sample(shader_data, wi, TEXTURE_ARGS, wo, pdf, color) :
        UberV2_MicrofacetGGX_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf, color);
}

float3 UberV2_Refraction_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const bool is_singular = (shader_data->refraction_roughness < ROUGHNESS_EPS);

    return is_singular ?
        UberV2_IdealRefract_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf) :
        UberV2_MicrofacetRefractionGGX_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
}

float3 UberV2_Coating_Sample(
    // Preprocessed shader input data
    UberV2ShaderData const* shader_data,
    // Incoming direction
    float3 wi,
    // Texture args
    TEXTURE_ARG_LIST,
    // Outgoing  direction
    float3* wo,
    // PDF at wo
    float* pdf
)
{
    const float3 ks = shader_data->coating_color.xyz;

    return UberV2_IdealReflect_Sample(shader_data, wi, TEXTURE_ARGS, wo, pdf, ks);
}

#endif


void UberV2_ApplyShadingNormal(
    // Geometry
    DifferentialGeometry* dg,
    // Prepared UberV2 shader inputs
    UberV2ShaderData const* shader_data
)
{
    const int layers = dg->mat.layers;

    if ((layers & kShadingNormalLayer) == kShadingNormalLayer)
    {
        dg->n = normalize(shader_data->shading_normal.z * dg->n + shader_data->shading_normal.x * dg->dpdu + shader_data->shading_normal.y * dg->dpdv);
        dg->dpdv = normalize(cross(dg->n, dg->dpdu));
        dg->dpdu = normalize(cross(dg->dpdv, dg->n));
    }
}

#endif


/// Emissive BRDF sampling
float3 Emissive_GetLe(
    // Geometry
    DifferentialGeometry const* dg,
    // Texture args
    TEXTURE_ARG_LIST,
    UberV2ShaderData const* shader_data)
{
    return shader_data->emission_color.xyz;
}

/// BxDF singularity check
bool Bxdf_IsSingular(DifferentialGeometry const* dg)
{
    return (dg->mat.flags & kBxdfFlagsSingular) == kBxdfFlagsSingular;
}

/// BxDF emission check
bool Bxdf_IsEmissive(DifferentialGeometry const* dg)
{
    return (dg->mat.flags & kBxdfFlagsEmissive) == kBxdfFlagsEmissive;
}

bool Bxdf_IsBtdf(DifferentialGeometry const* dg)
{
    return (dg->mat.flags & kBxdfFlagsBrdf) == 0;
}

bool Bxdf_IsReflection(DifferentialGeometry const* dg)
{
    return ((dg->mat.flags & kBxdfFlagsBrdf) == kBxdfFlagsBrdf) && ((dg->mat.flags & kBxdfFlagsDiffuse) == kBxdfFlagsDiffuse);
}

bool Bxdf_IsTransparency(DifferentialGeometry const* dg)
{
    return (dg->mat.flags & kBxdfFlagsTransparency) == kBxdfFlagsTransparency;
}

bool Bxdf_IsRefraction(DifferentialGeometry const* dg)
{
    return Bxdf_IsBtdf(dg) && !Bxdf_IsTransparency(dg);
}
void GetMaterialBxDFType8(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data)
{
	int bxdf_flags = 0;
	const float ndotwi = dot(dg->n, wi);
	bxdf_flags |= kBxdfFlagsBrdf;
	float top_ior = 1.0f;
		if (shader_data->reflection_roughness < ROUGHNESS_EPS)
		{
			bxdf_flags |= kBxdfFlagsSingular;
		}
		Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleReflection);
		Bxdf_SetFlags(dg, bxdf_flags);
		return;
	Bxdf_SetFlags(dg, bxdf_flags);
return;
}

float3 UberV2_Evaluate8(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend(1.0f, shader_data->reflection_ior, UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

float UberV2_GetPdf8(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend_F(1.0f, shader_data->reflection_ior, UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

void UberV2PrepareInputs8(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *data)
{
	int offset = dg->mat.offset + 1;
	data->reflection_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_roughness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy_rotation = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_ior = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_metalness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);

}
float3 UberV2_Sample8(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3* wo, float* pdf,UberV2ShaderData const* shader_data)
{
	const int sampledComponent = Bxdf_UberV2_GetSampledComponent(dg);
	float3 result;
	switch(sampledComponent)
	{
		case kBxdfUberV2SampleReflection: result = UberV2_Reflection_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
	}
	if (false)
	{
		*pdf = UberV2_GetPdf8(dg, wi, *wo, TEXTURE_ARGS, shader_data);
		return UberV2_Evaluate8(dg, wi, *wo, TEXTURE_ARGS, shader_data);
	}
	return result;
}

void GetMaterialBxDFType9(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data)
{
	int bxdf_flags = 0;
	const float ndotwi = dot(dg->n, wi);
	bxdf_flags |= kBxdfFlagsEmissive;
	bxdf_flags |= kBxdfFlagsBrdf;
	float top_ior = 1.0f;
		if (shader_data->reflection_roughness < ROUGHNESS_EPS)
		{
			bxdf_flags |= kBxdfFlagsSingular;
		}
		Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleReflection);
		Bxdf_SetFlags(dg, bxdf_flags);
		return;
	Bxdf_SetFlags(dg, bxdf_flags);
return;
}

float3 UberV2_Evaluate9(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend(1.0f, shader_data->reflection_ior, UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

float UberV2_GetPdf9(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend_F(1.0f, shader_data->reflection_ior, UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

void UberV2PrepareInputs9(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *data)
{
	int offset = dg->mat.offset + 1;
	data->emission_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_roughness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy_rotation = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_ior = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_metalness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);

}
float3 UberV2_Sample9(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3* wo, float* pdf,UberV2ShaderData const* shader_data)
{
	const int sampledComponent = Bxdf_UberV2_GetSampledComponent(dg);
	float3 result;
	switch(sampledComponent)
	{
		case kBxdfUberV2SampleReflection: result = UberV2_Reflection_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
	}
	if (false)
	{
		*pdf = UberV2_GetPdf9(dg, wi, *wo, TEXTURE_ARGS, shader_data);
		return UberV2_Evaluate9(dg, wi, *wo, TEXTURE_ARGS, shader_data);
	}
	return result;
}

void GetMaterialBxDFType12(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data)
{
	int bxdf_flags = 0;
	const float ndotwi = dot(dg->n, wi);
	bxdf_flags |= kBxdfFlagsBrdf;
	float top_ior = 1.0f;
	const float sample3 = Sampler_Sample1D(sampler, SAMPLER_ARGS);
	const float fresnel1 = CalculateFresnel(top_ior, shader_data->coating_ior, ndotwi);
	if (sample3 < fresnel1)
	{
		bxdf_flags |= kBxdfFlagsSingular;
		Bxdf_SetFlags(dg, bxdf_flags);
		Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleCoating);
		return;
	}
	top_ior = shader_data->coating_ior;
		if (shader_data->reflection_roughness < ROUGHNESS_EPS)
		{
			bxdf_flags |= kBxdfFlagsSingular;
		}
		Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleReflection);
		Bxdf_SetFlags(dg, bxdf_flags);
		return;
	Bxdf_SetFlags(dg, bxdf_flags);
return;
}

float3 UberV2_Evaluate12(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend(1.0f, shader_data->coating_ior, UberV2_IdealReflect_Evaluate(shader_data, wi, wo, TEXTURE_ARGS), (Fresnel_Blend(shader_data->coating_ior, shader_data->reflection_ior, UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS)), wi)), wi));
}

float UberV2_GetPdf12(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend_F(1.0f, shader_data->coating_ior, UberV2_IdealReflect_GetPdf(shader_data, wi, wo, TEXTURE_ARGS), (Fresnel_Blend_F(shader_data->coating_ior, shader_data->reflection_ior, UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS)), wi)), wi));
}

void UberV2PrepareInputs12(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *data)
{
	int offset = dg->mat.offset + 1;
	data->coating_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->coating_ior = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_roughness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy_rotation = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_ior = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_metalness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);

}
float3 UberV2_Sample12(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3* wo, float* pdf,UberV2ShaderData const* shader_data)
{
	const int sampledComponent = Bxdf_UberV2_GetSampledComponent(dg);
	float3 result;
	switch(sampledComponent)
	{
		case kBxdfUberV2SampleCoating: result = UberV2_Coating_Sample(shader_data, wi, TEXTURE_ARGS, wo, pdf);
			break;
		case kBxdfUberV2SampleReflection: result = UberV2_Reflection_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
	}
	if (false)
	{
		*pdf = UberV2_GetPdf12(dg, wi, *wo, TEXTURE_ARGS, shader_data);
		return UberV2_Evaluate12(dg, wi, *wo, TEXTURE_ARGS, shader_data);
	}
	return result;
}

void GetMaterialBxDFType16(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data)
{
	int bxdf_flags = 0;
	Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleDiffuse);
	bxdf_flags |= kBxdfFlagsBrdf;
	Bxdf_SetFlags(dg, bxdf_flags);
return;
}

float3 UberV2_Evaluate16(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (UberV2_Lambert_Evaluate(shader_data, wi, wo, TEXTURE_ARGS));
}

float UberV2_GetPdf16(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (UberV2_Lambert_GetPdf(shader_data, wi, wo, TEXTURE_ARGS));
}

void UberV2PrepareInputs16(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *data)
{
	int offset = dg->mat.offset + 1;
	data->diffuse_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);

}
float3 UberV2_Sample16(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3* wo, float* pdf,UberV2ShaderData const* shader_data)
{
	const int sampledComponent = Bxdf_UberV2_GetSampledComponent(dg);
	float3 result;
	switch(sampledComponent)
	{
		case kBxdfUberV2SampleDiffuse: result = UberV2_Lambert_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
	}
	if (false)
	{
		*pdf = UberV2_GetPdf16(dg, wi, *wo, TEXTURE_ARGS, shader_data);
		return UberV2_Evaluate16(dg, wi, *wo, TEXTURE_ARGS, shader_data);
	}
	return result;
}

void GetMaterialBxDFType24(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data)
{
	int bxdf_flags = 0;
	const float ndotwi = dot(dg->n, wi);
	bxdf_flags |= kBxdfFlagsBrdf;
	float top_ior = 1.0f;
	const float fresnel2 = CalculateFresnel(top_ior, shader_data->reflection_ior, ndotwi);
	const float sample4 = Sampler_Sample1D(sampler, SAMPLER_ARGS);
	if (sample4 < fresnel2)
	{
		if (shader_data->reflection_roughness < ROUGHNESS_EPS)
		{
			bxdf_flags |= kBxdfFlagsSingular;
		}
		Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleReflection);
		Bxdf_SetFlags(dg, bxdf_flags);
		return;
	}
	Bxdf_UberV2_SetSampledComponent(dg, kBxdfUberV2SampleDiffuse);
	bxdf_flags |= kBxdfFlagsBrdf;
	Bxdf_SetFlags(dg, bxdf_flags);
return;
}

float3 UberV2_Evaluate24(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend(1.0f, shader_data->reflection_ior, UberV2_Reflection_Evaluate(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Lambert_Evaluate(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

float UberV2_GetPdf24(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	return (Fresnel_Blend_F(1.0f, shader_data->reflection_ior, UberV2_Reflection_GetPdf(shader_data, wi, wo, TEXTURE_ARGS), (UberV2_Lambert_GetPdf(shader_data, wi, wo, TEXTURE_ARGS)), wi));
}

void UberV2PrepareInputs24(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *data)
{
	int offset = dg->mat.offset + 1;
	data->reflection_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_roughness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_anisotropy_rotation = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_ior = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->reflection_metalness = GetInputMapFloat(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);
	data->diffuse_color = GetInputMapFloat4(material_attributes[offset++], dg, input_map_values, TEXTURE_ARGS);

}
float3 UberV2_Sample24(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3* wo, float* pdf,UberV2ShaderData const* shader_data)
{
	const int sampledComponent = Bxdf_UberV2_GetSampledComponent(dg);
	float3 result;
	switch(sampledComponent)
	{
		case kBxdfUberV2SampleReflection: result = UberV2_Reflection_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
		case kBxdfUberV2SampleDiffuse: result = UberV2_Lambert_Sample(shader_data, wi, TEXTURE_ARGS, sample, wo, pdf);
			break;
	}
	if (false)
	{
		*pdf = UberV2_GetPdf24(dg, wi, *wo, TEXTURE_ARGS, shader_data);
		return UberV2_Evaluate24(dg, wi, *wo, TEXTURE_ARGS, shader_data);
	}
	return result;
}

void UberV2PrepareInputs(DifferentialGeometry const* dg, GLOBAL InputMapData const* restrict input_map_values,GLOBAL int const* restrict material_attributes, TEXTURE_ARG_LIST, UberV2ShaderData *shader_data)
{
	switch(dg->mat.layers)
	{
		case 8:
			return UberV2PrepareInputs8(dg, input_map_values, material_attributes, TEXTURE_ARGS, shader_data);
		case 9:
			return UberV2PrepareInputs9(dg, input_map_values, material_attributes, TEXTURE_ARGS, shader_data);
		case 12:
			return UberV2PrepareInputs12(dg, input_map_values, material_attributes, TEXTURE_ARGS, shader_data);
		case 16:
			return UberV2PrepareInputs16(dg, input_map_values, material_attributes, TEXTURE_ARGS, shader_data);
		case 24:
			return UberV2PrepareInputs24(dg, input_map_values, material_attributes, TEXTURE_ARGS, shader_data);
	}
}
void GetMaterialBxDFType(float3 wi, Sampler* sampler, SAMPLER_ARG_LIST, DifferentialGeometry* dg, UberV2ShaderData const* shader_data){
	switch(dg->mat.layers)
	{
		case 8:
			return GetMaterialBxDFType8(wi, sampler, SAMPLER_ARGS, dg, shader_data);
			break;
		case 9:
			return GetMaterialBxDFType9(wi, sampler, SAMPLER_ARGS, dg, shader_data);
			break;
		case 12:
			return GetMaterialBxDFType12(wi, sampler, SAMPLER_ARGS, dg, shader_data);
			break;
		case 16:
			return GetMaterialBxDFType16(wi, sampler, SAMPLER_ARGS, dg, shader_data);
			break;
		case 24:
			return GetMaterialBxDFType24(wi, sampler, SAMPLER_ARGS, dg, shader_data);
			break;
	}
}
float UberV2_GetPdf(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	float3 wi_t = matrix_mul_vector3(dg->world_to_tangent, wi);
	float3 wo_t = matrix_mul_vector3(dg->world_to_tangent, wo);
	switch (dg->mat.layers)
	{
		case 8:
			return UberV2_GetPdf8(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 9:
			return UberV2_GetPdf9(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 12:
			return UberV2_GetPdf12(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 16:
			return UberV2_GetPdf16(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 24:
			return UberV2_GetPdf24(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
	}
	return 0.0f;
}
float3 UberV2_Evaluate(DifferentialGeometry const* dg, float3 wi, float3 wo, TEXTURE_ARG_LIST, UberV2ShaderData const* shader_data)
{
	float3 wi_t = matrix_mul_vector3(dg->world_to_tangent, wi);
	float3 wo_t = matrix_mul_vector3(dg->world_to_tangent, wo);
	switch (dg->mat.layers)
	{
		case 8:
			return UberV2_Evaluate8(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 9:
			return UberV2_Evaluate9(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 12:
			return UberV2_Evaluate12(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 16:
			return UberV2_Evaluate16(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
		case 24:
			return UberV2_Evaluate24(dg, wi_t, wo_t, TEXTURE_ARGS, shader_data);
	}
	return (float3)(0.0f);
}
float3 UberV2_Sample(DifferentialGeometry const* dg, float3 wi, TEXTURE_ARG_LIST, float2 sample, float3 *wo, float *pdf, UberV2ShaderData const* shader_data)
{
	float3 wi_t = matrix_mul_vector3(dg->world_to_tangent, wi);
	float3 wo_t;
	float3 res = 0.f;
	switch (dg->mat.layers)
	{
		case 8:
			res = UberV2_Sample8(dg, wi_t, TEXTURE_ARGS, sample, &wo_t, pdf, shader_data);
			break;
		case 9:
			res = UberV2_Sample9(dg, wi_t, TEXTURE_ARGS, sample, &wo_t, pdf, shader_data);
			break;
		case 12:
			res = UberV2_Sample12(dg, wi_t, TEXTURE_ARGS, sample, &wo_t, pdf, shader_data);
			break;
		case 16:
			res = UberV2_Sample16(dg, wi_t, TEXTURE_ARGS, sample, &wo_t, pdf, shader_data);
			break;
		case 24:
			res = UberV2_Sample24(dg, wi_t, TEXTURE_ARGS, sample, &wo_t, pdf, shader_data);
			break;
	}
	*wo = matrix_mul_vector3(dg->tangent_to_world, wo_t);	return res;
}


#endif // BXDF_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef LIGHT_CL
#define LIGHT_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef SCENE_CL
#define SCENE_CL


typedef struct
{
    // Vertices
    GLOBAL float3 const* restrict vertices;
    // Normals
    GLOBAL float3 const* restrict normals;
    // UVs
    GLOBAL float2 const* restrict uvs;
    // Indices
    GLOBAL int const* restrict indices;
    // Shapes
    GLOBAL Shape const* restrict shapes;
    // Material attributes
    GLOBAL int const* restrict material_attributes;
    // Input map values
    GLOBAL InputMapData const* restrict input_map_values;
    // Emissive objects
    GLOBAL Light const* restrict lights;
    // Envmap idx
    int env_light_idx;
    // Number of emissive objects
    int num_lights;
    // Light distribution 
    GLOBAL int const* restrict light_distribution;
} Scene;

// Get triangle vertices given scene, shape index and prim index
INLINE void Scene_GetTriangleVertices(Scene const* scene, int shape_idx, int prim_idx, float3* v0, float3* v1, float3* v2)
{
    // Extract shape data
    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch positions and transform to world space
    *v0 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i0]);
    *v1 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i1]);
    *v2 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i2]);
}

// Get triangle uvs given scene, shape index and prim index
INLINE void Scene_GetTriangleUVs(Scene const* scene, int shape_idx, int prim_idx, float2* uv0, float2* uv1, float2* uv2)
{
    // Extract shape data
    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch positions and transform to world space
    *uv0 = scene->uvs[shape.startvtx + i0];
    *uv1 = scene->uvs[shape.startvtx + i1];
    *uv2 = scene->uvs[shape.startvtx + i2];
}


// Interpolate position, normal and uv
INLINE void Scene_InterpolateAttributes(Scene const* scene, int shape_idx, int prim_idx, float2 barycentrics, float3* p, float3* n, float2* uv, float* area)
{
    // Extract shape data
    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch normals
    float3 n0 = scene->normals[shape.startvtx + i0];
    float3 n1 = scene->normals[shape.startvtx + i1];
    float3 n2 = scene->normals[shape.startvtx + i2];

    // Fetch positions and transform to world space
    float3 v0 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i0]);
    float3 v1 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i1]);
    float3 v2 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i2]);

    // Fetch UVs
    float2 uv0 = scene->uvs[shape.startvtx + i0];
    float2 uv1 = scene->uvs[shape.startvtx + i1];
    float2 uv2 = scene->uvs[shape.startvtx + i2];

    // Calculate barycentric position and normal
    *p = (1.f - barycentrics.x - barycentrics.y) * v0 + barycentrics.x * v1 + barycentrics.y * v2;
    *n = normalize(matrix_mul_vector3(shape.transform, (1.f - barycentrics.x - barycentrics.y) * n0 + barycentrics.x * n1 + barycentrics.y * n2));
    *uv = (1.f - barycentrics.x - barycentrics.y) * uv0 + barycentrics.x * uv1 + barycentrics.y * uv2;
    *area = 0.5f * length(cross(v2 - v0, v1 - v0));
}

// Interpolate position, normal and uv
INLINE void Scene_InterpolateVertices(Scene const* scene, int shape_idx, int prim_idx, float2 barycentrics, float3* p)
{
    // Extract shape data
    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch positions and transform to world space
    float3 v0 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i0]);
    float3 v1 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i1]);
    float3 v2 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i2]);

    // Calculate barycentric position and normal
    *p = (1.f - barycentrics.x - barycentrics.y) * v0 + barycentrics.x * v1 + barycentrics.y * v2;
}

// Interpolate position, normal and uv
INLINE void Scene_InterpolateVerticesFromIntersection(Scene const* scene, Intersection const* isect, float3* p)
{
    // Extract shape data
    int shape_idx = isect->shapeid - 1;
    int prim_idx = isect->primid;
    float2 barycentrics = isect->uvwt.xy;

    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch positions and transform to world space
    float3 v0 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i0]);
    float3 v1 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i1]);
    float3 v2 = matrix_mul_point3(shape.transform, scene->vertices[shape.startvtx + i2]);

    // Calculate barycentric position and normal
    *p = (1.f - barycentrics.x - barycentrics.y) * v0 + barycentrics.x * v1 + barycentrics.y * v2;
}

// Interpolate position, normal and uv
INLINE void Scene_InterpolateNormalsFromIntersection(Scene const* scene, Intersection const* isect, float3* n)
{
    // Extract shape data
    int shape_idx = isect->shapeid - 1;
    int prim_idx = isect->primid;
    float2 barycentrics = isect->uvwt.xy;

    Shape shape = scene->shapes[shape_idx];

    // Fetch indices starting from startidx and offset by prim_idx
    int i0 = scene->indices[shape.startidx + 3 * prim_idx];
    int i1 = scene->indices[shape.startidx + 3 * prim_idx + 1];
    int i2 = scene->indices[shape.startidx + 3 * prim_idx + 2];

    // Fetch normals
    float3 n0 = scene->normals[shape.startvtx + i0];
    float3 n1 = scene->normals[shape.startvtx + i1];
    float3 n2 = scene->normals[shape.startvtx + i2];

    // Calculate barycentric position and normal
    *n = normalize(matrix_mul_vector3(shape.transform, (1.f - barycentrics.x - barycentrics.y) * n0 + barycentrics.x * n1 + barycentrics.y * n2));
}

INLINE int Scene_GetVolumeIndex(Scene const* scene, int shape_idx)
{
    Shape shape = scene->shapes[shape_idx];
    return shape.volume_idx;
}

/// Fill DifferentialGeometry structure based on intersection info from RadeonRays
void Scene_FillDifferentialGeometry(// Scene
                              Scene const* scene,
                              // RadeonRays intersection
                              Intersection const* isect,
                              // Differential geometry
                              DifferentialGeometry* diffgeo
                              )
{
    // Determine shape and polygon
    int shape_idx = isect->shapeid - 1;
    int prim_idx = isect->primid;

    // Get barycentrics
    float2 barycentrics = isect->uvwt.xy;

    // Extract shape data
    Shape shape = scene->shapes[shape_idx];

    // Interpolate attributes
    float3 p;
    float3 n;
    float2 uv;
    float area;
    Scene_InterpolateAttributes(scene, shape_idx, prim_idx, barycentrics, &p, &n, &uv, &area);
    // Triangle area (for area lighting)
    diffgeo->area = area;

    // Calculate barycentric position and normal
    diffgeo->n = n;
    diffgeo->p = p;
    diffgeo->uv = uv;

    // Get vertices
    float3 v0, v1, v2;
    Scene_GetTriangleVertices(scene, shape_idx, prim_idx, &v0, &v1, &v2);

    // Calculate true normal
    diffgeo->ng = normalize(cross(v1 - v0, v2 - v0));

    // Get material at shading point
    diffgeo->mat = shape.material;

    // Get UVs
    float2 uv0, uv1, uv2;
    Scene_GetTriangleUVs(scene, shape_idx, prim_idx, &uv0, &uv1, &uv2);

    // Reverse geometric normal if shading normal points to different side
    if (dot(diffgeo->ng, diffgeo->n) < 0.f)
    {
        diffgeo->ng = -diffgeo->ng;
    }

    /// Calculate tangent basis
    /// From PBRT book
    float du1 = uv0.x - uv2.x;
    float du2 = uv1.x - uv2.x;
    float dv1 = uv0.y - uv2.y;
    float dv2 = uv1.y - uv2.y;
    float3 dp1 = v0 - v2;
    float3 dp2 = v1 - v2;
    float det = du1 * dv2 - dv1 * du2;

    if (0 && det != 0.f)
    {
        float invdet = 1.f / det;
        diffgeo->dpdu = normalize( (dv2 * dp1 - dv1 * dp2) * invdet );
        diffgeo->dpdv = normalize( (-du2 * dp1 + du1 * dp2) * invdet );
        diffgeo->dpdu -= dot(diffgeo->n, diffgeo->dpdu) * diffgeo->n;
        diffgeo->dpdu = normalize(diffgeo->dpdu);
        
        diffgeo->dpdv -= dot(diffgeo->n, diffgeo->dpdv) * diffgeo->n;
        diffgeo->dpdv -= dot(diffgeo->dpdu, diffgeo->dpdv) * diffgeo->dpdu;
        diffgeo->dpdv = normalize(diffgeo->dpdv);
    }
    else
    {
        diffgeo->dpdu = normalize(GetOrthoVector(diffgeo->n));
        diffgeo->dpdv = normalize(cross(diffgeo->n, diffgeo->dpdu));
    }
}


// Calculate tangent transform matrices inside differential geometry
INLINE void DifferentialGeometry_CalculateTangentTransforms(DifferentialGeometry* diffgeo)
{
    diffgeo->world_to_tangent = matrix_from_rows3(diffgeo->dpdu, diffgeo->n, diffgeo->dpdv);

    diffgeo->world_to_tangent.m0.w = -dot(diffgeo->dpdu, diffgeo->p);
    diffgeo->world_to_tangent.m1.w = -dot(diffgeo->n, diffgeo->p);
    diffgeo->world_to_tangent.m2.w = -dot(diffgeo->dpdv, diffgeo->p);

    diffgeo->tangent_to_world = matrix_from_cols3(diffgeo->world_to_tangent.m0.xyz, 
        diffgeo->world_to_tangent.m1.xyz, diffgeo->world_to_tangent.m2.xyz);

    diffgeo->tangent_to_world.m0.w = diffgeo->p.x;
    diffgeo->tangent_to_world.m1.w = diffgeo->p.y;
    diffgeo->tangent_to_world.m2.w = diffgeo->p.z;
}

#define POWER_SAMPLING

// Sample light index
INLINE int Scene_SampleLight(Scene const* scene, float sample, float* pdf)
{
#ifndef POWER_SAMPLING
    int num_lights = scene->num_lights;
    int light_idx = clamp((int)(sample * num_lights), 0, num_lights - 1);
    *pdf = 1.f / num_lights;
    return light_idx;
#else
    int num_lights = scene->num_lights;
    int light_idx = Distribution1D_SampleDiscrete(sample, scene->light_distribution, pdf);
    return light_idx;
#endif
}

#endif
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef PATH_CL
#define PATH_CL


typedef struct _Path
{
    float3 throughput;
    int volume;
    int flags;
    int active;
    int extra1;
} Path;

typedef enum _PathFlags
{
    kNone = 0x0,
    kKilled = 0x1,
    kScattered = 0x2,
    kOpaque = 0x4
} PathFlags;

INLINE bool Path_IsScattered(__global Path const* path)
{
    return path->flags & kScattered;
}

INLINE bool Path_IsAlive(__global Path const* path)
{
    return ((path->flags & kKilled) == 0);
}

INLINE bool Path_ContainsOpacity(__global Path const* path)
{
    return path->flags & kOpaque;
}

INLINE void Path_ClearScatterFlag(__global Path* path)
{
    path->flags &= ~kScattered;
}

INLINE void Path_SetScatterFlag(__global Path* path)
{
    path->flags |= kScattered;
}

INLINE void Path_SetOpacityFlag(__global Path* path)
{
    path->flags |= kOpaque;
}

INLINE void Path_ClearBxdfFlags(__global Path* path)
{
    path->flags &= (kKilled | kScattered | kOpaque);
}

INLINE int Path_GetBxdfFlags(__global Path const* path)
{
    return path->flags >> 8;
}

INLINE int Path_SetBxdfFlags(__global Path* path, int flags)
{
    return path->flags |= (flags << 8);
}

INLINE void Path_Restart(__global Path* path)
{
    path->flags = 0;
}

INLINE int Path_GetVolumeIdx(__global Path const* path)
{
    return path->volume;
}

INLINE void Path_SetVolumeIdx(__global Path* path, int volume_idx)
{
    path->volume = volume_idx;
}

INLINE float3 Path_GetThroughput(__global Path const* path)
{
    float3 t = path->throughput;
    return t;
}

INLINE void Path_MulThroughput(__global Path* path, float3 mul)
{
    path->throughput *= mul;
}

INLINE void Path_Kill(__global Path* path)
{
    path->flags |= kKilled;
}

INLINE void Path_AddContribution(__global Path* path, __global float3* output, int idx, float3 val)
{
    output[idx] += Path_GetThroughput(path) * val;
}

INLINE bool Path_IsSpecular(__global Path const* path)
{
    int flags = Path_GetBxdfFlags(path);
    return (flags & kBxdfFlagsSingular) == kBxdfFlagsSingular;
}

INLINE void Path_SetFlags(DifferentialGeometry* diffgeo, GLOBAL Path* restrict path)
{
    Path_ClearBxdfFlags(path);
    Path_SetBxdfFlags(path, Bxdf_GetFlags(diffgeo));
}

#endif


enum LightInteractionType
{
    kLightInteractionUnknown,
    kLightInteractionSurface,
    kLightInteractionVolume
};

INLINE
bool IntersectTriangle(ray const* r, float3 v1, float3 v2, float3 v3, float* a, float* b)
{
    const float3 e1 = v2 - v1;
    const float3 e2 = v3 - v1;
    const float3 s1 = cross(r->d.xyz, e2);
    const float  invd = native_recip(dot(s1, e1));
    const float3 d = r->o.xyz - v1;
    const float  b1 = dot(d, s1) * invd;
    const float3 s2 = cross(d, e1);
    const float  b2 = dot(r->d.xyz, s2) * invd;
    const float temp = dot(e2, s2) * invd;

    if (b1 < 0.f || b1 > 1.f || b2 < 0.f || b1 + b2 > 1.f)
    {
        return false;
    }
    else
    {
        *a = b1;
        *b = b2;
        return true;
    }
}

INLINE int EnvironmentLight_GetTexture(Light const* light, int bxdf_flags)
{
    int tex = light->tex;

    if ((bxdf_flags & kBxdfFlagsBrdf) && (light->tex_reflection != -1) && ((bxdf_flags & kBxdfFlagsDiffuse) != kBxdfFlagsDiffuse))
        tex = light->tex_reflection;

    if (((bxdf_flags & kBxdfFlagsBrdf) == 0) && light->tex_refraction != -1)
        tex = light->tex_refraction;

    if ((bxdf_flags & kBxdfFlagsTransparency) && light->tex_transparency != -1)
        tex = light->tex_transparency;

    return tex;
}

INLINE int EnvironmentLight_GetBackgroundTexture(Light const* light)
{
    return light->tex_background == -1 ? light->tex : light->tex_background;
}

/*
 Environment light
 */
/// Get intensity for a given direction
float3 EnvironmentLight_GetLe(// Light
                              Light const* light,
                              // Scene
                              Scene const* scene,
                              // Geometry
                              DifferentialGeometry const* dg,
                              // Path flags
                              int bxdf_flags,
                              // Light inteaction type
                              int interaction_type,
                              // Direction to light source
                              float3* wo,
                              // Textures
                              TEXTURE_ARG_LIST
                              )
{
    // Sample envmap
    *wo *= CRAZY_HIGH_DISTANCE;

    int tex = EnvironmentLight_GetTexture(light, bxdf_flags);

    if (tex == -1)
    {
        return 0.f;
    }

    return light->multiplier * Texture_SampleEnvMap(normalize(*wo), TEXTURE_ARGS_IDX(tex), light->ibl_mirror_x);
}

/// Sample direction to the light
float3 EnvironmentLight_Sample(// Light
                               Light const* light,
                               // Scene
                               Scene const* scene,
                               // Geometry
                               DifferentialGeometry const* dg,
                               // Textures
                               TEXTURE_ARG_LIST,
                               // Sample
                               float2 sample,
                               // Path flags
                               int bxdf_flags,
                               // Light inteaction type
                               int interaction_type,
                               // Direction to light source
                               float3* wo,
                               // PDF
                               float* pdf
                              )
{
    float3 d;

    if (interaction_type != kLightInteractionVolume)
    {
        d = Sample_MapToHemisphere(sample, dg->n, 0.f);
        *pdf = 1.f / (2.f * PI);
    }
    else
    {
        d = Sample_MapToSphere(sample);
        *pdf = 1.f / (4.f * PI);
    }

    // Generate direction
    *wo = CRAZY_HIGH_DISTANCE * d;

    int tex = EnvironmentLight_GetTexture(light, bxdf_flags);

    if (tex == -1)
    {
        *pdf = 0.f;
        return 0.f;
    }

    // Sample envmap
    return light->multiplier * Texture_SampleEnvMap(d, TEXTURE_ARGS_IDX(tex), light->ibl_mirror_x);
}

/// Get PDF for a given direction
float EnvironmentLight_GetPdf(
                              // Light
                              Light const* light,
                              // Scene
                              Scene const* scene,
                              // Geometry
                              DifferentialGeometry const* dg,
                              // Path flags
                              int bxdf_flags,
                              // Light inteaction type
                              int interaction_type,
                              // Direction to light source
                              float3 wo,
                              // Textures
                              TEXTURE_ARG_LIST
                              )
{
    if (interaction_type != kLightInteractionVolume)
    {
        return 1.f / (2.f * PI);
    }
    else
    {
        return 1.f / (4.f * PI);
    }
}


/*
 Area light
 */
// Get intensity for a given direction
float3 AreaLight_GetLe(// Emissive object
                       Light const* light,
                       // Scene
                       Scene const* scene,
                       // Geometry
                       DifferentialGeometry const* dg,
                       // Direction to light source
                       float3* wo,
                       // Textures
                       TEXTURE_ARG_LIST
                       )
{
    ray r;
    r.o.xyz = dg->p;
    r.d.xyz = *wo;

    int shapeidx = light->shapeidx;
    int primidx = light->primidx;

    float3 v0, v1, v2;
    Scene_GetTriangleVertices(scene, shapeidx, primidx, &v0, &v1, &v2);

    float a, b;
    if (IntersectTriangle(&r, v0, v1, v2, &a, &b))
    {
        float3 n;
        float3 p;
        float2 tx;
        float area;
        Scene_InterpolateAttributes(scene, shapeidx, primidx, make_float2(a, b), &p, &n, &tx, &area);

        float3 d = p - dg->p;
        *wo = d;

        int material_offset = scene->shapes[shapeidx].material.offset;

        const float3 ke = GetUberV2EmissionColor(material_offset, dg, scene->input_map_values, scene->material_attributes, TEXTURE_ARGS).xyz;
        
        return ke;
    }
    else
    {
        return make_float3(0.f, 0.f, 0.f);
    }
}

/// Sample direction to the light
float3 AreaLight_Sample(// Emissive object
                        Light const* light,
                        // Scene
                        Scene const* scene,
                        // Geometry
                        DifferentialGeometry const* dg,
                        // Textures
                        TEXTURE_ARG_LIST,
                        // Sample
                        float2 sample,
                        // Direction to light source
                        float3* wo,
                        // PDF
                        float* pdf)
{
    int shapeidx = light->shapeidx;
    int primidx = light->primidx;

    // Generate sample on triangle
    float r0 = sample.x;
    float r1 = sample.y;

    // Convert random to barycentric coords
    float2 uv;

    uv.x = 1.f - native_sqrt(r0);
    uv.y = native_sqrt(r0) * r1;

    float3 n;
    float3 p;
    float2 tx;
    float area;
    Scene_InterpolateAttributes(scene, shapeidx, primidx, uv, &p, &n, &tx, &area);

    *wo = p - dg->p;

    int material_offset = scene->shapes[shapeidx].material.offset;

    const float3 ke = GetUberV2EmissionColor(material_offset, dg, scene->input_map_values, scene->material_attributes, TEXTURE_ARGS).xyz;
    float3 v = -normalize(*wo);

    float ndotv = dot(n, v);

    if (ndotv > 0.f)
    {
        float dist2 = dot(*wo, *wo);
        float denom = fabs(ndotv) * area;
        *pdf = denom > 0.f ? dist2 / denom : 0.f;
        return ke;
    }
    else
    {
        *pdf = 0.f;
        return 0.f;
    }
}

/// Get PDF for a given direction
float AreaLight_GetPdf(// Emissive object
                       Light const* light,
                       // Scene
                       Scene const* scene,
                       // Geometry
                       DifferentialGeometry const* dg,
                       // Direction to light source
                       float3 wo,
                       // Textures
                       TEXTURE_ARG_LIST
                       )
{
    ray r;
    r.o.xyz = dg->p;
    r.d.xyz = wo;

    int shapeidx = light->shapeidx;
    int primidx = light->primidx;

    float3 v0, v1, v2;
    Scene_GetTriangleVertices(scene, shapeidx, primidx, &v0, &v1, &v2);

    // Intersect ray against this area light
    float a, b;
    if (IntersectTriangle(&r, v0, v1, v2, &a, &b))
    {
        float3 n;
        float3 p;
        float2 tx;
        float area;
        Scene_InterpolateAttributes(scene, shapeidx, primidx, make_float2(a, b), &p, &n, &tx, &area);

        float3 d = p - dg->p;
        float dist2 = dot(d, d) ;
        float denom = (fabs(dot(-normalize(d), n)) * area);

        return denom > 0.f ? dist2 / denom : 0.f;
    }
    else
    {
        return 0.f;
    }
}

float3 AreaLight_SampleVertex(
    // Emissive object
    Light const* light,
    // Scene
    Scene const* scene,
    // Textures
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample0,
    float2 sample1,
    // Direction to light source
    float3* p,
    float3* n,
    float3* wo,
    // PDF
    float* pdf)
{
    int shapeidx = light->shapeidx;
    int primidx = light->primidx;

    // Generate sample on triangle
    float r0 = sample0.x;
    float r1 = sample0.y;

    // Convert random to barycentric coords
    float2 uv;
    uv.x = native_sqrt(r0) * (1.f - r1);
    uv.y = native_sqrt(r0) * r1;

    float2 tx;
    float area;
    Scene_InterpolateAttributes(scene, shapeidx, primidx, uv, p, n, &tx, &area);

    int material_offset = scene->shapes[shapeidx].material.offset;

    /*const float3 ke = GetUberV2EmissionColor(material_offset, dg, scene->input_map_values, scene->material_attributes, TEXTURE_ARGS).xyz;*/
    const float3 ke = make_float3(0.f, 0.f, 0.f);
    *wo = Sample_MapToHemisphere(sample1, *n, 1.f);
    *pdf = (1.f / area) * fabs(dot(*n, *wo)) / PI;

    return ke;
}

/*
Directional light
*/
// Get intensity for a given direction
float3 DirectionalLight_GetLe(// Emissive object
    Light const* light,
    // Scene
    Scene const* scene,
    // Geometry
    DifferentialGeometry const* dg,
    // Direction to light source
    float3* wo,
    // Textures
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

/// Sample direction to the light
float3 DirectionalLight_Sample(// Emissive object
    Light const* light,
    // Scene
    Scene const* scene,
    // Geometry
    DifferentialGeometry const* dg,
    // Textures
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample,
    // Direction to light source
    float3* wo,
    // PDF
    float* pdf)
{
    *wo = CRAZY_HIGH_DISTANCE * -light->d;
    *pdf = 1.f;
    return light->intensity;
}

/// Get PDF for a given direction
float DirectionalLight_GetPdf(// Emissive object
    Light const* light,
    // Scene
    Scene const* scene,
    // Geometry
    DifferentialGeometry const* dg,
    // Direction to light source
    float3 wo,
    // Textures
    TEXTURE_ARG_LIST
)
{
    return 0.f;
}

/*
 Point light
 */
// Get intensity for a given direction
float3 PointLight_GetLe(// Emissive object
                              Light const* light,
                              // Scene
                              Scene const* scene,
                              // Geometry
                              DifferentialGeometry const* dg,
                              // Direction to light source
                              float3* wo,
                              // Textures
                              TEXTURE_ARG_LIST
                              )
{
    return 0.f;
}

/// Sample direction to the light
float3 PointLight_Sample(// Emissive object
                               Light const* light,
                               // Scene
                               Scene const* scene,
                               // Geometry
                               DifferentialGeometry const* dg,
                               // Textures
                               TEXTURE_ARG_LIST,
                               // Sample
                               float2 sample,
                               // Direction to light source
                               float3* wo,
                               // PDF
                               float* pdf)
{
    *wo = light->p - dg->p;
    *pdf = 1.f;
    return light->intensity / dot(*wo, *wo);
}

/// Get PDF for a given direction
float PointLight_GetPdf(// Emissive object
                              Light const* light,
                              // Scene
                              Scene const* scene,
                              // Geometry
                              DifferentialGeometry const* dg,
                              // Direction to light source
                              float3 wo,
                              // Textures
                              TEXTURE_ARG_LIST
                              )
{
    return 0.f;
}

/// Sample vertex on the light
float3 PointLight_SampleVertex(
    // Light object
    Light const* light,
    // Scene
    Scene const* scene,
    // Textures
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample0,
    float2 sample1,
    // Direction to light source
    float3* p,
    float3* n,
    float3* wo,
    // PDF
    float* pdf)
{
    *p = light->p;
    *n = make_float3(0.f, 1.f, 0.f);
    *wo = Sample_MapToSphere(sample0);
    *pdf = 1.f / (4.f * PI);
    return light->intensity;
}

/*
 Spot light
 */
// Get intensity for a given direction
float3 SpotLight_GetLe(// Emissive object
                        Light const* light,
                        // Scene
                        Scene const* scene,
                        // Geometry
                        DifferentialGeometry const* dg,
                        // Direction to light source
                        float3* wo,
                        // Textures
                        TEXTURE_ARG_LIST
                        )
{
    return 0.f;
}

/// Sample direction to the light
float3 SpotLight_Sample(// Emissive object
                         Light const* light,
                         // Scene
                         Scene const* scene,
                         // Geometry
                         DifferentialGeometry const* dg,
                         // Textures
                         TEXTURE_ARG_LIST,
                         // Sample
                         float2 sample,
                         // Direction to light source
                         float3* wo,
                         // PDF
                         float* pdf)
{
    *wo = light->p - dg->p;
    float ddotwo = dot(-normalize(*wo), light->d);
    
    if (ddotwo > light->oa)
    {
        float3 intensity = light->intensity / dot(*wo, *wo);
        *pdf = 1.f;
        return ddotwo > light->ia ? intensity : intensity * (1.f - (light->ia - ddotwo) / (light->ia - light->oa));
    }
    else
    {
        *pdf = 0.f;
        return 0.f;
    }
}

/// Get PDF for a given direction
float SpotLight_GetPdf(// Emissive object
                        Light const* light,
                        // Scene
                        Scene const* scene,
                        // Geometry
                        DifferentialGeometry const* dg,
                        // Direction to light source
                        float3 wo,
                        // Textures
                        TEXTURE_ARG_LIST
                        )
{
    return 0.f;
}


/*
 Dispatch calls
 */

/// Get intensity for a given direction
float3 Light_GetLe(// Light index
                   int idx,
                   // Scene
                   Scene const* scene,
                   // Geometry
                   DifferentialGeometry const* dg,
                   // Path flags
                    int bxdf_flags,
                   // Light inteaction type
                   int interaction_type,
                   // Direction to light source
                   float3* wo,
                   // Textures
                   TEXTURE_ARG_LIST
                   )
{
    Light light = scene->lights[idx];

    switch(light.type)
    {
        case kIbl:
            return EnvironmentLight_GetLe(&light, scene, dg, bxdf_flags, interaction_type, wo, TEXTURE_ARGS);
        case kArea:
            return AreaLight_GetLe(&light, scene, dg, wo, TEXTURE_ARGS);
        case kDirectional:
            return DirectionalLight_GetLe(&light, scene, dg, wo, TEXTURE_ARGS);
        case kPoint:
            return PointLight_GetLe(&light, scene, dg, wo, TEXTURE_ARGS);
        case kSpot:
            return SpotLight_GetLe(&light, scene, dg, wo, TEXTURE_ARGS);
    }

    return make_float3(0.f, 0.f, 0.f);
}

/// Sample direction to the light
float3 Light_Sample(// Light index
                    int idx,
                    // Scene
                    Scene const* scene,
                    // Geometry
                    DifferentialGeometry const* dg,
                    // Textures
                    TEXTURE_ARG_LIST,
                    // Sample
                    float2 sample,
                    // Path flags
                    int bxdf_flags,
                    // Light inteaction type
                    int interaction_type,
                    // Direction to light source
                    float3* wo,
                    // PDF
                    float* pdf)
{
    Light light = scene->lights[idx];

    switch(light.type)
    {
        case kIbl:
            return EnvironmentLight_Sample(&light, scene, dg, TEXTURE_ARGS, sample, bxdf_flags, interaction_type, wo, pdf);
        case kArea:
            return AreaLight_Sample(&light, scene, dg, TEXTURE_ARGS, sample, wo, pdf);
        case kDirectional:
            return DirectionalLight_Sample(&light, scene, dg, TEXTURE_ARGS, sample, wo, pdf);
        case kPoint:
            return PointLight_Sample(&light, scene, dg, TEXTURE_ARGS, sample, wo, pdf);
        case kSpot:
            return SpotLight_Sample(&light, scene, dg, TEXTURE_ARGS, sample, wo, pdf);
    }

    *pdf = 0.f;
    return make_float3(0.f, 0.f, 0.f);
}

/// Get PDF for a given direction
float Light_GetPdf(// Light index
                   int idx,
                   // Scene
                   Scene const* scene,
                   // Geometry
                   DifferentialGeometry const* dg,
                    // Path flags
                    int bxdf_flags,
                    // Light inteaction type
                    int interaction_type,
                   // Direction to light source
                   float3 wo,
                   // Textures
                   TEXTURE_ARG_LIST
                   )
{
    Light light = scene->lights[idx];

    switch(light.type)
    {
        case kIbl:
            return EnvironmentLight_GetPdf(&light, scene, dg, bxdf_flags, interaction_type, wo, TEXTURE_ARGS);
        case kArea:
            return AreaLight_GetPdf(&light, scene, dg, wo, TEXTURE_ARGS);
        case kDirectional:
            return DirectionalLight_GetPdf(&light, scene, dg, wo, TEXTURE_ARGS);
        case kPoint:
            return PointLight_GetPdf(&light, scene, dg, wo, TEXTURE_ARGS);
        case kSpot:
            return SpotLight_GetPdf(&light, scene, dg, wo, TEXTURE_ARGS);
    }

    return 0.f;
}

/// Sample vertex on the light
float3 Light_SampleVertex(
    // Light index
    int idx,
    // Scene
    Scene const* scene,
    // Textures
    TEXTURE_ARG_LIST,
    // Sample
    float2 sample0,
    float2 sample1,
    // Point on the light
    float3* p,
    // Normal at light vertex
    float3* n,
    // Direction
    float3* wo,
    // PDF
    float* pdf)
{
    Light light = scene->lights[idx];

    switch (light.type)
    {
        case kArea:
            return AreaLight_SampleVertex(&light, scene, TEXTURE_ARGS, sample0, sample1, p, n, wo, pdf);
        case kPoint:
            return PointLight_SampleVertex(&light, scene, TEXTURE_ARGS, sample0, sample1, p, n, wo, pdf);
    }

    *pdf = 0.f;
    return make_float3(0.f, 0.f, 0.f);
}

/// Check if the light is singular
bool Light_IsSingular(__global Light const* light)
{
    return light->type == kPoint ||
        light->type == kSpot ||
        light->type == kDirectional;
}

#endif // LIGHT_CLnv
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef VOLUMETRICS_CL
#define VOLUMETRICS_CL


#define FAKE_SHAPE_SENTINEL 0xFFFFFF

float PhaseFunctionHG(float3 wi, float3 wo, float g)
{
    float costheta = dot(wi, wo);
    return 1.f / (4.f * PI) *
        (1.f - g*g) / native_powr(1.f + g*g - 2.f * g * costheta, 1.5f);
}

// See PBRT for derivation
float PhaseFunctionHG_Sample(float3 wi, float g, float2 sample, float3* wo)
{
    float costheta = 0.f;
    if (fabs(g) < 1e-5)
    {
        costheta = 1.f - 2.f * sample.x;
    }
    else
    {
        float temp = (1.f - g * g) / (1.f - g + 2.f * g * sample.x);
        costheta = (1 + g * g - temp * temp) / (2.f * g);
    }

    float phi = 2.f * PI * sample.y;

    float3 u = GetOrthoVector(-wi);
    float3 v = normalize(cross(-wi, u));
    *wo = u * native_cos(phi) + v * native_sin(phi) - wi * costheta;

    return PhaseFunctionHG(wi, *wo, g);
}

// Evaluate volume transmittance along the ray [0, dist] segment
float3 Volume_Transmittance(GLOBAL Volume const* volume, GLOBAL ray const* ray, float dist)
{
    switch (volume->type)
    {
        case kHomogeneous:
        {
            // For homogeneous it is e(-sigma * dist)
            float3 sigma_t = TEXTURED_INPUT_GET_COLOR(volume->sigma_a) +
                             TEXTURED_INPUT_GET_COLOR(volume->sigma_s);
            return native_exp(-sigma_t * dist);
        }
    }
    
    return 1.f;
}

// Evaluate volume selfemission along the ray [0, dist] segment
float3 Volume_Emission(GLOBAL Volume const* volume, GLOBAL ray const* ray, float dist)
{
    switch (volume->type)
    {
        case kHomogeneous:
        {
            return TEXTURED_INPUT_GET_COLOR(volume->sigma_e) * dist;
        }
    }
    
    return 0.f;
}

// Sample volume in order to find next scattering event
float Volume_SampleDistance(GLOBAL Volume const* volume, GLOBAL ray const* ray, float maxdist, float2 sample, float* pdf)
{
    // Sample component
    float3 sigma_s = TEXTURED_INPUT_GET_COLOR(volume->sigma_s);
    float sigma = sample.x < 0.33f ? sigma_s.x :
                  sample.x < 0.66f ? sigma_s.y : sigma_s.z;

    switch (volume->type)
    {
        case kHomogeneous:
        {
            
            float d = sigma > 0.f ? (-native_log(sample.y) / sigma) : -1.f;
            float temp = (1.f / 3.f) * (sigma_s.x * native_exp(-sigma_s.x * d)
                + sigma_s.y * native_exp(-sigma_s.y * d)
                + sigma_s.z * native_exp(-sigma_s.z * d));
            *pdf = sigma > 0.f ? temp : 0.f;
            return d;
        }
    }
    
    return -1.f;
}

// Sample volume in order to find next scattering event
float Volume_GetDistancePdf(GLOBAL Volume const* volume, float dist)
{
    switch (volume->type)
    {
    case kHomogeneous:
    {
        float3 sigma_s = TEXTURED_INPUT_GET_COLOR(volume->sigma_s);
        return (1.f / 3.f) * (native_exp(-sigma_s.x * dist)
                            + native_exp(-sigma_s.y * dist)
                            + native_exp(-sigma_s.z * dist));
    }
    }

    return 0.f;
}

// Apply volume effects (absorbtion and emission) and scatter if needed.
// The rays we handling here might intersect something or miss, 
// since scattering can happen even for missed rays.
// That's why this function is called prior to ray compaction.
// In case ray has missed geometry (has shapeid < 0) and has been scattered,
// we put FAKE_SHAPE_SENTINEL into shapeid to prevent ray from being compacted away.
//
KERNEL void SampleVolume(
    // Ray batch
    GLOBAL ray const* rays,
    // Pixel indices
    GLOBAL int const* pixelindices,
    // Output indices
    GLOBAL int const* output_indices,
    // Number of rays
    GLOBAL int const* numrays,
    // Volumes
    GLOBAL Volume const* volumes,
    // Textures
    TEXTURE_ARG_LIST,
    // RNG seed
    uint rngseed,
    // Sampler state
    GLOBAL uint* random,
    // Sobol matrices
    GLOBAL uint const* sobol_mat,
    // Current bounce 
    int bounce,
    // Current frame
    int frame,
    // Intersection data
    GLOBAL Intersection* isects,
    // Current paths
    GLOBAL Path* paths,
    // Output
    GLOBAL float3* output
    )
{
    int globalid = get_global_id(0);

    // Only handle active rays
    if (globalid < *numrays)
    {
        int pixelidx = pixelindices[globalid];
        
        GLOBAL Path* path = paths + pixelidx;

        // Path can be dead here since compaction step has not 
        // yet been applied
        if (!Path_IsAlive(path))
            return;

        int volidx = Path_GetVolumeIdx(path);

        // Check if we are inside some volume
        if (volidx != -1)
        {
            Sampler sampler;
#if SAMPLER == SOBOL
            uint scramble = random[pixelidx] * 0x1fe3434f;
            Sampler_Init(&sampler, frame, SAMPLE_DIM_SURFACE_OFFSET + bounce * SAMPLE_DIMS_PER_BOUNCE + SAMPLE_DIM_VOLUME_APPLY_OFFSET, scramble);
#elif SAMPLER == RANDOM
            uint scramble = pixelidx * rngseed;
            Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
            uint rnd = random[pixelidx];
            uint scramble = rnd * 0x1fe3434f * ((frame + 71 * rnd) / (CMJ_DIM * CMJ_DIM));
            Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_SURFACE_OFFSET + bounce * SAMPLE_DIMS_PER_BOUNCE + SAMPLE_DIM_VOLUME_APPLY_OFFSET, scramble);
#endif

            // Try sampling volume for a next scattering event
            float pdf = 0.f;
            float maxdist = Intersection_GetDistance(isects + globalid);
            float2 sample = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
            float2 sample1 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
            float d = Volume_SampleDistance(&volumes[volidx], &rays[globalid], maxdist, make_float2(sample.x, sample1.y), &pdf);
            
            // Check if we shall skip the event (it is either outside of a volume or not happened at all)
            bool skip = d < 0 || d > maxdist || pdf <= 0.f;

            if (skip)
            {
                // In case we skip we just need to apply volume absorbtion and emission for the segment we went through
                // and clear scatter flag
                Path_ClearScatterFlag(path);
                // And finally update the throughput
                Path_MulThroughput(path, Volume_Transmittance(&volumes[volidx], &rays[globalid], maxdist) * Volume_GetDistancePdf(&volumes[volidx], maxdist));
                // Emission contribution accounting for a throughput we have so far
                Path_AddContribution(path, output, output_indices[pixelidx], Volume_Emission(&volumes[volidx], &rays[globalid], maxdist));
            }
            else
            {
                // Set scattering flag to notify ShadeVolume kernel to handle this path
                Path_SetScatterFlag(path);
                // Update the throughput
                float3 sigma_s = TEXTURED_INPUT_GET_COLOR(volumes[volidx].sigma_s);
                Path_MulThroughput(path, sigma_s * (Volume_Transmittance(&volumes[volidx], &rays[globalid], d) / pdf));
                // Emission contribution accounting for a throughput we have so far
                Path_AddContribution(path, output, output_indices[pixelidx], Volume_Emission(&volumes[volidx], &rays[globalid], d) / pdf);
                // Put fake shape to prevent from being compacted away
                isects[globalid].shapeid = FAKE_SHAPE_SENTINEL;
                // And keep scattering distance around as well
                isects[globalid].uvwt.w = d;
            }
        }
    }
}

#endif // VOLUMETRICS_CL
/**********************************************************************
Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************/
#ifndef VERTEX_CL
#define VERTEX_CL

// Path vertex type
enum PathVertexType
{
    kCamera,
    kSurface,
    kVolume,
    kLight
};

// Path vertex descriptor
typedef struct _PathVertex
{
    float3 position;
    float3 shading_normal;
    float3 geometric_normal;
    float2 uv;
    float pdf_forward;
    float pdf_backward;
    float3 flow;
    float3 unused;
    int type;
    int material_index;
    int flags;
    int padding;
} PathVertex;

// Initialize path vertex
INLINE
void PathVertex_Init(
    PathVertex* v, 
    float3 p, 
    float3 n,
    float3 ng, 
    float2 uv, 
    float pdf_fwd,
    float pdf_bwd, 
    float3 flow, 
    int type,
    int matidx
)
{
    v->position = p;
    v->shading_normal = n;
    v->geometric_normal = ng;
    v->uv = uv;
    v->pdf_forward = pdf_fwd;
    v->pdf_backward = pdf_bwd;
    v->flow = flow;
    v->type = type;
    v->material_index = matidx;
    v->flags = 0;
    v->unused = 0.f;
}

#endif


// Fill AOVs
KERNEL void FillAOVsUberV2(
    // Ray batch
    GLOBAL ray const* restrict rays,
    // Intersection data
    GLOBAL Intersection const* restrict isects,
    // Pixel indices
    GLOBAL int const* restrict pixel_idx,
    // Number of pixels
    GLOBAL int const* restrict num_items,
    // Vertices
    GLOBAL float3 const*restrict  vertices,
    // Normals
    GLOBAL float3 const* restrict normals,
    // UVs
    GLOBAL float2 const* restrict uvs,
    // Indices
    GLOBAL int const* restrict indices,
    // Shapes
    GLOBAL Shape const* restrict shapes,
    GLOBAL ShapeAdditionalData const* restrict shapes_additional,
    // Materials
    GLOBAL int const* restrict material_attributes,
    // Textures
    TEXTURE_ARG_LIST,
    // Environment texture index
    int env_light_idx,
    // Background texture index
    int background_idx,
    // Output size
    int width,
    int height,
    // Emissives
    GLOBAL Light const* restrict lights,
    // Number of emissive objects
    int num_lights,
    // RNG seed
    uint rngseed,
    // Sampler states
    GLOBAL uint* restrict random,
    // Sobol matrices
    GLOBAL uint const* restrict sobol_mat, 
    // Frame
    int frame,
    GLOBAL InputMapData const* restrict input_map_values,
    // World position flag
    int world_position_enabled, 
    // World position AOV
    GLOBAL float4* restrict aov_world_position,
    // World normal flag
    int world_shading_normal_enabled,
    // World normal AOV
    GLOBAL float4* restrict aov_world_shading_normal,
    // World true normal flag
    int world_geometric_normal_enabled,
    // World true normal AOV
    GLOBAL float4* restrict aov_world_geometric_normal,
    // UV flag
    int uv_enabled,
    // UV AOV
    GLOBAL float4* restrict aov_uv,
    // Wireframe flag
    int wireframe_enabled,
    // Wireframe AOV
    GLOBAL float4* restrict aov_wireframe,
    // Albedo flag
    int albedo_enabled,
    // Wireframe AOV
    GLOBAL float4* restrict aov_albedo,
    // World tangent flag
    int world_tangent_enabled,
    // World tangent AOV
    GLOBAL float4* restrict aov_world_tangent,
    // World bitangent flag
    int world_bitangent_enabled,
    // World bitangent AOV
    GLOBAL float4* restrict aov_world_bitangent,
    // Gloss enabled flag
    int gloss_enabled,
    // Specularity map
    GLOBAL float4* restrict aov_gloss,
    // Mesh_id enabled flag
    int mesh_id_enabled,
    // Mesh_id AOV
    GLOBAL float4* restrict mesh_id,
    // Group id enabled flag
    int group_id_enabled,
    // Group id AOV
    GLOBAL float4* restrict group_id,
    // Background enabled flag
    int background_enabled,
    // Background aov
    GLOBAL float4* restrict aov_background,
    // Depth enabled flag
    int depth_enabled,
    // Depth map
    GLOBAL float4* restrict aov_depth,
    // Shape id map enabled flag
    int shape_ids_enabled,
    // Shape id map stores shape ud in every pixel
    // And negative number if there is no any shape in the pixel
    GLOBAL float4* restrict aov_shape_ids
)
{
    int global_id = get_global_id(0);

    Scene scene =
    {
        vertices,
        normals,
        uvs,
        indices,
        shapes,
        material_attributes,
        input_map_values,
        lights,
        env_light_idx,
        num_lights
    };

    // Only applied to active rays after compaction
    if (global_id < *num_items)
    {
        Intersection isect = isects[global_id];
        int idx = pixel_idx[global_id];

        if (shape_ids_enabled)
            aov_shape_ids[idx].x = -1;

        if (background_enabled)
        {
            if (background_idx != -1)
            {
                float x = (float)(idx % width) / (float)width;
                float y = (float)(idx / width) / (float)height;
                float2 uv = make_float2(x, y);
                aov_background[idx].xyz += Texture_Sample2D(uv, TEXTURE_ARGS_IDX(background_idx)).xyz;
            }
            else if (env_light_idx != -1)
            {
                Light light = lights[env_light_idx];
                int tex = EnvironmentLight_GetBackgroundTexture(&light);
                if (tex != -1)
                {
                    aov_background[idx].xyz += light.multiplier * Texture_SampleEnvMap(rays[global_id].d.xyz, TEXTURE_ARGS_IDX(tex), light.ibl_mirror_x);
                }
            }
            aov_background[idx].w += 1.0f;
        }

        if (isect.shapeid > -1)
        {
            // Fetch incoming ray direction
            float3 wi = -normalize(rays[global_id].d.xyz);

            Sampler sampler;
#if SAMPLER == SOBOL 
            uint scramble = random[global_id] * 0x1fe3434f;
            Sampler_Init(&sampler, frame, SAMPLE_DIM_SURFACE_OFFSET, scramble);
#elif SAMPLER == RANDOM
            uint scramble = global_id * rngseed;
            Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
            uint rnd = random[global_id];
            uint scramble = rnd * 0x1fe3434f * ((frame + 331 * rnd) / (CMJ_DIM * CMJ_DIM));
            Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_SURFACE_OFFSET, scramble);
#endif

            // Fill surface data
            DifferentialGeometry diffgeo;
            Scene_FillDifferentialGeometry(&scene, &isect, &diffgeo);

            if (world_position_enabled)
            {
                aov_world_position[idx].xyz += diffgeo.p;
                aov_world_position[idx].w += 1.f;
            }

            if (world_shading_normal_enabled)
            {
                float ngdotwi = dot(diffgeo.ng, wi);
                bool backfacing = ngdotwi < 0.f;

                // Select BxDF
                UberV2ShaderData uber_shader_data;
                UberV2PrepareInputs(&diffgeo, input_map_values, material_attributes, TEXTURE_ARGS, &uber_shader_data);
                GetMaterialBxDFType(wi, &sampler, SAMPLER_ARGS, &diffgeo, &uber_shader_data);

                float s = Bxdf_IsBtdf(&diffgeo) ? (-sign(ngdotwi)) : 1.f;
                if (backfacing && !Bxdf_IsBtdf(&diffgeo))
                {
                    //Reverse normal and tangents in this case
                    //but not for BTDFs, since BTDFs rely
                    //on normal direction in order to arrange   
                    //indices of refraction
                    diffgeo.n = -diffgeo.n;
                    diffgeo.dpdu = -diffgeo.dpdu;
                    diffgeo.dpdv = -diffgeo.dpdv;
                }
                UberV2_ApplyShadingNormal(&diffgeo, &uber_shader_data);
                DifferentialGeometry_CalculateTangentTransforms(&diffgeo);

                aov_world_shading_normal[idx].xyz += diffgeo.n;
                aov_world_shading_normal[idx].w += 1.f;
            }

            if (world_geometric_normal_enabled)
            {
                aov_world_geometric_normal[idx].xyz += diffgeo.ng;
                aov_world_geometric_normal[idx].w += 1.f;
            }

            if (wireframe_enabled)
            {
                bool hit = (isect.uvwt.x < 1e-3) || (isect.uvwt.y < 1e-3) || (1.f - isect.uvwt.x - isect.uvwt.y < 1e-3);
                float3 value = hit ? make_float3(1.f, 1.f, 1.f) : make_float3(0.f, 0.f, 0.f);
                aov_wireframe[idx].xyz += value;
                aov_wireframe[idx].w += 1.f;
            }

            if (uv_enabled)
            {
                aov_uv[idx].xy += diffgeo.uv.xy;
                aov_uv[idx].w += 1.f;
            }

            if (albedo_enabled)
            {
                float ngdotwi = dot(diffgeo.ng, wi);
                bool backfacing = ngdotwi < 0.f;

                // Select BxDF
                UberV2ShaderData uber_shader_data;
                UberV2PrepareInputs(&diffgeo, input_map_values, material_attributes, TEXTURE_ARGS, &uber_shader_data);

                const float3 kd = ((diffgeo.mat.layers & kDiffuseLayer) == kDiffuseLayer) ?
                    uber_shader_data.diffuse_color.xyz : (float3)(0.0f);

                aov_albedo[idx].xyz += kd;
                aov_albedo[idx].w += 1.f;
            }

            if (world_tangent_enabled)
            {
                float ngdotwi = dot(diffgeo.ng, wi);
                bool backfacing = ngdotwi < 0.f;

                // Select BxDF
                UberV2ShaderData uber_shader_data;
                UberV2PrepareInputs(&diffgeo, input_map_values, material_attributes, TEXTURE_ARGS, &uber_shader_data);
                GetMaterialBxDFType(wi, &sampler, SAMPLER_ARGS, &diffgeo, &uber_shader_data);

                float s = Bxdf_IsBtdf(&diffgeo) ? (-sign(ngdotwi)) : 1.f;
                if (backfacing && !Bxdf_IsBtdf(&diffgeo))
                {
                    //Reverse normal and tangents in this case
                    //but not for BTDFs, since BTDFs rely
                    //on normal direction in order to arrange
                    //indices of refraction
                    diffgeo.n = -diffgeo.n;
                    diffgeo.dpdu = -diffgeo.dpdu;
                    diffgeo.dpdv = -diffgeo.dpdv;
                }

                UberV2_ApplyShadingNormal(&diffgeo, &uber_shader_data);
                DifferentialGeometry_CalculateTangentTransforms(&diffgeo);

                aov_world_tangent[idx].xyz += diffgeo.dpdu;
                aov_world_tangent[idx].w += 1.f;
            }

            if (world_bitangent_enabled)
            {
                float ngdotwi = dot(diffgeo.ng, wi);
                bool backfacing = ngdotwi < 0.f;

                // Select BxDF
                UberV2ShaderData uber_shader_data;
                UberV2PrepareInputs(&diffgeo, input_map_values, material_attributes, TEXTURE_ARGS, &uber_shader_data);
                GetMaterialBxDFType(wi, &sampler, SAMPLER_ARGS, &diffgeo, &uber_shader_data);

                float s = Bxdf_IsBtdf(&diffgeo) ? (-sign(ngdotwi)) : 1.f;
                if (backfacing && !Bxdf_IsBtdf(&diffgeo))
                {
                    //Reverse normal and tangents in this case
                    //but not for BTDFs, since BTDFs rely
                    //on normal direction in order to arrange
                    //indices of refraction
                    diffgeo.n = -diffgeo.n;
                    diffgeo.dpdu = -diffgeo.dpdu;
                    diffgeo.dpdv = -diffgeo.dpdv;
                }

                UberV2_ApplyShadingNormal(&diffgeo, &uber_shader_data);
                DifferentialGeometry_CalculateTangentTransforms(&diffgeo);

                aov_world_bitangent[idx].xyz += diffgeo.dpdv;
                aov_world_bitangent[idx].w += 1.f;
            }

            if (gloss_enabled)
            {
                float ngdotwi = dot(diffgeo.ng, wi);
                bool backfacing = ngdotwi < 0.f;

                // Select BxDF
                UberV2ShaderData uber_shader_data;
                UberV2PrepareInputs(&diffgeo, input_map_values, material_attributes, TEXTURE_ARGS, &uber_shader_data);
                GetMaterialBxDFType(wi, &sampler, SAMPLER_ARGS, &diffgeo, &uber_shader_data);

                float gloss = 0.f;
                if ((diffgeo.mat.layers & kCoatingLayer) == kCoatingLayer)
                {
                    gloss = 1.0f;
                }
                else if ((diffgeo.mat.layers & kReflectionLayer) == kReflectionLayer)
                {
                    gloss = 1.0f - uber_shader_data.reflection_roughness;
                    if ((diffgeo.mat.layers & kRefractionLayer) == kRefractionLayer)
                    {
                        gloss = max(gloss, 1.0f - uber_shader_data.refraction_roughness);
                    }
                }
                else if ((diffgeo.mat.layers & kRefractionLayer) == kRefractionLayer)
                {
                    gloss = 1.0f - uber_shader_data.refraction_roughness;
                }

                aov_gloss[idx].xyz += gloss;
                aov_gloss[idx].w += 1.f;
            }
            
            if (mesh_id_enabled)
            {
                Sampler shapeid_sampler;
                shapeid_sampler.index = shapes[isect.shapeid - 1].id;
                mesh_id[idx].xyz += clamp(make_float3(UniformSampler_Sample1D(&shapeid_sampler),
                    UniformSampler_Sample1D(&shapeid_sampler),
                    UniformSampler_Sample1D(&shapeid_sampler)), 0.0f, 1.0f);
                mesh_id[idx].w += 1.0f;
            }

            if (group_id_enabled)
            {
                Sampler groupid_sampler;
                groupid_sampler.index = shapes_additional[isect.shapeid - 1].group_id;
                group_id[idx].xyz += clamp(make_float3(UniformSampler_Sample1D(&groupid_sampler),
                    UniformSampler_Sample1D(&groupid_sampler),
                    UniformSampler_Sample1D(&groupid_sampler)), 0.0f, 1.0f);
                group_id[idx].w += 1.0f;
            }

            if (depth_enabled)
            {
                float w = aov_depth[idx].w;
                if (w == 0.f)
                {
                    aov_depth[idx].xyz = isect.uvwt.w;
                    aov_depth[idx].w = 1.f;
                }
                else
                {
                    aov_depth[idx].xyz += isect.uvwt.w;
                    aov_depth[idx].w += 1.f;
                }
            }

            if (shape_ids_enabled)
            {
                aov_shape_ids[idx].x = shapes[isect.shapeid - 1].id;
            }
        }
    }
}



#endif // MONTE_CARLO_RENDERER_CL
