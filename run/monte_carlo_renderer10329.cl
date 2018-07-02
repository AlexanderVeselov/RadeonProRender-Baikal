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
#ifndef VOLUMETRICS_CL
#define VOLUMETRICS_CL
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


// Pinhole camera implementation.
// This kernel is being used if aperture value = 0.
KERNEL
void PerspectiveCamera_GeneratePaths(
    // Camera
    GLOBAL Camera const* restrict camera, 
    // Image resolution
    int output_width,
    int output_height,
    // Pixel domain buffer
    GLOBAL int const* restrict pixel_idx,
    // Size of pixel domain buffer
    GLOBAL int const* restrict num_pixels,
    // RNG seed value
    uint rng_seed,
    // Current frame
    uint frame,
    // Rays to generate
    GLOBAL ray* restrict rays,
    // RNG data
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat
)
{
    int global_id = get_global_id(0);

    // Check borders
    if (global_id < *num_pixels)
    {
        int idx = pixel_idx[global_id];
        int y = idx / output_width;
        int x = idx % output_width;

        // Get pointer to ray & path handles
        GLOBAL ray* my_ray = rays + global_id;

        // Initialize sampler
        Sampler sampler;
#if SAMPLER == SOBOL
        uint scramble = random[x + output_width * y] * 0x1fe3434f;

        if (frame & 0xF)
        {
            random[x + output_width * y] = WangHash(scramble);
        }

        Sampler_Init(&sampler, frame, SAMPLE_DIM_CAMERA_OFFSET, scramble);
#elif SAMPLER == RANDOM
        uint scramble = x + output_width * y * rng_seed;
        Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
        uint rnd = random[x + output_width * y];
        uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
        Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_CAMERA_OFFSET, scramble);
#endif

        // Generate sample
#ifndef BAIKAL_GENERATE_SAMPLE_AT_PIXEL_CENTER
        float2 sample0 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
#else
        float2 sample0 = make_float2(0.5f, 0.5f);
#endif

        // Calculate [0..1] image plane sample
        float2 img_sample;
        img_sample.x = (float)x / output_width + sample0.x / output_width;
        img_sample.y = (float)y / output_height + sample0.y / output_height;

        // Transform into [-0.5, 0.5]
        float2 h_sample = img_sample - make_float2(0.5f, 0.5f);
        // Transform into [-dim/2, dim/2]
        float2 c_sample = h_sample * camera->dim;

        // Calculate direction to image plane
        my_ray->d.xyz = normalize(camera->focal_length * camera->forward + c_sample.x * camera->right + c_sample.y * camera->up);
        // Origin == camera position + nearz * d
        my_ray->o.xyz = camera->p + camera->zcap.x * my_ray->d.xyz;
        // Max T value = zfar - znear since we moved origin to znear
        my_ray->o.w = camera->zcap.y - camera->zcap.x;
        // Generate random time from 0 to 1
        my_ray->d.w = sample0.x;
        // Set ray max
        my_ray->extra.x = 0xFFFFFFFF;
        my_ray->extra.y = 0xFFFFFFFF;
        Ray_SetExtra(my_ray, 1.f);
        Ray_SetMask(my_ray, VISIBILITY_MASK_PRIMARY);
    }
}

// Physical camera implemenation.
// This kernel is being used if aperture > 0.
KERNEL void PerspectiveCameraDof_GeneratePaths(
    // Camera
    GLOBAL Camera const* restrict camera,
    // Image resolution
    int output_width,
    int output_height,
    // Pixel domain buffer
    GLOBAL int const* restrict pixel_idx,
    // Size of pixel domain buffer
    GLOBAL int const* restrict num_pixels,
    // RNG seed value
    uint rng_seed,
    // Current frame
    uint frame,
    // Rays to generate
    GLOBAL ray* restrict rays,
    // RNG data
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat
)
{
    int global_id = get_global_id(0);

    // Check borders
    if (global_id < *num_pixels)
    {
        int idx = pixel_idx[global_id];
        int y = idx / output_width;
        int x = idx % output_width;

        // Get pointer to ray & path handles
        GLOBAL ray* my_ray = rays + global_id;

        // Initialize sampler
        Sampler sampler;
#if SAMPLER == SOBOL
        uint scramble = random[x + output_width * y] * 0x1fe3434f;

        if (frame & 0xF)
        {
            random[x + output_width * y] = WangHash(scramble);
        }

        Sampler_Init(&sampler, frame, SAMPLE_DIM_CAMERA_OFFSET, scramble);
#elif SAMPLER == RANDOM
        uint scramble = x + output_width * y * rng_seed;
        Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
        uint rnd = random[x + output_width * y];
        uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
        Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_CAMERA_OFFSET, scramble);
#endif

        // Generate pixel and lens samples
#ifndef BAIKAL_GENERATE_SAMPLE_AT_PIXEL_CENTER
        float2 sample0 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
#else
        float2 sample0 = make_float2(0.5f, 0.5f);
#endif
        float2 sample1 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);

        // Calculate [0..1] image plane sample
        float2 img_sample;
        img_sample.x = (float)x / output_width + sample0.x / output_width;
        img_sample.y = (float)y / output_height + sample0.y / output_height;

        // Transform into [-0.5, 0.5]
        float2 h_sample = img_sample - make_float2(0.5f, 0.5f);
        // Transform into [-dim/2, dim/2]
        float2 c_sample = h_sample * camera->dim;

        // Generate sample on the lens
        float2 lens_sample = camera->aperture * Sample_MapToDiskConcentric(sample1);
        // Calculate position on focal plane
        float2 focal_plane_sample = c_sample * camera->focus_distance / camera->focal_length;
        // Calculate ray direction
        float2 camera_dir = focal_plane_sample - lens_sample;

        // Calculate direction to image plane
        my_ray->d.xyz = normalize(camera->forward * camera->focus_distance + camera->right * camera_dir.x + camera->up * camera_dir.y);
        // Origin == camera position + nearz * d
        my_ray->o.xyz = camera->p + lens_sample.x * camera->right + lens_sample.y * camera->up;
        // Max T value = zfar - znear since we moved origin to znear
        my_ray->o.w = camera->zcap.y - camera->zcap.x;
        // Generate random time from 0 to 1
        my_ray->d.w = sample0.x;
        // Set ray max
        my_ray->extra.x = 0xFFFFFFFF;
        my_ray->extra.y = 0xFFFFFFFF;
        Ray_SetExtra(my_ray, 1.f);
        Ray_SetMask(my_ray, VISIBILITY_MASK_PRIMARY);
    }
}


KERNEL
void PerspectiveCamera_GenerateVertices(
    // Camera
    GLOBAL Camera const* restrict camera,
    // Image resolution
    int output_width,
    int output_height,
    // Pixel domain buffer
    GLOBAL int const* restrict pixel_idx,
    // Size of pixel domain buffer
    GLOBAL int const* restrict num_pixels,
    // RNG seed value
    uint rng_seed,
    // Current frame
    uint frame,
    // RNG data
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat,
    // Rays to generate
    GLOBAL ray* restrict rays,
    // Eye subpath vertices
    GLOBAL PathVertex* restrict eye_subpath,
    // Eye subpath length
    GLOBAL int* restrict eye_subpath_length,
    // Path buffer
    GLOBAL Path* restrict paths
)

{
    int global_id = get_global_id(0);

    // Check borders
    if (global_id < *num_pixels)
    {
        int idx = pixel_idx[global_id];
        int y = idx / output_width;
        int x = idx % output_width;

        // Get pointer to ray & path handles
        GLOBAL ray* my_ray = rays + global_id;

        GLOBAL PathVertex* my_vertex = eye_subpath + BDPT_MAX_SUBPATH_LEN * idx;
        GLOBAL int* my_count = eye_subpath_length + idx;
        GLOBAL Path* my_path = paths + idx;

        // Initialize sampler
        Sampler sampler;
#if SAMPLER == SOBOL
        uint scramble = random[x + output_width * y] * 0x1fe3434f;

        if (frame & 0xF)
        {
            random[x + output_width * y] = WangHash(scramble);
        }

        Sampler_Init(&sampler, frame, SAMPLE_DIM_CAMERA_OFFSET, scramble);
#elif SAMPLER == RANDOM
        uint scramble = x + output_width * y * rng_seed;
        Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
        uint rnd = random[x + output_width * y];
        uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
        Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_CAMERA_OFFSET, scramble);
#endif

        // Generate sample
        float2 sample0 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);

        // Calculate [0..1] image plane sample
        float2 img_sample;
        img_sample.x = (float)x / output_width + sample0.x / output_width;
        img_sample.y = (float)y / output_height + sample0.y / output_height;

        // Transform into [-0.5, 0.5]
        float2 h_sample = img_sample - make_float2(0.5f, 0.5f);
        // Transform into [-dim/2, dim/2]
        float2 c_sample = h_sample * camera->dim;

        // Calculate direction to image plane
        my_ray->d.xyz = normalize(camera->focal_length * camera->forward + c_sample.x * camera->right + c_sample.y * camera->up);
        // Origin == camera position + nearz * d
        my_ray->o.xyz = camera->p + camera->zcap.x * my_ray->d.xyz;
        // Max T value = zfar - znear since we moved origin to znear
        my_ray->o.w = camera->zcap.y - camera->zcap.x;
        // Generate random time from 0 to 1
        my_ray->d.w = sample0.x;
        // Set ray max
        my_ray->extra.x = 0xFFFFFFFF;
        my_ray->extra.y = 0xFFFFFFFF;
        Ray_SetExtra(my_ray, 1.f);

        PathVertex v;
        PathVertex_Init(&v,
            camera->p,
            camera->forward,
            camera->forward,
            0.f,
            1.f,
            1.f,
            1.f,
            kCamera,
            -1);

        *my_count = 1;
        *my_vertex = v;

        // Initlize path data
        my_path->throughput = make_float3(1.f, 1.f, 1.f);
        my_path->volume = -1;
        my_path->flags = 0;
        my_path->active = 0xFF;
    }
}

KERNEL
void PerspectiveCameraDof_GenerateVertices(
    // Camera
    GLOBAL Camera const* restrict camera,
    // Image resolution
    int output_width,
    int output_height,
    // Pixel domain buffer
    GLOBAL int const* restrict pixel_idx,
    // Size of pixel domain buffer
    GLOBAL int const* restrict num_pixels,
    // RNG seed value
    uint rng_seed,
    // Current frame
    uint frame,
    // RNG data
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat,
    // Rays to generate
    GLOBAL ray* restrict rays,
    // Eye subpath vertices
    GLOBAL PathVertex* restrict eye_subpath,
    // Eye subpath length
    GLOBAL int* restrict eye_subpath_length,
    // Path buffer
    GLOBAL Path* restrict paths
)

{
    int global_id = get_global_id(0);

    // Check borders
    if (global_id < *num_pixels)
    {
        int idx = pixel_idx[global_id];
        int y = idx / output_width;
        int x = idx % output_width;

        // Get pointer to ray & path handles
        GLOBAL ray* my_ray = rays + global_id;
        GLOBAL PathVertex* my_vertex = eye_subpath + BDPT_MAX_SUBPATH_LEN * (y * output_width + x);
        GLOBAL int* my_count = eye_subpath_length + y * output_width + x;
        GLOBAL Path* my_path = paths + y * output_width + x;

        // Initialize sampler
        Sampler sampler;
#if SAMPLER == SOBOL
        uint scramble = random[x + output_width * y] * 0x1fe3434f;

        if (frame & 0xF)
        {
            random[x + output_width * y] = WangHash(scramble);
        }

        Sampler_Init(&sampler, frame, SAMPLE_DIM_CAMERA_OFFSET, scramble);
#elif SAMPLER == RANDOM
        uint scramble = x + output_width * y * rng_seed;
        Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
        uint rnd = random[x + output_width * y];
        uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
        Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_CAMERA_OFFSET, scramble);
#endif

        // Generate pixel and lens samples
        float2 sample0 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
        float2 sample1 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);

        // Calculate [0..1] image plane sample
        float2 img_sample;
        img_sample.x = (float)x / output_width + sample0.x / output_width;
        img_sample.y = (float)y / output_height + sample0.y / output_height;

        // Transform into [-0.5, 0.5]
        float2 h_sample = img_sample - make_float2(0.5f, 0.5f);
        // Transform into [-dim/2, dim/2]
        float2 c_sample = h_sample * camera->dim;

        // Generate sample on the lens
        float2 lens_sample = camera->aperture * Sample_MapToDiskConcentric(sample1);
        // Calculate position on focal plane
        float2 focal_plane_sample = c_sample * camera->focus_distance / camera->focal_length;
        // Calculate ray direction
        float2 camera_dir = focal_plane_sample - lens_sample;

        // Calculate direction to image plane
        my_ray->d.xyz = normalize(camera->forward * camera->focus_distance + camera->right * camera_dir.x + camera->up * camera_dir.y);
        // Origin == camera position + nearz * d
        my_ray->o.xyz = camera->p + lens_sample.x * camera->right + lens_sample.y * camera->up;
        // Max T value = zfar - znear since we moved origin to znear
        my_ray->o.w = camera->zcap.y - camera->zcap.x;
        // Generate random time from 0 to 1
        my_ray->d.w = sample0.x;
        // Set ray max
        my_ray->extra.x = 0xFFFFFFFF;
        my_ray->extra.y = 0xFFFFFFFF;
        Ray_SetExtra(my_ray, 1.f);

        PathVertex v;
        PathVertex_Init(&v,
            camera->p,
            camera->forward,
            camera->forward,
            0.f,
            1.f,
            1.f,
            1.f,
            kCamera,
            -1);

        *my_count = 1;
        *my_vertex = v;

        // Initlize path data
        my_path->throughput = make_float3(1.f, 1.f, 1.f);
        my_path->volume = -1;
        my_path->flags = 0;
        my_path->active = 0xFF;
    }
}

uint Part1By1(uint x)
{
    x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

uint Morton2D(uint x, uint y)
{
    return (Part1By1(y) << 1) + Part1By1(x);
}

KERNEL void GenerateTileDomain(
    int output_width,
    int output_height,
    int offset_x,
    int offset_y,
    int width,
    int height,
    uint rng_seed,
    uint frame,
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat,
    GLOBAL int* restrict indices,
    GLOBAL int* restrict count
)
{
    int2 global_id;
    global_id.x = get_global_id(0);
    global_id.y = get_global_id(1);

    int2 local_id;
    local_id.x = get_local_id(0);
    local_id.y = get_local_id(1);

    int2 group_id;
    group_id.x = get_group_id(0);
    group_id.y = get_group_id(1);

    int2 tile_size;
    tile_size.x = get_local_size(0);
    tile_size.y = get_local_size(1);

    int num_tiles_x = output_width / tile_size.x;
    int num_tiles_y = output_height / tile_size.y;

    int start_idx = output_width * offset_y + offset_x;

    if (global_id.x < width && global_id.y < height)
    {
        int idx = start_idx +
            (group_id.y * tile_size.y + local_id.y) * output_width +
            (group_id.x * tile_size.x + local_id.x);

        indices[global_id.y * width + global_id.x] = idx;
    }

    if (global_id.x == 0 && global_id.y == 0)
    {
        *count = width * height;
    }
}

KERNEL void GenerateTileDomain_Adaptive(
    int output_width,
    int output_height,
    int offset_x,
    int offset_y,
    int width,
    int height,
    uint rng_seed,
    uint frame,
    GLOBAL uint* restrict random,
    GLOBAL uint const* restrict sobol_mat,
    GLOBAL int const* restrict tile_distribution,
    GLOBAL int* restrict indices,
    GLOBAL int* restrict count
)
{
    int2 global_id;
    global_id.x = get_global_id(0);
    global_id.y = get_global_id(1);

    int2 local_id;
    local_id.x = get_local_id(0);
    local_id.y = get_local_id(1);

    int2 group_id;
    group_id.x = get_group_id(0);
    group_id.y = get_group_id(1);

    int2 tile_size;
    tile_size.x = get_local_size(0);
    tile_size.y = get_local_size(1);


    // Initialize sampler  
    Sampler sampler;
    int x = global_id.x;
    int y = global_id.y;
#if SAMPLER == SOBOL
    uint scramble = random[x + output_width * y] * 0x1fe3434f;

    if (frame & 0xF)
    {
        random[x + output_width * y] = WangHash(scramble);
    }

    Sampler_Init(&sampler, frame, SAMPLE_DIM_IMG_PLANE_EVALUATE_OFFSET, scramble);
#elif SAMPLER == RANDOM
    uint scramble = x + output_width * y * rng_seed;
    Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
    uint rnd = random[group_id.x + output_width *group_id.y];
    uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
    Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_IMG_PLANE_EVALUATE_OFFSET, scramble);
#endif

    float2 sample = Sampler_Sample2D(&sampler, SAMPLER_ARGS);

    float pdf;
    int tile = Distribution1D_SampleDiscrete(sample.x, tile_distribution, &pdf);

    int num_tiles_x = output_width / tile_size.x;
    int num_tiles_y = output_height / tile_size.y;

    int tile_y = clamp(tile / num_tiles_x , 0, num_tiles_y - 1);
    int tile_x = clamp(tile % num_tiles_x, 0, num_tiles_x - 1);

    int start_idx = output_width * offset_y + offset_x;

    if (global_id.x < width && global_id.y < height)
    {
        int idx = start_idx +
            (tile_y * tile_size.y + local_id.y) * output_width +
            (tile_x * tile_size.x + local_id.x);

        indices[global_id.y * width + global_id.x] = idx;
    }

    if (global_id.x == 0 && global_id.y == 0)
    {
        *count = width * height;
    }
}

// Copy data to interop texture if supported
KERNEL void AccumulateData(
    GLOBAL float4 const* src_data,
    int num_elements,
    GLOBAL float4* dst_data
)
{
    int global_id = get_global_id(0);

    if (global_id < num_elements)
    {
        float4 v = src_data[global_id];
        dst_data[global_id] += v;
    }
}

//#define ADAPTIVITY_DEBUG
// Copy data to interop texture if supported
KERNEL void ApplyGammaAndCopyData(
    GLOBAL float4 const* data,
    int img_width,
    int img_height,
    float gamma,
    write_only image2d_t img
)
{
    int global_id = get_global_id(0);

    int global_idx = global_id % img_width;
    int global_idy = global_id / img_width;

    if (global_idx < img_width && global_idy < img_height)
    {
        float4 v = data[global_id];
#ifdef ADAPTIVITY_DEBUG
        float a = v.w < 1024 ? min(1.f, v.w / 1024.f) : 0.f;
        float4 mul_color = make_float4(1.f, 1.f - a, 1.f - a, 1.f);
        v *= mul_color;
#endif

        float4 val = clamp(native_powr(v / v.w, 1.f / gamma), 0.f, 1.f);
        write_imagef(img, make_int2(global_idx, global_idy), val);
    }
} 

KERNEL void AccumulateSingleSample(
    GLOBAL float4 const* restrict src_sample_data,
    GLOBAL float4* restrict dst_accumulation_data,
    GLOBAL int* restrict scatter_indices,
    int num_elements
)
{
    int global_id = get_global_id(0);

    if (global_id < num_elements)
    {
        int idx = scatter_indices[global_id];
        float4 sample = src_sample_data[global_id];
        dst_accumulation_data[idx].xyz += sample.xyz;
        dst_accumulation_data[idx].w += 1.f;
    }
}

INLINE void group_reduce_add(__local float* lds, int size, int lid)
{
    for (int offset = (size >> 1); offset > 0; offset >>= 1)
    {
        if (lid < offset)
        {
            lds[lid] += lds[lid + offset];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
}

INLINE void group_reduce_min(__local float* lds, int size, int lid)
{
    for (int offset = (size >> 1); offset > 0; offset >>= 1)
    {
        if (lid < offset)
        {
            lds[lid] = min(lds[lid], lds[lid + offset]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
}


INLINE void group_reduce_max(__local float* lds, int size, int lid)
{
    for (int offset = (size >> 1); offset > 0; offset >>= 1)
    {
        if (lid < offset)
        {
            lds[lid] = max(lds[lid], lds[lid + offset]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
}


KERNEL void EstimateVariance(
    GLOBAL float4 const* restrict image_buffer,
    GLOBAL float* restrict variance_buffer,
    int width,
    int height
)
{
    __local float lds[256];

    int x = get_global_id(0);
    int y = get_global_id(1);
    int lx = get_local_id(0);
    int ly = get_local_id(1);
    int gx = get_group_id(0);
    int gy = get_group_id(1);
    int wx = get_local_size(0);
    int wy = get_local_size(1);
    int num_tiles = (width + wx - 1) / wx;
    int lid = ly * wx + lx;

    float value = 0.f;
    if (x < width && y < height)
    {
        float4 rw = image_buffer[y * width + x]; rw /= rw.w;
        value = 4*luminance(clamp(rw.xyz, 0.f, 1.f));
        rw = y + 1 < height ? image_buffer[(y + 1) * width + x] : image_buffer[y * width + x]; rw /= rw.w;
        value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        rw = y - 1 >= 0 ? image_buffer[(y - 1) * width + x] : image_buffer[y * width + x]; rw /= rw.w;
        value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        rw = x + 1 < width ? image_buffer[y * width + x + 1] : image_buffer[y * width + x]; rw /= rw.w;
        value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        rw = x - 1 >= 0 ? image_buffer[y * width + x - 1] : image_buffer[y * width + x]; rw /= rw.w;
        value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        //rw = y + 1 < height && x + 1 < width ? image_buffer[(y + 1) * width + x + 1] : image_buffer[y * width + x]; rw /= rw.w;
        //value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        //rw = y - 1 >= 0 && x - 1 >= 0 ? image_buffer[(y - 1) * width + x - 1] : image_buffer[y * width + x]; rw /= rw.w;
        //value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        //rw = y + 1 < height && x - 1 >= 0 ? image_buffer[(y + 1) * width + x - 1] : image_buffer[y * width + x]; rw /= rw.w;
        //value -= luminance(clamp(rw.xyz, 0.f, 1.f));
        //rw = y - 1 >= 0 && x + 1 < width ? image_buffer[(y - 1) * width + x + 1] : image_buffer[y * width + x]; rw /= rw.w;
        //value -= luminance(clamp(rw.xyz, 0.f, 1.f));
    }

    value = fabs(value);
    lds[lid] = value;
    barrier(CLK_LOCAL_MEM_FENCE);

    group_reduce_add(lds, 256, lid);

    float mean = lds[0] / (wx * wy);
    barrier(CLK_LOCAL_MEM_FENCE);

    /*lds[lid] = (mean - value) * (mean - value);
    barrier(CLK_LOCAL_MEM_FENCE);

    group_reduce_add(lds, 256, lid);*/

    if (x < width && y < height)
    {
        if (lx == 0 && ly == 0)
        {
            //float dev = lds[0] / (wx * wy - 1);
            variance_buffer[gy * num_tiles + gx] = mean;
        }
    }
}

KERNEL
void  OrthographicCamera_GeneratePaths(
                                     // Camera
                                     GLOBAL Camera const* restrict camera,
                                     // Image resolution
                                     int output_width,
                                     int output_height,
                                     // Pixel domain buffer
                                     GLOBAL int const* restrict pixel_idx,
                                     // Size of pixel domain buffer
                                     GLOBAL int const* restrict num_pixels,
                                     // RNG seed value
                                     uint rng_seed,
                                     // Current frame
                                     uint frame,
                                     // Rays to generate
                                     GLOBAL ray* restrict rays,
                                     // RNG data
                                     GLOBAL uint* restrict random,
                                     GLOBAL uint const* restrict sobol_mat
                                     )
{
    int global_id = get_global_id(0);
    
    // Check borders
    if (global_id < *num_pixels)
    {
        int idx = pixel_idx[global_id];
        int y = idx / output_width;
        int x = idx % output_width;
        
        // Get pointer to ray & path handles
        GLOBAL ray* my_ray = rays + global_id;
        
        // Initialize sampler
        Sampler sampler;
#if SAMPLER == SOBOL
        uint scramble = random[x + output_width * y] * 0x1fe3434f;
        
        if (frame & 0xF)
        {
            random[x + output_width * y] = WangHash(scramble);
        }
        
        Sampler_Init(&sampler, frame, SAMPLE_DIM_CAMERA_OFFSET, scramble);
#elif SAMPLER == RANDOM
        uint scramble = x + output_width * y * rng_seed;
        Sampler_Init(&sampler, scramble);
#elif SAMPLER == CMJ
        uint rnd = random[x + output_width * y];
        uint scramble = rnd * 0x1fe3434f * ((frame + 133 * rnd) / (CMJ_DIM * CMJ_DIM));
        Sampler_Init(&sampler, frame % (CMJ_DIM * CMJ_DIM), SAMPLE_DIM_CAMERA_OFFSET, scramble);
#endif
        
        // Generate sample
#ifndef BAIKAL_GENERATE_SAMPLE_AT_PIXEL_CENTER
        float2 sample0 = Sampler_Sample2D(&sampler, SAMPLER_ARGS);
#else
        float2 sample0 = make_float2(0.5f, 0.5f);
#endif
        
        // Calculate [0..1] image plane sample
        float2 img_sample;
        img_sample.x = (float)x / output_width + sample0.x / output_width;
        img_sample.y = (float)y / output_height + sample0.y / output_height;
        
        // Transform into [-0.5, 0.5]
        float2 h_sample = img_sample - make_float2(0.5f, 0.5f);
        // Transform into [-dim/2, dim/2]
        float2 c_sample = h_sample * camera->dim;
        
        // Calculate direction to image plane
        my_ray->d.xyz = normalize(camera->forward);
        // Origin == camera position + nearz * d
        my_ray->o.xyz = camera->p + c_sample.x * camera->right + c_sample.y * camera->up;
        // Max T value = zfar - znear since we moved origin to znear
        my_ray->o.w = camera->zcap.y - camera->zcap.x;
        // Generate random time from 0 to 1
        my_ray->d.w = sample0.x;
        // Set ray max
        my_ray->extra.x = 0xFFFFFFFF;
        my_ray->extra.y = 0xFFFFFFFF;
        Ray_SetExtra(my_ray, 1.f);
        Ray_SetMask(my_ray, VISIBILITY_MASK_PRIMARY);
    }
}

///< Illuminate missing rays
KERNEL void ShadeBackgroundImage(
    // Ray batch
    GLOBAL ray const* restrict rays,
    // Intersection data
    GLOBAL Intersection const* restrict isects,
    // Pixel indices
    GLOBAL int const* restrict pixel_indices,
    // Output indices
    GLOBAL int const*  restrict output_indices,
    // Number of rays
    int num_rays,
    int background_idx,
    // Output size
    int width,
    int height,
    // Textures
    TEXTURE_ARG_LIST,
    // Output values
    GLOBAL float4* restrict output
)
{
    int global_id = get_global_id(0);

    if (global_id < num_rays)
    {
        int pixel_idx = pixel_indices[global_id];
        int output_index = output_indices[pixel_idx];

        float x = (float)(output_index % width) / (float)width;
        float y = (float)(output_index / width) / (float)height;

        float4 v = make_float4(0.f, 0.f, 0.f, 1.f);

        // In case of a miss
        if (isects[global_id].shapeid < 0)
        {
            float2 uv = make_float2(x, y);
            v.xyz = Texture_Sample2D(uv, TEXTURE_ARGS_IDX(background_idx)).xyz;
        }
        
        ADD_FLOAT4(&output[output_index], v);
    }
}


#endif // MONTE_CARLO_RENDERER_CL
