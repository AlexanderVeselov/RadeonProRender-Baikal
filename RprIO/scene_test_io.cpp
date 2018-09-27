#define _USE_MATH_DEFINES
#include <cmath>

#include "scene_io.h"
#include <vector>
#include <cstring>

class SceneIoTest : public SceneIo::Loader
{
public:
    SceneIoTest();
    ~SceneIoTest() = default;

    rpr_int LoadScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const override;

    rpr_int SaveScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const override;
private:
    rpr_int AddSphere(rpr_context context, rpr_scene scene, rpr_material_node material, rpr_char const* name,
        rpr_uint lat, rpr_uint lon, rpr_float x, rpr_float y, rpr_float z, rpr_float r, rpr_shape* out_sphere) const;

    rpr_int AddPlane(rpr_context context, rpr_scene scene, rpr_material_node material, rpr_char const* name,
        rpr_float x, rpr_float y, rpr_float z, rpr_float width, rpr_float height, rpr_float nx, rpr_float ny, rpr_float nz, rpr_shape* out_plane) const;

    rpr_int AddIbl(rpr_context context, rpr_scene scene, rpr_char const* image_filename, rpr_light* out_light) const;

    rpr_int AddSpotLight(rpr_context context, rpr_scene scene, rpr_float pitch, rpr_float yaw, rpr_float x, rpr_float y, rpr_float z,
        rpr_float r, rpr_float g, rpr_float b, rpr_float iangle, rpr_float oangle, rpr_light* out_light) const;

    rpr_int CreateDiffuseMaterial(rpr_material_system material_system, rpr_char const* name,
        rpr_float r, rpr_float g, rpr_float b, rpr_material_node* out_material) const;

    rpr_int CreateEmissiveMaterial(rpr_material_system material_system, rpr_char const* name,
        rpr_float r, rpr_float g, rpr_float b, rpr_material_node* out_material) const;

};

// Create an instance of test loader
static SceneIoTest test_loader;

inline SceneIoTest::SceneIoTest()
    : SceneIo::Loader("test", this)
{}

rpr_int SceneIoTest::CreateDiffuseMaterial(rpr_material_system material_system, rpr_char const* name,
    rpr_float r, rpr_float g, rpr_float b, rpr_material_node* out_material) const
{
    assert(out_material);

    rpr_material_node material = nullptr;
    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &material);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_DIFFUSE);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, r, g, b, 1.0f);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(material, name);
    RETURN_IF_FAILED(status);

    *out_material = material;

    return RPR_SUCCESS;
}

rpr_int SceneIoTest::CreateEmissiveMaterial(rpr_material_system material_system, rpr_char const* name,
    rpr_float r, rpr_float g, rpr_float b, rpr_material_node* out_material) const
{
    assert(out_material);

    rpr_material_node material = nullptr;
    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &material);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_EMISSION);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_EMISSION_COLOR, r, g, b, 1.0f);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(material, name);
    RETURN_IF_FAILED(status);

    *out_material = material;

    return RPR_SUCCESS;
}

rpr_int SceneIoTest::AddSphere(rpr_context context, rpr_scene scene, rpr_material_node material, rpr_char const* name,
    rpr_uint lat, rpr_uint lon, rpr_float x, rpr_float y, rpr_float z, rpr_float r, rpr_shape* out_sphere) const
{
    std::size_t num_verts = (lat - 2) * lon + 2;
    std::size_t num_tris = (lat - 2) * (lon - 1) * 2;

    std::vector<rpr_float> vertices(num_verts * 3);
    std::vector<rpr_float> normals(num_verts * 3);
    std::vector<rpr_float> uvs(num_verts * 2);
    std::vector<rpr_int> indices(num_tris * 3);

    std::size_t t = 0;
    for (std::uint32_t j = 1; j < lat - 1; j++)
    {
        for (std::uint32_t i = 0; i < lon; i++)
        {
            float theta = float(j) / (lat - 1) * (float)M_PI;
            float phi = float(i) / (lon - 1) * (float)M_PI * 2.0f;
            vertices[t * 3] = r * sinf(theta) * cosf(phi) + x;
            vertices[t * 3 + 1] = r * cosf(theta) + y;
            vertices[t * 3 + 2] = r * -sinf(theta) * sinf(phi) + z;

            normals[t * 3] = sinf(theta) * cosf(phi);
            normals[t * 3 + 1] = cosf(theta);
            normals[t * 3 + 2] = -sinf(theta) * sinf(phi);

            uvs[t * 2] = phi / (2 * (float)M_PI);
            uvs[t * 2 + 1] = theta / ((float)M_PI);

            ++t;
        }
    }

    // North pole of the sphere
    {
        vertices[t * 3] = x;
        vertices[t * 3 + 1] = y + r;
        vertices[t * 3 + 2] = z;

        normals[t * 3] = 0.0f;
        normals[t * 3 + 1] = 1.0f;
        normals[t * 3 + 2] = 0.0f;

        uvs[t * 2] = 0.0f;
        uvs[t * 2 + 1] = 0.0f;
        ++t;
    }

    // South pole of the sphere
    {
        vertices[t * 3] = x;
        vertices[t * 3 + 1] = y - r;
        vertices[t * 3 + 2] = z;

        normals[t * 3] = 0.0f;
        normals[t * 3 + 1] = -1.0f;
        normals[t * 3 + 2] = 0.0f;

        uvs[t * 2] = 1.0f;
        uvs[t * 2 + 1] = 1.0f;
        ++t;
    }

    t = 0U;
    for (std::uint32_t j = 0U; j < lat - 3; j++)
    {
        for (std::uint32_t i = 0U; i < lon - 1; i++)
        {
            indices[t++] = j * lon + i;
            indices[t++] = (j + 1) * lon + i + 1;
            indices[t++] = j * lon + i + 1;
            indices[t++] = j * lon + i;
            indices[t++] = (j + 1) * lon + i;
            indices[t++] = (j + 1) * lon + i + 1;
        }
    }

    for (std::uint32_t i = 0U; i < lon - 1; i++)
    {
        indices[t++] = (lat - 2) * lon;
        indices[t++] = i;
        indices[t++] = i + 1;
        indices[t++] = (lat - 2) * lon + 1;
        indices[t++] = (lat - 3) * lon + i + 1;
        indices[t++] = (lat - 3) * lon + i;
    }

    std::vector<rpr_int> faces(indices.size() / 3, 3);

    rpr_shape sphere = nullptr;
    rpr_int status = rprContextCreateMesh(context,
        vertices.data(), vertices.size(), sizeof(rpr_float) * 3,
        normals.data(), normals.size(), sizeof(rpr_float) * 3,
        uvs.data(), uvs.size(), sizeof(rpr_float) * 2,
        indices.data(), sizeof(rpr_int),
        indices.data(), sizeof(rpr_int),
        indices.data(), sizeof(rpr_int),
        faces.data(), faces.size(), &sphere);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(sphere, name);
    RETURN_IF_FAILED(status);

    status = rprShapeSetMaterial(sphere, material);
    RETURN_IF_FAILED(status);

    status = rprSceneAttachShape(scene, sphere);
    RETURN_IF_FAILED(status);

    if (out_sphere)
    {
        *out_sphere = sphere;
    }

    return RPR_SUCCESS;

}

rpr_int SceneIoTest::AddPlane(rpr_context context, rpr_scene scene, rpr_material_node material, rpr_char const* name,
    rpr_float x, rpr_float y, rpr_float z, rpr_float width, rpr_float height, rpr_float nx, rpr_float ny, rpr_float nz, rpr_shape* out_plane) const
{
    struct Vertex
    {
        rpr_float x, y, z;
        rpr_float nx, ny, nz;
        rpr_float u, v;
    };

    auto InvLength = [](rpr_float x, rpr_float y, rpr_float z)
    {
        return 1.0f / (x * x + y * y + z * z);
    };

    // Normalize plane normal
    rpr_float normal_inv_length = InvLength(nx, ny, nz);
    nx *= normal_inv_length;
    ny *= normal_inv_length;
    nz *= normal_inv_length;

    // Get an axiliary vector to build a tangent vector
    rpr_float axis_x = fabs(nx) > 0.001f ? 0.0f : 1.0f;
    rpr_float axis_y = 0.0f;
    rpr_float axis_z = fabs(nx) > 0.001f ? 1.0f : 0.0f;

    // Get plane tangent vector = cross(axis, normal)
    rpr_float tx = axis_y * nz - ny * axis_z;
    rpr_float ty = nx * axis_z - axis_x * nz;
    rpr_float tz = axis_x * ny - axis_y * nx;

    // Normalize tangent vector
    rpr_float t_inv_length = InvLength(tx, ty, tz);
    tx *= t_inv_length;
    ty *= t_inv_length;
    tz *= t_inv_length;

    // Get plane bitangent vector = cross(normal, tangent)
    rpr_float sx = ny * tz - ty * nz;
    rpr_float sy = tx * nz - nx * tz;
    rpr_float sz = nx * ty - ny * tx;

    Vertex vertices[4] =
    {
        {
            -sx * width - tx * height + x,
            -sy * width - ty * height + y,
            -sz * width - tz * height + z,
            nx, ny, nz, 0.0f, 0.0f
        },
        {
            sx * width - tx * height + x,
            sy * width - ty * height + y,
            sz * width - tz * height + z,
            nx, ny, nz, 1.0f, 0.0f
        },
        {
            sx * width + tx * height + x,
            sy * width + ty * height + y,
            sz * width + tz * height + z,
            nx, ny, nz, 1.0f, 1.0f
        },
        {
            -sx * width + tx * height + x,
            -sy * width + ty * height + y,
            -sz * width + tz * height + z,
            nx, ny, nz, 0.0f, 1.0f
        }
    };

    rpr_int indices[] =
    {
        3, 1, 0,
        2, 1, 3
    };

    rpr_int num_face_vertices[] =
    {
        3, 3
    };

    rpr_shape plane = nullptr;

    rpr_int status = rprContextCreateMesh(context,
        (rpr_float const*)&vertices[0], 4, sizeof(Vertex),
        (rpr_float const*)((char*)&vertices[0] + sizeof(rpr_float) * 3), 4, sizeof(Vertex),
        (rpr_float const*)((char*)&vertices[0] + sizeof(rpr_float) * 6), 4, sizeof(Vertex),
        indices, sizeof(rpr_int),
        indices, sizeof(rpr_int),
        indices, sizeof(rpr_int),
        num_face_vertices, 2, &plane);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(plane, name);
    RETURN_IF_FAILED(status);

    status = rprShapeSetMaterial(plane, material);
    RETURN_IF_FAILED(status);

    status = rprSceneAttachShape(scene, plane);
    RETURN_IF_FAILED(status);

    if (out_plane)
    {
        *out_plane = plane;
    }

    return RPR_SUCCESS;

}

rpr_int SceneIoTest::AddIbl(rpr_context context, rpr_scene scene, rpr_char const* image_filename, rpr_light* out_light) const
{
    rpr_light light = nullptr;
    rpr_int status = rprContextCreateEnvironmentLight(context, &light);
    RETURN_IF_FAILED(status);

    rpr_image image = nullptr;
    status = rprContextCreateImageFromFile(context, image_filename, &image);
    RETURN_IF_FAILED(status);

    status = rprEnvironmentLightSetImage(light, image);
    RETURN_IF_FAILED(status);

    status = rprSceneAttachLight(scene, light);
    RETURN_IF_FAILED(status);

    if (out_light)
    {
        *out_light = light;
    }

    return RPR_SUCCESS;

}

rpr_int SceneIoTest::AddSpotLight(rpr_context context, rpr_scene scene, rpr_float pitch, rpr_float yaw,
    rpr_float x, rpr_float y, rpr_float z, rpr_float r, rpr_float g, rpr_float b,
    rpr_float iangle, rpr_float oangle, rpr_light* out_light) const
{
    rpr_light light = nullptr;
    rpr_int status = rprContextCreateSpotLight(context, &light);
    RETURN_IF_FAILED(status);

    rpr_float const light_matrix[] =
    {
        std::cos(yaw), 0.0f, std::sin(yaw), x,
        std::sin(pitch) * std::sin(yaw), std::cos(pitch), -std::cos(yaw) * std::sin(pitch), y,
        -std::cos(pitch) * std::sin(yaw), std::sin(pitch), std::cos(pitch) * std::cos(yaw), z,
        0.0f, 0.0f, 0.0f, 1.0f
    };

    status = rprLightSetTransform(light, true, light_matrix);
    RETURN_IF_FAILED(status);

    status = rprSpotLightSetRadiantPower3f(light, r, g, b);
    RETURN_IF_FAILED(status);

    status = rprSpotLightSetConeShape(light, iangle, oangle);
    RETURN_IF_FAILED(status);

    status = rprSceneAttachLight(scene, light);
    RETURN_IF_FAILED(status);

    if (out_light)
    {
        *out_light = light;
    }

    return RPR_SUCCESS;

}

rpr_int SceneIoTest::LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const
{
    assert(out_scene);

    rpr_int status = RPR_SUCCESS;

    rpr_scene scene = nullptr;
    status = rprContextCreateScene(context, &scene);
    RETURN_IF_FAILED(status);

    std::string fname(filename);
    fname = fname.substr(0, fname.rfind(".test"));
    fname = fname.substr(fname.find(basepath) + strlen(basepath));

    rpr_char const* ibl_filename = "../Resources/Textures/studio015.hdr";

    if (fname == "plane+spot")
    {
        rpr_material_node plane_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "plane_material", 0.8f, 0.8f, 0.8f, &plane_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        status = AddSpotLight(context, scene, -(float)M_PI / 4.0f, 0.0f, 0.0f, 4.0f, 4.0f,
            100.0f, 50.0f, 50.0f, std::cos((float)M_PI / 6), std::cos((float)M_PI / 4), nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "plane+ibl")
    {
        rpr_material_node plane_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "plane_material", 0.8f, 0.8f, 0.8f, &plane_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        status = AddIbl(context, scene, ibl_filename, nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "sphere+ibl")
    {
        rpr_material_node sphere_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "sphere_material", 0.8f, 0.8f, 0.8f, &sphere_material);
        RETURN_IF_FAILED(status);

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 0.0f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        status = AddIbl(context, scene, ibl_filename, nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "sphere+plane")
    {
        rpr_material_node sphere_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "sphere_material", 0.8f, 0.8f, 0.8f, &sphere_material);
        RETURN_IF_FAILED(status);

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 2.5f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        rpr_material_node plane_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "plane_material", 0.8f, 0.8f, 0.8f, &plane_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "sphere+plane+area")
    {
        rpr_material_node sphere_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "sphere_material", 0.8f, 0.8f, 0.8f, &sphere_material);
        RETURN_IF_FAILED(status);

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 2.5f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        rpr_material_node plane_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "plane_material", 0.8f, 0.8f, 0.8f, &plane_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        rpr_material_node emissive_material = nullptr;
        status = CreateEmissiveMaterial(material_system, "emissive_material", 3.1f, 3.0f, 2.8f, &emissive_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, emissive_material, "arealight", 0.0f, 6.0f, 0.0f, 2.0f, 2.0f, 0.0f, -1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "sphere+plane+ibl")
    {
        rpr_material_node sphere_material = nullptr;
        {
            // Create refraction material for sphere
            rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &sphere_material);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_REFRACTION);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_COLOR, 0.7f, 1.0f, 0.7f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_IOR, 1.5f, 1.5f, 1.5f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_ROUGHNESS, 0.1f, 0.1f, 0.1f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprObjectSetName(sphere_material, "sphere_material");
            RETURN_IF_FAILED(status);

        }

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 2.5f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        rpr_material_node plane_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "plane_material", 0.8f, 0.8f, 0.8f, &plane_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        status = AddIbl(context, scene, ibl_filename, nullptr);
        RETURN_IF_FAILED(status);
    }
    else if (fname == "uberv2_test_spheres")
    {
        const rpr_float roughness = 0.05f;

        rpr_material_node sphere_material = nullptr;
        {
            rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &sphere_material);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, 1.0f, 1.0f, 1.0f, 1.0f);
            RETURN_IF_FAILED(status);
            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_COATING_COLOR,
                1.0f, 0.0f, 0.0f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
                roughness, roughness, roughness, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_COLOR,
                0.0f, 1.0f, 0.0f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_COLOR,
                0.0f, 0.0f, 1.0f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_ROUGHNESS,
                roughness, roughness, roughness, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_COATING |
                RPR_UBER_MATERIAL_LAYER_REFLECTION | RPR_UBER_MATERIAL_LAYER_REFRACTION);
            RETURN_IF_FAILED(status);

        }

        rpr_shape sphere = nullptr;
        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 0.0f, 0.0f, 0.9f, &sphere);
        RETURN_IF_FAILED(status);

        rpr_float const translation[] =
        {
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, -10.f,
            0.0f, 0.0f, 0.0f, 1.0f
        };

        status = rprShapeSetTransform(sphere, true, translation);
        RETURN_IF_FAILED(status);

        for (int i = 0; i < 5; ++i)
        {
            for (int j = 0; j < 10; ++j)
            {
                rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &sphere_material);
                RETURN_IF_FAILED(status);

                status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, 1.0f, 0.0f, 0.0f, 1.0f);
                RETURN_IF_FAILED(status);

                switch (i)
                {
                case 0:
                    status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                        RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_COATING);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_COATING_IOR,
                        1.0f + (float)j / 5.f, 1.0f + (float)j / 5.f, 1.0f + (float)j / 5.f, 1.0f);
                    RETURN_IF_FAILED(status);

                    break;
                case 1:
                    status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                        RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFLECTION);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
                        (float)j / 10.f, (float)j / 10.f, (float)j / 10.f, 1.0f);
                    RETURN_IF_FAILED(status);

                    break;
                case 2:
                    status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                        RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFLECTION);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
                        roughness, roughness, roughness, 1.0f);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_IOR,
                        1.0f + (float)j / 5.f, 1.0f + (float)j / 5.f, 1.0f + (float)j / 5.f, 1.0f);
                    RETURN_IF_FAILED(status);

                    break;
                case 3:
                    status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                        RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFRACTION);
                    RETURN_IF_FAILED(status);


                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFRACTION_ROUGHNESS,
                        (float)j / 10.f, (float)j / 10.f, (float)j / 10.f, 1.0f);
                    RETURN_IF_FAILED(status);

                    break;
                case 4:
                    status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFLECTION | RPR_UBER_MATERIAL_LAYER_TRANSPARENCY);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_TRANSPARENCY,
                        (float)j / 9.f, (float)j / 9.f, (float)j / 9.f, 1.0f);
                    RETURN_IF_FAILED(status);

                    status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
                        roughness, roughness, roughness, 1.0f);
                    RETURN_IF_FAILED(status);

                    break;
                }

                rpr_shape sphere_instance = nullptr;

                // TODO: Use instances here
                // We can't do it now because instances don't use their custom material, only material from the base mesh
                status = AddSphere(context, scene, sphere_material, "sphere_instance", 64, 32, 0.0f, 0.0f, 0.0f, 0.9f, &sphere_instance);
                RETURN_IF_FAILED(status);

                rpr_float const translation[] =
                {
                    1.0f, 0.0f, 0.0f, j * 2.f - 9.f,
                    0.0f, 1.0f, 0.0f, i * 2.f - 3.f,
                    0.0f, 0.0f, 1.0f, -10.f,
                    0.0f, 0.0f, 0.0f, 1.0f
                };

                status = rprShapeSetTransform(sphere_instance, true, translation);
                RETURN_IF_FAILED(status);

            }
        }

        AddIbl(context, scene, ibl_filename, nullptr);
    }
    else if (fname == "sphere+uberv2+ibl")
    {
        rpr_material_node texture_node = nullptr;
        rpr_int status = CreateTextureNode("test_albedo1.jpg", "../Resources/Textures",
            context, material_system, true, &texture_node);
        RETURN_IF_FAILED(status);

        rpr_material_node sphere_material = nullptr;
        {
            status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &sphere_material);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFLECTION);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(sphere_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_COLOR, texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS, 0.05f, 0.05f, 0.05f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprObjectSetName(sphere_material, "sphere_material");
            RETURN_IF_FAILED(status);

        }

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 0.0f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        status = AddIbl(context, scene, ibl_filename, nullptr);
        RETURN_IF_FAILED(status);

    }
    else if (fname == "sphere+plane_uberv2+ibl+normalmap")
    {
        rpr_material_node sphere_material = nullptr;
        {
            rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &sphere_material);
            RETURN_IF_FAILED(status);

            status = rprObjectSetName(sphere_material, "sphere_material");
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(sphere_material, RPR_UBER_MATERIAL_LAYERS,
                RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_REFLECTION);
            RETURN_IF_FAILED(status);

            rpr_material_node albedo_texture_node = nullptr;
            status = CreateTextureNode("test_albedo1.jpg", "../Resources/Textures", context, material_system, true, &albedo_texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(sphere_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, albedo_texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_COLOR, albedo_texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(sphere_material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS, 0.05f, 0.05f, 0.05f, 1.0f);
            RETURN_IF_FAILED(status);

        }

        status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 2.5f, 0.0f, 2.0f, nullptr);
        RETURN_IF_FAILED(status);

        rpr_material_node plane_material = nullptr;
        {
            rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &plane_material);
            RETURN_IF_FAILED(status);

            status = rprObjectSetName(plane_material, "plane_material");
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(plane_material, RPR_UBER_MATERIAL_LAYERS,
                RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_SHADING_NORMAL);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(plane_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, 0.8f, 0.8f, 0.8f, 1.0f);
            RETURN_IF_FAILED(status);

            rpr_material_node normalmap_texture_node = nullptr;
            status = CreateTextureNode("test_normal.jpg", "../Resources/Textures", context, material_system, false, &normalmap_texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(plane_material, RPR_UBER_MATERIAL_NORMAL, normalmap_texture_node);
            RETURN_IF_FAILED(status);

        }

        status = AddPlane(context, scene, plane_material, "plane", 0.0f, 0.0f, 0.0f, 8.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        AddIbl(context, scene, ibl_filename, nullptr);
    }
    else if (fname == "transparent_planes")
    {
        rpr_material_node transparent_material = nullptr;
        rpr_int status = RPR_SUCCESS;
        {
            // Create transparent material
            status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &transparent_material);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputU_ext(transparent_material, RPR_UBER_MATERIAL_LAYERS,
                RPR_UBER_MATERIAL_LAYER_DIFFUSE | RPR_UBER_MATERIAL_LAYER_TRANSPARENCY);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(transparent_material, RPR_UBER_MATERIAL_DIFFUSE_COLOR,
                1.0f, 0.0f, 0.0f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputF_ext(transparent_material, RPR_UBER_MATERIAL_TRANSPARENCY,
                0.9f, 0.9f, 0.9f, 1.0f);
            RETURN_IF_FAILED(status);

            status = rprObjectSetName(transparent_material, "transparent_material");
            RETURN_IF_FAILED(status);

        }

        for (int i = 0; i < 8; ++i)
        {
            status = AddPlane(context, scene, transparent_material, (std::string("plane") + std::to_string(i)).c_str(),
                i * 2.0f - 8.0f, 4.0f, 0.0f, 8.0f, 4.0f, 1.0f, 0.0f, 0.0f, nullptr);
            RETURN_IF_FAILED(status);

        }

        rpr_material_node floor_material = nullptr;
        status = CreateDiffuseMaterial(material_system, "floor_material", 0.8f, 0.8f, 0.8f, &floor_material);
        RETURN_IF_FAILED(status);

        status = AddPlane(context, scene, floor_material, "floor_plane", 0.0f, 0.0f, 0.0f, 16.0f, 8.0f, 0.0f, 1.0f, 0.0f, nullptr);
        RETURN_IF_FAILED(status);

        AddIbl(context, scene, ibl_filename, nullptr);
    }
    else if (fname == "400materials")
    {
        for (int i = 0; i < 400; ++i)
        {
            rpr_material_node sphere_material = nullptr;
            rpr_int status = CreateDiffuseMaterial(material_system, "sphere_material",
                std::cos((float)(i * M_PI) / 20.0f) * 0.5f + 0.5f,
                (float)i / 400.0f,
                std::sin((float)(i * M_PI) / 20.0f) * 0.5f + 0.5f,
                &sphere_material);
            RETURN_IF_FAILED(status);

            // TODO: Use instances here
            rpr_shape sphere = nullptr;
            status = AddSphere(context, scene, sphere_material, "sphere", 64, 32, 0.0f, 0.0f, 0.0f, 0.25f, &sphere);
            RETURN_IF_FAILED(status);

            rpr_float const translation[] =
            {
                1.0f, 0.0f, 0.0f, std::cos((float)(i * M_PI) / 20.0f) * 4.0f,
                0.0f, 1.0f, 0.0f, (float)i / 25.0f,
                0.0f, 0.0f, 1.0f, std::sin((float)(i * M_PI) / 20.0f) * 4.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            };

            status = rprShapeSetTransform(sphere, true, translation);
            RETURN_IF_FAILED(status);
        }

        AddIbl(context, scene, ibl_filename, nullptr);
    }
    else
    {
        // Test scene not found
        return RPR_ERROR_INVALID_PARAMETER;
    }

    status = rprContextSetScene(context, scene);
    RETURN_IF_FAILED(status);

    *out_scene = scene;

    return RPR_SUCCESS;

}

rpr_int SceneIoTest::SaveScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const
{
    return RPR_ERROR_UNSUPPORTED;
}
