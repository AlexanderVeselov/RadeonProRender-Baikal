#include "scene_io.h"
#include "tiny_obj_loader.h"
#include <vector>
#include <set>
#include <cassert>

class SceneIoOBJ : public SceneIo::Loader
{
public:
    SceneIoOBJ();
    ~SceneIoOBJ();

    rpr_int LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const override;

    rpr_int SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const override;

private:
    rpr_int TranslateMaterial(const char* basepath, tinyobj::material_t in_material,
        rpr_context context, rpr_material_system material_system, rpr_material_node* out_material) const;

    mutable std::map<std::string, rpr_material_node> m_material_cache;
};

// Create an instance of gltf loader
static SceneIoOBJ obj_loader;

inline SceneIoOBJ::SceneIoOBJ()
    : SceneIo::Loader("obj", this)
{
    SceneIo::RegisterLoader("objm", this);
}

inline SceneIoOBJ::~SceneIoOBJ()
{
    SceneIo::UnregisterLoader("objm");
}

rpr_int SceneIoOBJ::TranslateMaterial(const char* basepath, tinyobj::material_t in_material,
    rpr_context context, rpr_material_system material_system, rpr_material_node* out_material) const
{
    assert(out_material);

    auto iter = m_material_cache.find(in_material.name);

    if (iter != m_material_cache.cend())
    {
        *out_material = iter->second;
        return RPR_SUCCESS;
    }

    rpr_int status = RPR_SUCCESS;
    rpr_material_node material = nullptr;

    status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &material);
    RETURN_IF_FAILED(status);

    auto SquareLength = [](tinyobj::real_t const* vec)
    {
        return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    };

    // Special case: emissive material
    if (SquareLength(in_material.emission) > 0.0f)
    {
        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_EMISSION_COLOR,
            in_material.emission[0], in_material.emission[1], in_material.emission[2], 1.0f);
        RETURN_IF_FAILED(status);

        // Set material layers
        status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_EMISSION);
        RETURN_IF_FAILED(status);

        // Set material name
        status = rprObjectSetName(material, in_material.name.c_str());
        RETURN_IF_FAILED(status);

        m_material_cache.emplace(std::make_pair(in_material.name, material));

        *out_material = material;

        return RPR_SUCCESS;

    }

    rpr_uint layers = 0;
    bool apply_gamma = true;

    constexpr rpr_float default_ior = 3.0f;
    constexpr rpr_float default_roughness = 0.01f;

    // Refraction layer
    if (SquareLength(in_material.transmittance) > 0.0f)
    {
        layers |= RPR_UBER_MATERIAL_LAYER_REFRACTION;

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_COLOR,
            in_material.transmittance[0], in_material.transmittance[1], in_material.transmittance[2], 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_ROUGHNESS,
            default_roughness, default_roughness, default_roughness, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_IOR,
            default_ior, default_ior, default_ior, 1.0f);
        RETURN_IF_FAILED(status);

    }

    // Reflection layer
    if (SquareLength(in_material.specular) > 0.0f)
    {
        layers |= RPR_UBER_MATERIAL_LAYER_REFLECTION;

        if (!in_material.reflection_texname.empty())
        {
            rpr_material_node texture_node = nullptr;
            status = CreateTextureNode(in_material.reflection_texname.c_str(), basepath,
                context, material_system, apply_gamma, &texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_REFLECTION_COLOR, texture_node);
            RETURN_IF_FAILED(status);
        }
        else
        {
            status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_COLOR,
                in_material.specular[0], in_material.specular[1], in_material.specular[2], 1.0f);
            RETURN_IF_FAILED(status);
        }

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
            default_roughness, default_roughness, default_roughness, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_IOR,
            default_ior, default_ior, default_ior, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_METALNESS,
            1.0f, 1.0f, 1.0f, 1.0f);
        RETURN_IF_FAILED(status);

    }

    // Set a bump map if we have one
    if (!in_material.bump_texname.empty())
    {
        layers |= RPR_UBER_MATERIAL_LAYER_SHADING_NORMAL;

        rpr_material_node texture_node = nullptr;
        status = CreateTextureNode(in_material.bump_texname.c_str(), basepath,
            context, material_system, apply_gamma, &texture_node);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_NORMAL, texture_node);
        RETURN_IF_FAILED(status);
    }

    // Diffuse layer
    {
        layers |= RPR_UBER_MATERIAL_LAYER_DIFFUSE;

        if (!in_material.diffuse_texname.empty())
        {
            rpr_material_node texture_node = nullptr;
            status = CreateTextureNode(in_material.diffuse_texname.c_str(), basepath,
                context, material_system, apply_gamma, &texture_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, texture_node);
            RETURN_IF_FAILED(status);
        }
        else
        {
            status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_DIFFUSE_COLOR,
                in_material.diffuse[0], in_material.diffuse[1], in_material.diffuse[2], 1.0f);
            RETURN_IF_FAILED(status);
        }

    }

    // Set material layers
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, layers);
    RETURN_IF_FAILED(status);

    // Set material name
    status = rprObjectSetName(material, in_material.name.c_str());
    RETURN_IF_FAILED(status);

    m_material_cache.emplace(std::make_pair(in_material.name, material));

    *out_material = material;

    return RPR_SUCCESS;
}

rpr_int SceneIoOBJ::LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const
{
    assert(out_scene);

    using namespace tinyobj;
    // Loader data
    std::vector<shape_t> objshapes;
    std::vector<material_t> objmaterials;
    attrib_t attrib;

    // Try loading file
    std::string err;
    bool res = LoadObj(&attrib, &objshapes, &objmaterials, &err, filename, basepath, true);
    if (!res)
    {
        return RPR_ERROR_IO_ERROR;
    }

    rpr_int status = RPR_SUCCESS;

    rpr_scene scene = nullptr;
    status = rprContextCreateScene(context, &scene);
    RETURN_IF_FAILED(status);

    std::vector<rpr_material_node> rpr_materials(objmaterials.size());
    for (std::size_t i = 0; i < objmaterials.size(); ++i)
    {
        rpr_material_node material = nullptr;
        status = TranslateMaterial(basepath, objmaterials[i], context, material_system, &material);
        RETURN_IF_FAILED(status);

        rpr_materials[i] = material;
    }

    // Enumerate all shapes in the scene
    for (int s = 0; s < (int)objshapes.size(); ++s)
    {
        shape_t const& shape = objshapes[s];

        // Find all materials used by this shape.
        std::set<int> used_materials(std::begin(shape.mesh.material_ids), std::end(shape.mesh.material_ids));

        // Split the mesh into multiple meshes, each with only one material.
        for (int used_material : used_materials)
        {
            // Map from old index to new index.
            auto index_comp = [](index_t const& a, index_t const& b)
            {
                return (a.vertex_index < b.vertex_index)
                    || (a.vertex_index == b.vertex_index && a.normal_index < b.normal_index)
                    || (a.vertex_index == b.vertex_index && a.normal_index == b.normal_index
                        && a.texcoord_index < b.texcoord_index);
            };
            std::map<index_t, rpr_int, decltype(index_comp)> used_indices(index_comp);

            // Remapped indices.
            std::vector<rpr_int> indices;

            // Collected vertex/normal/texcoord data.
            std::vector<rpr_float> vertices, normals, texcoords;

            // Go through each face in the mesh.
            for (std::size_t i = 0; i < shape.mesh.material_ids.size(); ++i)
            {
                // Skip faces which don't use the current material.
                if (shape.mesh.material_ids[i] != used_material)
                {
                    continue;
                }

                const int num_face_vertices = shape.mesh.num_face_vertices[i];
                assert(num_face_vertices == 3 && "expected triangles");

                // For each vertex index of this face.
                for (int j = 0; j < num_face_vertices; ++j)
                {
                    index_t old_index = shape.mesh.indices[num_face_vertices * i + j];

                    // Collect vertex/normal/texcoord data. Avoid inserting the same data twice.
                    auto result = used_indices.emplace(old_index, (unsigned int)(vertices.size() / 3));
                    if (result.second) // Did insert?
                    {
                        // Push the new data.
                        for (int k = 0; k < 3; ++k)
                        {
                            vertices.push_back(attrib.vertices[3 * old_index.vertex_index + k]);
                        }

                        for (int k = 0; k < 3; ++k)
                        {
                            normals.push_back(attrib.normals[3 * old_index.normal_index + k]);
                        }

                        for (int k = 0; k < 2; ++k)
                        {
                            // If an uv is present
                            if (old_index.texcoord_index != -1)
                            {
                                texcoords.push_back(attrib.texcoords[2 * old_index.texcoord_index + k]);
                            }
                            else
                            {
                                texcoords.push_back(0.0f);
                            }
                        }

                    }

                    const unsigned int new_index = result.first->second;
                    indices.push_back(new_index);
                }
            }

            // Each face is a triangle
            std::vector<rpr_int> faces(indices.size() / 3, 3);

            // Create mesh
            rpr_shape mesh = nullptr;
            status = rprContextCreateMesh(context,
                vertices.data(),  vertices.size(),  sizeof(rpr_float) * 3,
                normals.data(),   normals.size(),   sizeof(rpr_float) * 3,
                texcoords.data(), texcoords.size(), sizeof(rpr_float) * 2,
                indices.data(),   sizeof(rpr_int),
                indices.data(),   sizeof(rpr_int),
                indices.data(),   sizeof(rpr_int),
                faces.data(), faces.size(), &mesh);
            RETURN_IF_FAILED(status);

            // Set material
            if (used_material >= 0)
            {
                status = rprShapeSetMaterial(mesh, rpr_materials[used_material]);
                RETURN_IF_FAILED(status);
            }

            // Attach to the scene
            status = rprSceneAttachShape(scene, mesh);
            RETURN_IF_FAILED(status);

        }
    }

    // Add directional and environment lights
    rpr_light directional_light = nullptr;
    {
        status = rprContextCreateDirectionalLight(context, &directional_light);
        RETURN_IF_FAILED(status);

        rpr_float lightm[] =
        {
             0.5f,  0.0f,  0.87f, 0.0f,
             0.75f, 0.5f, -0.43f, 0.0f,
            -0.43f, 0.87f, 0.25f, 0.0f,
             0.0f,  0.0f,  0.0f,  1.0f
        };

        status = rprLightSetTransform(directional_light, false, lightm);
        RETURN_IF_FAILED(status);
        status = rprDirectionalLightSetRadiantPower3f(directional_light, 2.0f, 1.9f, 1.7f);
        RETURN_IF_FAILED(status);
        status = rprSceneAttachLight(scene, directional_light);
        RETURN_IF_FAILED(status);
    }

    rpr_light env_light = nullptr;
    {
        status = rprContextCreateEnvironmentLight(context, &env_light);
        RETURN_IF_FAILED(status);
        rpr_image image = nullptr;
        status = rprContextCreateImageFromFile(context, "../Resources/Textures/studio015.hdr", &image);
        RETURN_IF_FAILED(status);
        status = rprEnvironmentLightSetImage(env_light, image);
        RETURN_IF_FAILED(status);
        status = rprSceneAttachLight(scene, env_light);
        RETURN_IF_FAILED(status);
    }

    status = rprContextSetScene(context, scene);
    RETURN_IF_FAILED(status);

    *out_scene = scene;

    return RPR_SUCCESS;
}

rpr_int SceneIoOBJ::SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const
{
    return RPR_ERROR_UNIMPLEMENTED;
}
