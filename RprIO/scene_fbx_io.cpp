#define _USE_MATH_DEFINES
#define FBXSDK_NEW_API
#define KFBX_DLLINFO
#include <fbxsdk.h>

#include <cmath>

#include "scene_io.h"
#include <vector>
#include <set>
#include <cassert>
#include <stack>

using namespace fbxsdk;

class SceneIoFBX : public SceneIo::Loader
{
public:
    SceneIoFBX();
    ~SceneIoFBX();

    rpr_int LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const override;

    rpr_int SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const override;

private:
    rpr_int AddMesh(FbxNode* node, const char* basepath, rpr_context context, rpr_material_system material_system, rpr_scene scene) const;
    rpr_int AddLight(FbxNode* node, const char* basepath, rpr_context context, rpr_material_system material_system, rpr_scene scene) const;

    rpr_char const* FindTextureName(FbxProperty const& prop) const;

    rpr_int TranslateMaterial(const char* basepath, FbxSurfaceMaterial const* in_material,
        rpr_context context, rpr_material_system material_system, rpr_material_node* out_material) const;

    mutable std::map<FbxSurfaceMaterial const*, rpr_material_node> m_material_cache;
};

// Create an instance of fbx loader
static SceneIoFBX fbx_loader;

inline SceneIoFBX::SceneIoFBX()
    : SceneIo::Loader("fbx", this)
{
}

inline SceneIoFBX::~SceneIoFBX()
{
    SceneIo::UnregisterLoader("fbx");
}

rpr_char const* SceneIoFBX::FindTextureName(FbxProperty const& prop) const
{
    for (int i = 0; i < prop.GetSrcObjectCount<FbxFileTexture>(); i++)
    {
        FbxFileTexture* texture = prop.GetSrcObject<FbxFileTexture>(i);

        if (texture)
        {
            return texture->GetRelativeFileName();
        }
    }

    return nullptr;
}

rpr_int SceneIoFBX::TranslateMaterial(const char* basepath, FbxSurfaceMaterial const* in_material,
    rpr_context context, rpr_material_system material_system, rpr_material_node* out_material) const
{
    assert(out_material);

    auto iter = m_material_cache.find(in_material);

    if (iter != m_material_cache.cend())
    {
        *out_material = iter->second;
        return RPR_SUCCESS;
    }

    rpr_material_node material = nullptr;

    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &material);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(material, in_material->GetName());
    RETURN_IF_FAILED(status);

    auto SquareLength = [](FbxDouble const* vec)
    {
        return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    };

    FbxDouble3 emission = in_material->FindProperty(FbxSurfaceMaterial::sEmissive).Get<FbxDouble3>();

    // Special case: emissive material
    if (SquareLength(emission.mData) > 0.0)
    {
        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_EMISSION_COLOR,
            static_cast<rpr_float>(emission[0]),
            static_cast<rpr_float>(emission[1]),
            static_cast<rpr_float>(emission[2]),
            1.0f);
        RETURN_IF_FAILED(status);

        // Set material layers
        status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, RPR_UBER_MATERIAL_LAYER_EMISSION);
        RETURN_IF_FAILED(status);

        m_material_cache.emplace(std::make_pair(in_material, material));

        *out_material = material;

        return RPR_SUCCESS;

    }

    rpr_uint layers = 0;
    bool apply_gamma = true;

    constexpr rpr_float default_ior = 3.0f;
    rpr_float roughness = static_cast<rpr_float>(in_material->FindProperty(FbxSurfaceMaterial::sShininess).Get<FbxDouble>());
    // Convert FBX shininess to roughness
    roughness = 1.0f - std::log2f(roughness) / 10.0f;

    // Refraction layer
    FbxDouble transparency_factor = in_material->FindProperty(FbxSurfaceMaterial::sTransparencyFactor).Get<FbxDouble>();
    if (transparency_factor > 0.0)
    {
        layers |= RPR_UBER_MATERIAL_LAYER_REFRACTION;

        FbxDouble3 transparency_color = in_material->FindProperty(FbxSurfaceMaterial::sTransparentColor).Get<FbxDouble3>();
        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_COLOR,
            static_cast<rpr_float>(transparency_color[0]),
            static_cast<rpr_float>(transparency_color[1]),
            static_cast<rpr_float>(transparency_color[2]), 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_ROUGHNESS,
            roughness, roughness, roughness, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFRACTION_IOR,
            default_ior, default_ior, default_ior, 1.0f);
        RETURN_IF_FAILED(status);

    }

    // Reflection layer
    FbxDouble specular_factor = in_material->FindProperty(FbxSurfaceMaterial::sSpecularFactor).Get<FbxDouble>();
    if (specular_factor > 0.0)
    {
        layers |= RPR_UBER_MATERIAL_LAYER_REFLECTION;

        rpr_char const* specular_texname = FindTextureName(in_material->FindProperty(FbxSurfaceMaterial::sSpecular));

        if (specular_texname)
        {
            rpr_material_node texture_node = nullptr;
            status = CreateTextureNode(specular_texname, basepath, context, material_system, apply_gamma, &texture_node);
            RETURN_IF_FAILED(status);
            status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_REFLECTION_COLOR, texture_node);
            RETURN_IF_FAILED(status);
        }
        else
        {
            FbxDouble3 reflection_color = in_material->FindProperty(FbxSurfaceMaterial::sSpecular).Get<FbxDouble3>();
            status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_COLOR,
                static_cast<rpr_float>(reflection_color[0]),
                static_cast<rpr_float>(reflection_color[1]),
                static_cast<rpr_float>(reflection_color[2]), 1.0f);
            RETURN_IF_FAILED(status);
        }

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_ROUGHNESS,
            roughness, roughness, roughness, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_IOR,
            default_ior, default_ior, default_ior, 1.0f);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_REFLECTION_METALNESS,
            1.0f, 1.0f, 1.0f, 1.0f);
        RETURN_IF_FAILED(status);

    }

    rpr_char const* bump_texname = FindTextureName(in_material->FindProperty(FbxSurfaceMaterial::sBump));

    // Set a bump map if we have one
    if (bump_texname)
    {
        layers |= RPR_UBER_MATERIAL_LAYER_SHADING_NORMAL;

        rpr_material_node bump_texture_node = nullptr;
        status = CreateTextureNode(bump_texname, basepath, context, material_system, apply_gamma, &bump_texture_node);
        RETURN_IF_FAILED(status);
        status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_NORMAL, bump_texture_node);
        RETURN_IF_FAILED(status);
    }

    // Diffuse layer
    {
        layers |= RPR_UBER_MATERIAL_LAYER_DIFFUSE;

        rpr_char const* diffuse_texname = FindTextureName(in_material->FindProperty(FbxSurfaceMaterial::sDiffuse));

        if (diffuse_texname)
        {
            rpr_material_node texture_node = nullptr;
            status = CreateTextureNode(diffuse_texname, basepath, context, material_system, apply_gamma, &texture_node);
            RETURN_IF_FAILED(status);
            status = rprMaterialNodeSetInputN_ext(material, RPR_UBER_MATERIAL_DIFFUSE_COLOR, texture_node);
            RETURN_IF_FAILED(status);
        }
        else
        {
            FbxDouble3 diffuse_color = in_material->FindProperty(FbxSurfaceMaterial::sDiffuse).Get<FbxDouble3>();
            status = rprMaterialNodeSetInputF_ext(material, RPR_UBER_MATERIAL_DIFFUSE_COLOR,
                static_cast<rpr_float>(diffuse_color[0]),
                static_cast<rpr_float>(diffuse_color[1]),
                static_cast<rpr_float>(diffuse_color[2]), 1.0f);
            RETURN_IF_FAILED(status);
        }

    }

    // Set material layers
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, layers);
    RETURN_IF_FAILED(status);

    // Set material name
    status = rprObjectSetName(material, in_material->GetName());
    RETURN_IF_FAILED(status);

    m_material_cache.emplace(std::make_pair(in_material, material));

    *out_material = material;

    return RPR_SUCCESS;
}


rpr_int SceneIoFBX::AddMesh(FbxNode* node, const char* basepath, rpr_context context, rpr_material_system material_system, rpr_scene scene) const
{
    FbxMesh const* fbx_mesh = node->GetMesh();
    FbxMatrix const& transform = node->EvaluateGlobalTransform();
    FbxAMatrix invtransp = node->EvaluateGlobalTransform().Inverse().Transpose();

    {
        auto num_triangles = fbx_mesh->GetPolygonCount();
        auto vertex_ptr = fbx_mesh->GetControlPoints();

        // 3 positions per triangle, each have x, y, z component
        std::vector<rpr_float> positions(num_triangles * 3 * 3);
        // 3 normals per triangle, each have x, y, z component
        std::vector<rpr_float> normals(num_triangles * 3 * 3);
        // 3 texcoords per triangle, each have u, v component
        std::vector<rpr_float> texcoords(num_triangles * 3 * 2);

        std::vector<rpr_int> indices(num_triangles * 3);

        FbxStringList uv_list;
        fbx_mesh->GetUVSetNames(uv_list);

        for (int tri = 0; tri < num_triangles; ++tri)
        {
            assert(fbx_mesh->GetPolygonSize(tri) == 3);

            for (int v = 0; v < fbx_mesh->GetPolygonSize(tri); ++v)
            {
                // Get position of the vertex
                int vertex_index = fbx_mesh->GetPolygonVertex(tri, v);
                int current_index = tri * 3 + v;

                FbxVector4 vertex = vertex_ptr[vertex_index];
                vertex = transform.MultNormalize(vertex);
                for (int k = 0; k < 3; ++k)
                {
                    positions[current_index * 3 + k] = static_cast<rpr_float>(vertex[k]);
                }

                // Get normal of the vertex
                FbxVector4 n;
                fbx_mesh->GetPolygonVertexNormal(tri, v, n);
                n = invtransp.MultT(n);
                n.Normalize();

                for (int k = 0; k < 3; ++k)
                {
                    normals[current_index * 3 + k] = static_cast<rpr_float>(n[k]);
                }

                // Get texture coordinate of the vertex
                if (uv_list.GetCount() > 0)
                {
                    FbxVector2 uv;
                    bool unmapped;
                    fbx_mesh->GetPolygonVertexUV(tri, v, uv_list.GetStringAt(0), uv, unmapped);

                    for (int k = 0; k < 2; ++k)
                    {
                        texcoords[current_index * 2 + k] = unmapped ? 0.0f : static_cast<rpr_float>(uv[k]);
                    }
                }

                // Store indices of processed vertex
                indices[current_index] = current_index;

            }
        }

        // Each face is a triangle
        std::vector<rpr_int> faces(indices.size() / 3, 3);

        // Create mesh
        rpr_shape mesh = nullptr;
        rpr_int status = rprContextCreateMesh(context,
            positions.data(), positions.size(), sizeof(rpr_float) * 3,
            normals.data(),   normals.size(),   sizeof(rpr_float) * 3,
            texcoords.data(), texcoords.size(), sizeof(rpr_float) * 2,
            indices.data(),   sizeof(rpr_int),
            indices.data(),   sizeof(rpr_int),
            indices.data(),   sizeof(rpr_int),
            faces.data(), faces.size(), &mesh);
        RETURN_IF_FAILED(status);

        status = rprObjectSetName(mesh, node->GetName());
        RETURN_IF_FAILED(status);

        FbxLayerElementArrayTemplate<int>* material_indices = nullptr;
        fbx_mesh->GetMaterialIndices(&material_indices);

        if (material_indices)
        {
            FbxSurfaceMaterial const* fbx_material = node->GetMaterial(material_indices->GetAt(0));
            rpr_material_node material = nullptr;
            status = TranslateMaterial(basepath, fbx_material, context, material_system, &material);
            RETURN_IF_FAILED(status);

            status = rprShapeSetMaterial(mesh, material);
            RETURN_IF_FAILED(status);
        }

        // Attach to the scene
        status = rprSceneAttachShape(scene, mesh);
        RETURN_IF_FAILED(status);
    }

    return RPR_SUCCESS;
}

rpr_int SceneIoFBX::AddLight(FbxNode* node, const char* basepath, rpr_context context, rpr_material_system material_system, rpr_scene scene) const
{
    FbxLight const* fbx_light = node->GetLight();

    FbxDouble intensity = fbx_light->Intensity.Get();
    FbxDouble3 color = fbx_light->Color.Get();

    rpr_light light = nullptr;
    rpr_int status = RPR_SUCCESS;

    switch (fbx_light->LightType.Get())
    {
    case FbxLight::ePoint:
    {
        status = rprContextCreatePointLight(context, &light);
        RETURN_IF_FAILED(status);

        status = rprPointLightSetRadiantPower3f(light,
            static_cast<rpr_float>(intensity * color[0]),
            static_cast<rpr_float>(intensity * color[1]),
            static_cast<rpr_float>(intensity * color[2]));
        RETURN_IF_FAILED(status);

        break;
    }
    case FbxLight::eDirectional:
    {
        status = rprContextCreateDirectionalLight(context, &light);
        RETURN_IF_FAILED(status);

        status = rprDirectionalLightSetRadiantPower3f(light,
            static_cast<rpr_float>(intensity * color[0]),
            static_cast<rpr_float>(intensity * color[1]),
            static_cast<rpr_float>(intensity * color[2]));
        RETURN_IF_FAILED(status);
        break;
    }
    case FbxLight::eSpot:
    {
        status = rprContextCreateSpotLight(context, &light);
        RETURN_IF_FAILED(status);

        status = rprSpotLightSetRadiantPower3f(light,
            static_cast<rpr_float>(intensity * color[0]),
            static_cast<rpr_float>(intensity * color[1]),
            static_cast<rpr_float>(intensity * color[2]));
        RETURN_IF_FAILED(status);

        FbxDouble inner_angle = fbx_light->InnerAngle.Get();
        FbxDouble outer_angle = fbx_light->OuterAngle.Get();
        status = rprSpotLightSetConeShape(light,
            static_cast<rpr_float>(std::cos(inner_angle / 180.0 * M_PI)),
            static_cast<rpr_float>(std::cos(outer_angle / 180.0 * M_PI)));
        RETURN_IF_FAILED(status);
        break;
    }
    default:
        return RPR_ERROR_IO_ERROR;
    }

    // Get light transform matrix
    FbxAMatrix const& transform = node->EvaluateGlobalTransform();
    rpr_float light_transform[16];
    for (int i = 0; i < 16; ++i)
    {
        light_transform[i] = static_cast<rpr_float>(transform.Get(i / 4, i % 4));
    }

    // Swap Y and Z rows in rotation matrix
    std::swap_ranges(light_transform + 4, light_transform + 7, light_transform + 8);

    status = rprLightSetTransform(light, false, light_transform);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(light, node->GetName());
    RETURN_IF_FAILED(status);

    status = rprSceneAttachLight(scene, light);
    RETURN_IF_FAILED(status);

    return RPR_SUCCESS;
}

rpr_int SceneIoFBX::LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene* out_scene) const
{
    assert(out_scene);

    FbxManager* fbx_manager = FbxManager::Create();
    FbxImporter* fbx_importer = FbxImporter::Create(fbx_manager, "Baikal FBX Importer");
    FbxScene* fbx_scene = FbxScene::Create(fbx_manager, "Scene");

    if (!fbx_importer->Initialize(filename))
    {
        return RPR_ERROR_IO_ERROR;
    }
    if (!fbx_importer->Import(fbx_scene))
    {
        return RPR_ERROR_IO_ERROR;
    }

    auto fbx_root_node = fbx_scene->GetRootNode();
    assert(fbx_root_node);

    FbxGeometryConverter converter(fbx_manager);
    converter.Triangulate(fbx_scene, true);

    std::stack<FbxNode*> node_stack;

    rpr_scene scene = nullptr;
    rpr_int status = rprContextCreateScene(context, &scene);
    RETURN_IF_FAILED(status);

    node_stack.push(fbx_root_node);
    while (!node_stack.empty())
    {
        auto node = node_stack.top();
        node_stack.pop();

        for (auto c = 0; c < node->GetChildCount(); ++c)
        {
            node_stack.push(node->GetChild(c));
        }

        auto attribs = node->GetNodeAttribute();

        if (!attribs)
        {
            continue;
        }

        switch (attribs->GetAttributeType())
        {
        case FbxNodeAttribute::eMesh:
            status = AddMesh(node, basepath, context, material_system, scene);
            break;
        case FbxNodeAttribute::eLight:
            status = AddLight(node, basepath, context, material_system, scene);
            break;
        default:
            // Skip other attributes
            break;
        }

        RETURN_IF_FAILED(status);
    }

    fbx_scene->Destroy();
    fbx_importer->Destroy();
    fbx_manager->Destroy();

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
        status = rprDirectionalLightSetRadiantPower3f(directional_light, 10.0f, 10.0f, 10.0f);
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

rpr_int SceneIoFBX::SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system material_system, rprx_context /* uber_context */, rpr_scene scene) const
{
    return RPR_ERROR_UNIMPLEMENTED;
}
