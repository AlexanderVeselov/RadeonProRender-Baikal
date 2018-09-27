#pragma once

#include "RadeonProRenderIO.h"
#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <experimental/filesystem>

class SceneIo
{
public:
    class Loader
    {
    public:
        // Load the scene from file using resourse base path
        virtual rpr_int LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
            rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* out_scene) const = 0;

        virtual rpr_int SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
            rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene) const = 0;

        Loader(rpr_char const* ext, SceneIo::Loader *loader);
        virtual ~Loader();

        rpr_int LoadImage(rpr_char const* filename, rpr_char const* basepath,
            rpr_context context, rpr_image* out_texture) const;

        rpr_int CreateTextureNode(rpr_char const* filename, rpr_char const* basepath,
            rpr_context context, rpr_material_system material_system, bool apply_gamma, rpr_material_node* out_node) const;

    private:
        Loader(const Loader &) = delete;
        Loader& operator= (const Loader &) = delete;


        rpr_char const* m_ext;
        mutable std::map<std::string, rpr_image> m_image_cache;
    };

    // Registers extension handler
    static void RegisterLoader(rpr_char const* ext, SceneIo::Loader *loader);
    // Deregisters extension handler
    static void UnregisterLoader(rpr_char const* ext);

    static rpr_int LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* scene);

    static rpr_int SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
        rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene);

private:
    static SceneIo* GetInstance();
    static SceneIo::Loader* FindLoader(std::string const& filename);

    // Constructor
    SceneIo() = default;
    // Destructor
    virtual ~SceneIo() = default;

    // Disallow copying
    SceneIo(SceneIo const&) = delete;
    SceneIo& operator = (SceneIo const&) = delete;

    std::map<std::string, SceneIo::Loader*> m_loaders;
};

#define RETURN_IF_FAILED(status) \
    if ((status) != RPR_SUCCESS) \
    {                            \
        return (status);         \
    }

inline void SceneIo::RegisterLoader(rpr_char const* ext, SceneIo::Loader *loader)
{
    GetInstance()->m_loaders[ext] = loader;
}

inline void SceneIo::UnregisterLoader(rpr_char const* ext)
{
    GetInstance()->m_loaders.erase(ext);
}

inline SceneIo* SceneIo::GetInstance()
{
    static SceneIo instance;
    return &instance;
}

inline SceneIo::Loader* SceneIo::FindLoader(std::string const& filename)
{
    std::string ext = filename.substr(filename.rfind(".") + 1);

    SceneIo *instance = GetInstance();
    auto loader_it = instance->m_loaders.find(ext);
    if (loader_it == instance->m_loaders.end())
    {
        return nullptr;
    }

    return loader_it->second;
}

inline rpr_int SceneIo::LoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* scene)
{
    Loader* loader = FindLoader(filename);
    if (!loader)
    {
        return RPR_ERROR_IO_ERROR;
    }

    return loader->LoadScene(filename, basepath, context, materialSystem, uberMatContext, scene);
}

inline rpr_int SceneIo::SaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene)
{
    Loader* loader = FindLoader(filename);
    if (!loader)
    {
        return RPR_ERROR_IO_ERROR;
    }

    return loader->SaveScene(filename, basepath, context, materialSystem, uberMatContext, scene);
}

inline SceneIo::Loader::Loader(rpr_char const* ext, SceneIo::Loader *loader)
    : m_ext(ext)
{
    SceneIo::RegisterLoader(m_ext, loader);
}

inline SceneIo::Loader::~Loader()
{
    SceneIo::UnregisterLoader(m_ext);
}

inline rpr_int SceneIo::Loader::LoadImage(rpr_char const* filename,
    rpr_char const* basepath, rpr_context context, rpr_image* out_image) const
{
    assert(out_image);

    std::string fname;
    // Do not use basepath if filename is absolute
    if (std::experimental::filesystem::path(filename).is_relative())
    {
        fname.append(basepath);
    }

    fname.append(filename);
    std::replace(fname.begin(), fname.end(), '\\', '/');

    auto iter = m_image_cache.find(fname);

    if (iter != m_image_cache.cend())
    {
        *out_image = iter->second;
        return RPR_SUCCESS;
    }

    rpr_image image = nullptr;
    rpr_int status = rprContextCreateImageFromFile(context, fname.c_str(), &image);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(image, filename);
    RETURN_IF_FAILED(status);

    m_image_cache.emplace(fname, image);
    *out_image = image;

    return RPR_SUCCESS;

}

inline rpr_int SceneIo::Loader::CreateTextureNode(rpr_char const* filename, rpr_char const* basepath,
    rpr_context context, rpr_material_system material_system, bool apply_gamma, rpr_material_node* out_node) const
{
    assert(out_node);

    rpr_image image = nullptr;
    rpr_int status = LoadImage(filename, basepath, context, &image);
    RETURN_IF_FAILED(status);

    rpr_material_node texture_node = nullptr;
    status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_IMAGE_TEXTURE, &texture_node);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputImageData(texture_node, "data", image);
    RETURN_IF_FAILED(status);

    if (apply_gamma)
    {
        rpr_material_node gamma_power_node = nullptr;
        status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_ARITHMETIC, &gamma_power_node);
        RETURN_IF_FAILED(status);
        status = rprMaterialNodeSetInputU(gamma_power_node, "op", RPR_MATERIAL_NODE_OP_POW);
        RETURN_IF_FAILED(status);
        status = rprMaterialNodeSetInputN(gamma_power_node, "color0", texture_node);
        RETURN_IF_FAILED(status);
        static constexpr float gamma_power = 1.0f / 2.2f;
        status = rprMaterialNodeSetInputF(gamma_power_node, "color1", gamma_power, gamma_power, gamma_power, gamma_power);
        RETURN_IF_FAILED(status);

        *out_node = gamma_power_node;

        return RPR_SUCCESS;
    }

    *out_node = texture_node;

    return RPR_SUCCESS;

}
