#include "ProRenderGLTF.h"
#include "scene_io.h"

class SceneIoGLTF : public SceneIo::Loader
{
public:
    SceneIoGLTF();
    ~SceneIoGLTF() = default;

    rpr_int LoadScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
        rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* out_scene) const override;

    rpr_int SaveScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
        rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene) const override;
};

// Create an instance of gltf loader
static SceneIoGLTF gltf_loader;

inline SceneIoGLTF::SceneIoGLTF()
    : SceneIo::Loader("gltf", this)
{}

rpr_int SceneIoGLTF::LoadScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* out_scene) const
{
    return rprImportFromGLTF(filename, context, materialSystem, uberMatContext, out_scene);
}

rpr_int SceneIoGLTF::SaveScene(rpr_char const* filename, rpr_char const* /* basepath */, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene) const
{
    return rprExportToGLTF(filename, context, materialSystem, uberMatContext, &scene, 1);
}
