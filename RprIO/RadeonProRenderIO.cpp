#include "RadeonProRenderIO.h"
#include "scene_io.h"
#include "material_io.h"

rpr_int rprLoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* scene)
{
    return SceneIo::LoadScene(filename, basepath, context, materialSystem, uberMatContext, scene);
}

rpr_int rprSaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context,
    rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene)
{
    return SceneIo::SaveScene(filename, basepath, context, materialSystem, uberMatContext, scene);
}

rpr_int rprReplaceSceneMaterials(rpr_char const* materials_xml, rpr_char const* mapping_xml,
    rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rpr_scene scene)
{
    auto material_io = MaterialIo::CreateMaterialIoXML();
    MaterialIo::MaterialMap material_mapping = material_io->LoadMaterialMapping(mapping_xml);

    std::map<std::string, rpr_material_node> new_materials;
    rpr_int status = material_io->LoadMaterials(materials_xml, materialSystem, new_materials);
    RETURN_IF_FAILED(status);

    return material_io->ReplaceSceneMaterials(scene, new_materials, material_mapping);
}

rpr_int rprSaveSceneMaterials(rpr_char const* materials_xml, rpr_char const* mapping_xml,
    rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rpr_scene scene)
{
    auto material_io = MaterialIo::CreateMaterialIoXML();
    rpr_int status = material_io->SaveIdentityMapping(mapping_xml, scene);
    RETURN_IF_FAILED(status);

    return material_io->SaveMaterialsFromScene(materials_xml, scene);

}
