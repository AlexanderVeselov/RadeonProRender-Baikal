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
{/*
    std::ifstream in_materials(basepath + "materials.xml");
    std::ifstream in_mapping(basepath + "mapping.xml");*/
    return RPR_ERROR_UNIMPLEMENTED;
}

rpr_int rprSaveSceneMaterials(rpr_char const* materials_xml, rpr_char const* mapping_xml,
    rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rpr_scene scene)
{
    auto material_io = MaterialIo::CreateMaterialIoXML();
    return material_io->SaveMaterialsFromScene("materials.xml", scene);

}
