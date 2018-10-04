#pragma once

#include "RadeonProRender.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>

class MaterialIo
{
public:
    // Create XML based material IO
    static std::unique_ptr<MaterialIo> CreateMaterialIoXML(rpr_context context, rpr_material_system material_system);

    using MaterialMap = std::map<std::string, std::string>;

    // Constructor
    MaterialIo() = default;
    // Destructor
    virtual ~MaterialIo() = 0;

    virtual rpr_int LoadMaterials(rpr_char const* filename, std::map<std::string, rpr_material_node> & materials) = 0;
    virtual rpr_int SaveMaterials(rpr_char const* filename, std::set<rpr_material_node> const& materials) = 0;
    virtual rpr_int SaveIdentityMapping(rpr_char const* filename, rpr_scene scene) = 0;
    virtual MaterialMap LoadMaterialMapping(rpr_char const* filename) = 0;

    rpr_int ReplaceSceneMaterials(rpr_scene scene, std::map<std::string, rpr_material_node> new_materials, MaterialMap const& mapping);

    rpr_int SaveMaterialsFromScene(rpr_char const* filename, rpr_scene scene);

    // Disallow copying
    MaterialIo(MaterialIo const&) = delete;
    MaterialIo& operator = (MaterialIo const&) = delete;

};

inline MaterialIo::~MaterialIo()
{
}
