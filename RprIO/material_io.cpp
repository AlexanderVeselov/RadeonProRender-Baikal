#include "XML/tinyxml2.h"
#include "material_io.h"

#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <stack>
#include <string>
#include <cassert>


using namespace tinyxml2;

//#define RETURN_IF_FAILED(status) \
//    if ((status) != RPR_SUCCESS) \
//    {                            \
//        return (status);         \
//    }

void RETURN_IF_FAILED(rpr_int status)
{
    if (status != RPR_SUCCESS)
    printf("%d\n", status);

}

namespace Baikal
{
    // This enum should be identical to _enum class Baikal::InputMap::InputMapType_
    // declared in BaikalNext/SceneGraph/inputmaps.h
    enum InputMapType
    {
        kConstantFloat3 = 0, // Holds constant float3 value. 
        kConstantFloat, // Holds constant float value
        kSampler, // Samples value from provided texture
        kAdd, // a + b
        kSub, // a - b
        kMul, // a * b
        kDiv, // a / b
        kSin, // sin(arg)
        kCos, // cos(arg)
        kTan, // tan(arg)
        kSelect, // Selects single component from input arg
        kDot3, // dot(a.xyz, b.xyz)
        kCross3, // cross(a.xyz, b.xyz)
        kLength3, // length(a.xyz)
        kNormalize3, // normalize(a.xyz)
        kPow, // pow(a, b)
        kAcos, // acos(arg)
        kAsin, // asin(arg)
        kAtan, // atan(arg)
        kLerp, // mix (a, b, control)
        kMin, // min(a, b)
        kMax, // max(a, b)
        kFloor, // floor(arg)
        kMod, // fmod(a, b)
        kAbs, // fabs(arg)
        kShuffle, // shuffle(arg, mask)
        kShuffle2, // shuffle2(a, b, mask)
        kDot4, // dot(a, b)
        kCross4, // cross(a, b)
        kMatMul, // arg * mat4
        kRemap, // remaps from source range (x->y) to destination range (x->y)
                // calculated as mix(float3(dest.x), float3(dest.y), (val - src.x)/(src.y - src.x))
        kSamplerBumpmap // samples normal from bump map texture
    };
}

// Map for conversion from RPR input types to Baikal input types (to save compatibility)
static const std::map<rpr_int, unsigned> kArithmeticOpRpr2Baikal =
{
    { RPR_MATERIAL_NODE_OP_ADD, Baikal::InputMapType::kAdd },
    { RPR_MATERIAL_NODE_OP_SUB, Baikal::InputMapType::kSub },
    { RPR_MATERIAL_NODE_OP_MUL, Baikal::InputMapType::kMul },
    { RPR_MATERIAL_NODE_OP_DIV, Baikal::InputMapType::kDiv },
    { RPR_MATERIAL_NODE_OP_SIN, Baikal::InputMapType::kSin },
    { RPR_MATERIAL_NODE_OP_COS, Baikal::InputMapType::kCos },
    { RPR_MATERIAL_NODE_OP_TAN, Baikal::InputMapType::kTan },
    { RPR_MATERIAL_NODE_OP_SELECT_X, Baikal::InputMapType::kSelect },
    { RPR_MATERIAL_NODE_OP_SELECT_Y, Baikal::InputMapType::kSelect },
    { RPR_MATERIAL_NODE_OP_SELECT_Z, Baikal::InputMapType::kSelect },
    { RPR_MATERIAL_NODE_OP_COMBINE, Baikal::InputMapType::kShuffle2 },
    { RPR_MATERIAL_NODE_OP_DOT3, Baikal::InputMapType::kDot3 },
    { RPR_MATERIAL_NODE_OP_CROSS3, Baikal::InputMapType::kCross3 },
    { RPR_MATERIAL_NODE_OP_LENGTH3, Baikal::InputMapType::kLength3 },
    { RPR_MATERIAL_NODE_OP_NORMALIZE3, Baikal::InputMapType::kNormalize3 },
    { RPR_MATERIAL_NODE_OP_POW, Baikal::InputMapType::kPow },
    { RPR_MATERIAL_NODE_OP_ACOS, Baikal::InputMapType::kAcos },
    { RPR_MATERIAL_NODE_OP_ASIN, Baikal::InputMapType::kAsin },
    { RPR_MATERIAL_NODE_OP_ATAN, Baikal::InputMapType::kAtan },
    { RPR_MATERIAL_NODE_OP_AVERAGE_XYZ, Baikal::InputMapType::kLerp },
    { RPR_MATERIAL_NODE_OP_AVERAGE, Baikal::InputMapType::kLerp },
    { RPR_MATERIAL_NODE_OP_MIN, Baikal::InputMapType::kMin },
    { RPR_MATERIAL_NODE_OP_MAX, Baikal::InputMapType::kMax },
    { RPR_MATERIAL_NODE_OP_FLOOR, Baikal::InputMapType::kFloor },
    { RPR_MATERIAL_NODE_OP_MOD, Baikal::InputMapType::kMod },
    { RPR_MATERIAL_NODE_OP_ABS, Baikal::InputMapType::kAbs },
    { RPR_MATERIAL_NODE_OP_SHUFFLE_YZWX, Baikal::InputMapType::kShuffle },
    { RPR_MATERIAL_NODE_OP_SHUFFLE_ZWXY, Baikal::InputMapType::kShuffle },
    { RPR_MATERIAL_NODE_OP_SHUFFLE_WXYZ, Baikal::InputMapType::kShuffle },
    { RPR_MATERIAL_NODE_OP_MAT_MUL, Baikal::InputMapType::kMatMul },
    { RPR_MATERIAL_NODE_OP_SELECT_W, Baikal::InputMapType::kSelect },
    { RPR_MATERIAL_NODE_OP_DOT4, Baikal::InputMapType::kDot4 },
    //{ RPR_MATERIAL_NODE_OP_LOG, ??? },          // Not present in Baikal
};

// XML based material IO implememtation
class MaterialIoXML : public MaterialIo
{
public:
    rpr_int LoadMaterials(rpr_char const* filename, rpr_context context, rpr_material_system material_system, std::map<std::string, rpr_material_node> & new_materials) override;
    rpr_int SaveMaterials(rpr_char const* filename, std::set<rpr_material_node> const& materials) override;

private:
    // Write single material
    rpr_int WriteMaterial(XMLPrinter& printer, const rpr_material_node material);

    rpr_int GetInputMapId(const rpr_material_node material, rpr_uint input_index, std::int64_t* out_id);
    // Write single InputMap
    rpr_int WriteInputMap(XMLPrinter& printer, const rpr_material_node material, rpr_uint input_index, std::int64_t* out_input_id = nullptr);
    rpr_int WriteFloatInput(XMLPrinter& printer, const rpr_material_node base_node, rpr_uint input_index, std::int64_t input_id);
    rpr_int WriteTextureInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id);
    rpr_int WriteArithmeticInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id);

    // Load input
    rpr_int LoadInput(rpr_context context, rpr_material_system material_system, rpr_material_node material,
        rpr_char const* input_name, std::int64_t input_id, std::map<std::int64_t, XMLElement*> const& xml_inputs);

    rpr_int LoadNodeInput(rpr_context context, rpr_material_system material_system, XMLElement* xml_input,
        std::map<std::int64_t, XMLElement*> const& xml_inputs, rpr_material_node* out_node);
    // Load single material
    rpr_int LoadMaterial(rpr_context context, rpr_material_system material_system, XMLElement& element, std::map<std::int64_t, XMLElement*> const& xml_input_maps, std::map<std::string, rpr_material_node> & out_materials);
    rpr_int LoadOneArgInput(rpr_context context, rpr_material_system material_system, rpr_uint operation,
        XMLElement* xml_input, std::map<std::int64_t, XMLElement*> const& xml_inputs, rpr_material_node* out_node);
    rpr_int LoadTwoArgInput(rpr_context context, rpr_material_system material_system, rpr_uint operation,
        XMLElement* xml_input, std::map<std::int64_t, XMLElement*> const& xml_inputs, rpr_material_node* out_node);

    // Texture to name map
    std::map<rpr_image, std::string> m_tex2name;

    std::map<std::string, rpr_image> m_name2tex;
    std::set<std::int64_t> m_saved_inputs;
    std::map<std::int64_t, rpr_material_node> m_loaded_inputs;

    std::string m_base_path;
    std::uint32_t m_current_material_index;

};

std::unique_ptr<MaterialIo> MaterialIo::CreateMaterialIoXML()
{
    return std::make_unique<MaterialIoXML>();
}

static std::string Float4ToString(rpr_float value[4])
{
    std::ostringstream oss;
    oss << value[0] << " " << value[1] << " " << value[2] << " " << value[3];
    return oss.str();
}

static std::string Matrix4x4ToString(rpr_float const matrix[4][4])
{
    std::ostringstream oss;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            oss << matrix[i][j] << " ";
        }
    }
    return oss.str();
}

static std::string ArrayToString(rpr_uint const* array, rpr_uint size)
{
    std::ostringstream oss;
    for (auto a = 0u; a < size; ++a)
    {
        oss << array[a] << " ";
    }
    return oss.str();
}

rpr_int MaterialIoXML::WriteMaterial(XMLPrinter& printer, const rpr_material_node material)
{
    printer.OpenElement("Material");

    rpr_char material_name[256];
    rpr_int status = rprMaterialNodeGetInfo(material, RPR_OBJECT_NAME, 0, material_name, nullptr);
    RETURN_IF_FAILED(status);

    printer.PushAttribute("name", material_name);
    printer.PushAttribute("id", m_current_material_index++);

    rpr_uint is_thin = 0;
    status = rprMaterialNodeGetInputInfo(material, RPR_UBER_MATERIAL_REFRACTION_THIN_SURFACE, RPR_MATERIAL_NODE_INPUT_VALUE, 0, &is_thin, nullptr);
    RETURN_IF_FAILED(status);
    printer.PushAttribute("thin", is_thin);

    rpr_uint refraction_ior_mode = 0;
    status = rprMaterialNodeGetInputInfo(material, RPR_UBER_MATERIAL_REFRACTION_IOR_MODE, RPR_MATERIAL_NODE_INPUT_VALUE, 0, &refraction_ior_mode, nullptr);
    RETURN_IF_FAILED(status);
    printer.PushAttribute("refraction_link_ior", refraction_ior_mode == RPR_UBER_MATERIAL_REFRACTION_MODE_LINKED);

    rpr_uint emission_mode = 0;
    status = rprMaterialNodeGetInputInfo(material, RPR_UBER_MATERIAL_EMISSION_MODE, RPR_MATERIAL_NODE_INPUT_VALUE, 0, &emission_mode, nullptr);
    RETURN_IF_FAILED(status);
    printer.PushAttribute("emission_doublesided", emission_mode == RPR_UBER_MATERIAL_EMISSION_MODE_DOUBLESIDED);

    rpr_uint sss_multiscatter = 0;
    status = rprMaterialNodeGetInputInfo(material, RPR_UBER_MATERIAL_SSS_MULTISCATTER, RPR_MATERIAL_NODE_INPUT_VALUE, 0, &sss_multiscatter, nullptr);
    RETURN_IF_FAILED(status);
    printer.PushAttribute("sss_multyscatter", sss_multiscatter);

    rpr_uint layers = 0;
    status = rprMaterialNodeGetInputInfo(material, RPR_UBER_MATERIAL_LAYERS, RPR_MATERIAL_NODE_INPUT_VALUE, 0, &layers, nullptr);
    RETURN_IF_FAILED(status);
    printer.PushAttribute("layers", layers);

    std::uint64_t num_inputs = 0;
    status = rprMaterialNodeGetInfo(material, RPR_MATERIAL_NODE_INPUT_COUNT, 0, &num_inputs, nullptr);
    RETURN_IF_FAILED(status);

    for (rpr_int i = 0; i < num_inputs; ++i)
    {
        // Skip "uberv2.layers"
        rpr_char input_name[256];
        status = rprMaterialNodeGetInputInfo(material, i, RPR_MATERIAL_NODE_INPUT_NAME_STRING, 0, &input_name, nullptr);
        RETURN_IF_FAILED(status);
        if (strcmp(input_name, "uberv2.layers") == 0)
        {
            continue;
        }

        std::int64_t input_id = -1;
        status = GetInputMapId(material, i, &input_id);
        RETURN_IF_FAILED(status);

        printer.PushAttribute(input_name, input_id);
    }

    printer.CloseElement();

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::SaveMaterials(rpr_char const* filename, std::set<rpr_material_node> const& materials)
{
    //std::string fname(filename);
    //auto slash = fname.find_last_of('/');
    //if (slash == std::string::npos)
    //{
    //    slash = fname.find_last_of('\\');
    //}

    //if (slash != std::string::npos)
    //{
    //    m_base_path.assign(fname.cbegin(), fname.cbegin() + slash + 1);
    //}
    //else
    //{
    //    m_base_path.clear();
    //}

    XMLDocument doc;
    XMLPrinter printer;

    rpr_int status = RPR_SUCCESS;

    printer.OpenElement("Inputs");

    for (rpr_material_node material : materials)
    {
        std::uint64_t num_inputs = 0;
        status = rprMaterialNodeGetInfo(material, RPR_MATERIAL_NODE_INPUT_COUNT, 0, &num_inputs, nullptr);
        RETURN_IF_FAILED(status);

        for (auto i = 0; i < num_inputs; ++i)
        {
            status = WriteInputMap(printer, material, i);
            RETURN_IF_FAILED(status);
        }
    }
    printer.CloseElement();

    printer.OpenElement("Materials");
    for (rpr_material_node material : materials)
    {
        status = WriteMaterial(printer, material);
        RETURN_IF_FAILED(status);
    }
    printer.CloseElement();

    doc.Parse(printer.CStr());
    doc.SaveFile(filename);

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::LoadInput(rpr_context context, rpr_material_system material_system, rpr_material_node material, 
    rpr_char const* input_name, std::int64_t input_id, std::map<std::int64_t, XMLElement*> const& xml_inputs)
{
    auto it = xml_inputs.find(input_id);
    if (it == xml_inputs.cend())
    {
        // Failed to find input
        return RPR_ERROR_IO_ERROR;
    }

    rpr_int status = RPR_SUCCESS;
    XMLElement const* xml_input = it->second;
    switch (xml_input->UnsignedAttribute("type"))
    {
    case 0:
    {
        // float3/float4 input
        rpr_float value[4];
        std::istringstream iss(xml_input->Attribute("value"));
        iss >> value[0] >> value[1] >> value[2] >> value[3];
        status = rprMaterialNodeSetInputF(material, input_name, value[0], value[1], value[2], value[3]);
        RETURN_IF_FAILED(status);
        break;
    }
    case 1:
    {
        // float input
        rpr_float value = xml_input->FloatAttribute("value");
        status = rprMaterialNodeSetInputF(material, input_name, value, value, value, value);
        RETURN_IF_FAILED(status);
        break;
    }
    default:
    {
        // Other input types - node input
        rpr_material_node input = nullptr;
        status = LoadNodeInput(context, material_system, it->second, xml_inputs, &input);
        RETURN_IF_FAILED(status);

        status = rprMaterialNodeSetInputN(material, input_name, input);
        RETURN_IF_FAILED(status);
        break;
    }
    }

    return RPR_SUCCESS;

}

rpr_int MaterialIoXML::LoadMaterial(rpr_context context, rpr_material_system material_system, XMLElement& element, std::map<std::int64_t, XMLElement*> const& xml_input_maps, std::map<std::string, rpr_material_node> & out_materials)
{
    rpr_material_node material = nullptr;
    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_UBERV2, &material);
    RETURN_IF_FAILED(status);

    status = rprObjectSetName(material, element.Attribute("name"));
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_LAYERS, element.UnsignedAttribute("layers"));
    RETURN_IF_FAILED(status);
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_REFRACTION_THIN_SURFACE, element.UnsignedAttribute("thin"));
    RETURN_IF_FAILED(status);
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_REFRACTION_IOR_MODE, element.BoolAttribute("refraction_link_ior"));
    RETURN_IF_FAILED(status);
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_EMISSION_MODE, element.BoolAttribute("emission_doublesided"));
    RETURN_IF_FAILED(status);
    status = rprMaterialNodeSetInputU_ext(material, RPR_UBER_MATERIAL_SSS_MULTISCATTER, element.UnsignedAttribute("sss_multyscatter"));
    RETURN_IF_FAILED(status);

    // Get uberv2 layer attributes
    for (XMLAttribute const* layer_attribute = element.FirstAttribute(); layer_attribute; layer_attribute = layer_attribute->Next())
    {
        char const* input_name = layer_attribute->Name();
        if (strncmp(input_name, "uberv2", 6) != 0)
        {
            continue;
        }

        std::int64_t child_id = layer_attribute->Int64Value();
        status = LoadInput(context, material_system, material, input_name, child_id, xml_input_maps);
        RETURN_IF_FAILED(status);
    }

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::LoadMaterials(rpr_char const* filename, rpr_context context, rpr_material_system material_system, std::map<std::string, rpr_material_node> & new_materials)
{
//    m_id2mat.clear();
//    m_name2tex.clear();
//
//    auto slash = file_name.find_last_of('/');
//    if (slash == std::string::npos) slash = file_name.find_last_of('\\');
//    if (slash != std::string::npos)
//        m_base_path.assign(file_name.cbegin(), file_name.cbegin() + slash + 1);
//    else
//        m_base_path.clear();
//
    XMLDocument doc;
    doc.LoadFile(filename);

    std::map<std::int64_t, XMLElement*> input_map_cache;
    auto inputs = doc.FirstChildElement("Inputs");
    for (auto element = inputs->FirstChildElement(); element; element = element->NextSiblingElement())
    {
        std::int64_t id = element->Int64Attribute("id");
        input_map_cache.insert(std::make_pair(id, element));
    }

    auto materials_node = doc.FirstChildElement("Materials");
    for (auto element = materials_node->FirstChildElement(); element; element = element->NextSiblingElement())
    {
        rpr_int status = LoadMaterial(context, material_system, *element, input_map_cache, new_materials);
        RETURN_IF_FAILED(status);
    }

    return RPR_SUCCESS;
}

rpr_int MaterialIo::SaveMaterialsFromScene(rpr_char const* filename, rpr_scene scene)
{
    std::size_t shape_count = 0;
    rpr_int status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_COUNT, sizeof(shape_count), &shape_count, nullptr);
    RETURN_IF_FAILED(status);

    std::vector<rpr_shape> shapes(shape_count);
    status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_LIST, sizeof(rpr_shape) * shape_count, shapes.data(), nullptr);
    RETURN_IF_FAILED(status);

    std::set<rpr_material_node> materials;

    for (rpr_shape shape : shapes)
    {
        rpr_material_node material = nullptr;
        status = rprShapeGetInfo(shape, RPR_SHAPE_MATERIAL, sizeof(material), &material, nullptr);
        RETURN_IF_FAILED(status);
        materials.insert(material);
    }

    status = SaveMaterials(filename, materials);
    RETURN_IF_FAILED(status);

    return RPR_SUCCESS;
}

rpr_int MaterialIo::ReplaceSceneMaterials(rpr_scene scene, std::map<std::string, rpr_material_node> new_materials, MaterialMap const& mapping)
{
    std::size_t shape_count = 0;
    rpr_int status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_COUNT, sizeof(shape_count), &shape_count, nullptr);
    RETURN_IF_FAILED(status);

    std::vector<rpr_shape> shapes(shape_count);
    status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_LIST, sizeof(rpr_shape) * shape_count, shapes.data(), nullptr);
    RETURN_IF_FAILED(status);

    for (rpr_shape shape : shapes)
    {
        rpr_material_node material = nullptr;
        status = rprShapeGetInfo(shape, RPR_SHAPE_MATERIAL, sizeof(material), &material, nullptr);
        RETURN_IF_FAILED(status);

        rpr_char material_name[256];
        status = rprMaterialNodeGetInfo(material, RPR_OBJECT_NAME, 0, material_name, nullptr);
        RETURN_IF_FAILED(status);

        auto it = mapping.find(material_name);
        if (it != mapping.cend())
        {
            auto mat_it = new_materials.find(it->second);

            if (mat_it != new_materials.cend())
            {
                status = rprShapeSetMaterial(shape, mat_it->second);
                RETURN_IF_FAILED(status);
            }
        }
    }

    return RPR_SUCCESS;
}

MaterialIo::MaterialMap MaterialIo::LoadMaterialMapping(rpr_char const* filename)
{
    MaterialMap map;

    XMLDocument doc;
    doc.LoadFile(filename);

    for (auto element = doc.FirstChildElement(); element; element = element->NextSiblingElement())
    {
        std::string from(element->Attribute("from"));
        std::string to(element->Attribute("to"));
        map.emplace(from, to);
    }

    return map;

}

rpr_int MaterialIo::SaveIdentityMapping(rpr_char const* filename, rpr_scene scene)
{
    std::size_t shape_count = 0;
    rpr_int status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_COUNT, sizeof(shape_count), &shape_count, nullptr);
    RETURN_IF_FAILED(status);

    std::vector<rpr_shape> shapes(shape_count);
    status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_LIST, sizeof(rpr_shape) * shape_count, shapes.data(), nullptr);
    RETURN_IF_FAILED(status);

    XMLPrinter printer;
    std::set<rpr_material_node> serialized_mats;
    for (rpr_shape shape : shapes)
    {
        rpr_material_node material = nullptr;
        status = rprShapeGetInfo(shape, RPR_SHAPE_MATERIAL, sizeof(material), &material, nullptr);
        RETURN_IF_FAILED(status);

        if (material && serialized_mats.find(material) == serialized_mats.cend())
        {
            rpr_char material_name[256];
            status = rprMaterialNodeGetInfo(material, RPR_OBJECT_NAME, 0, material_name, nullptr);
            RETURN_IF_FAILED(status);

            printer.OpenElement("Mapping");
            printer.PushAttribute("from", material_name);
            printer.PushAttribute("to", material_name);
            printer.CloseElement();
            serialized_mats.emplace(material);
        }
    }

    XMLDocument doc;
    doc.Parse(printer.CStr());
    doc.SaveFile(filename);

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::GetInputMapId(const rpr_material_node material, rpr_uint input_index, std::int64_t* out_id)
{
    rpr_uint input_type;
    rpr_int status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_TYPE, sizeof(input_type), &input_type, nullptr);
    RETURN_IF_FAILED(status);

    switch (input_type)
    {
    case RPR_MATERIAL_NODE_INPUT_TYPE_FLOAT4:
    {
        rpr_float value[4];
        status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(value), value, nullptr);
        RETURN_IF_FAILED(status);

        // Hash input value to get id
        std::ostringstream oss;
        oss << value[0] << value[1] << value[2] << value[3];
        std::hash<std::string> hash_fn;
        *out_id = static_cast<std::int64_t>(hash_fn(oss.str()) & 0x7FFF'FFFF'FFFF'FFFF);

        break;
    }
    case RPR_MATERIAL_NODE_INPUT_TYPE_NODE:
    {
        rpr_material_node input_node = nullptr;
        status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(input_node), &input_node, nullptr);
        RETURN_IF_FAILED(status);

        // Convert pointee address to id
        *out_id = reinterpret_cast<std::int64_t>(input_node);
        break;
    }
    default:
        assert(!"Material input type may only be node or float4!");
    }

    return RPR_SUCCESS;

}

rpr_int MaterialIoXML::WriteInputMap(XMLPrinter& printer, const rpr_material_node material, rpr_uint input_index, std::int64_t* out_input_id)
{
    // Skip "uberv2.layers"
    rpr_char input_name[256];
    rpr_int status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_NAME_STRING, sizeof(input_name), &input_name, nullptr);
    RETURN_IF_FAILED(status);

    if (strcmp(input_name, "uberv2.layers") == 0)
    {
        return RPR_SUCCESS;
    }

    // Get input id
    std::int64_t inputmap_id = -1;
    status = GetInputMapId(material, input_index, &inputmap_id);
    RETURN_IF_FAILED(status);

    // Skip if we have serialized this input before
    if (m_saved_inputs.find(inputmap_id) != m_saved_inputs.end())
    {
        if (out_input_id)
        {
            *out_input_id = inputmap_id;
        }

        return RPR_SUCCESS;
    }

    m_saved_inputs.insert(inputmap_id);

    // Get node type: is it a float4 or node input?
    rpr_uint input_type;
    status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_TYPE, sizeof(input_type), &input_type, nullptr);
    RETURN_IF_FAILED(status);

    if (input_type == RPR_MATERIAL_NODE_INPUT_TYPE_FLOAT4)
    {
        status = WriteFloatInput(printer, material, input_index, inputmap_id);
        RETURN_IF_FAILED(status);

    }
    else if (input_type == RPR_MATERIAL_NODE_INPUT_TYPE_NODE)
    {
        rpr_material_node input_node = nullptr;
        status = rprMaterialNodeGetInputInfo(material, input_index, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(input_node), &input_node, nullptr);
        RETURN_IF_FAILED(status);

        // Get type of node
        rpr_material_node_type node_type;
        status = rprMaterialNodeGetInfo(input_node, RPR_MATERIAL_NODE_TYPE, sizeof(node_type), &node_type, nullptr);
        RETURN_IF_FAILED(status);

        switch (node_type)
        {
        case RPR_MATERIAL_NODE_NORMAL_MAP:
        case RPR_MATERIAL_NODE_IMAGE_TEXTURE:
        case RPR_MATERIAL_NODE_NOISE2D_TEXTURE:
        case RPR_MATERIAL_NODE_DOT_TEXTURE:
        case RPR_MATERIAL_NODE_GRADIENT_TEXTURE:
        case RPR_MATERIAL_NODE_CHECKER_TEXTURE:
        case RPR_MATERIAL_NODE_CONSTANT_TEXTURE:
        case RPR_MATERIAL_NODE_BUMP_MAP:
            status = WriteTextureInput(printer, input_node, inputmap_id);
            RETURN_IF_FAILED(status);
            break;

        case RPR_MATERIAL_NODE_ARITHMETIC:
            status = WriteArithmeticInput(printer, input_node, inputmap_id);
            RETURN_IF_FAILED(status);
            break;

        case RPR_MATERIAL_NODE_UBERV2:
            assert(!"UberV2 material node cannot be set as input");
            break;

        default:
            assert(!"Something went wrong, other materials cannot be created");
            break;
        }
    }
    else
    {
        assert(!"Material input type may only be node or float4!");
    }

    if (out_input_id)
    {
        *out_input_id = inputmap_id;
    }

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::WriteFloatInput(XMLPrinter& printer, const rpr_material_node base_node, rpr_uint input_index, std::int64_t input_id)
{
    rpr_float value[4];
    rpr_int status = rprMaterialNodeGetInputInfo(base_node, input_index, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(value), value, nullptr);
    RETURN_IF_FAILED(status);

    printer.OpenElement("Input");
    {
        // Can't get name of float input in rpr layer
        printer.PushAttribute("name", "");
        printer.PushAttribute("id", input_id);

        // InputMapType::kConstantFloat3 = 0
        printer.PushAttribute("type", 0);
        printer.PushAttribute("value", Float4ToString(value).c_str());
    }
    printer.CloseElement();

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::WriteTextureInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id)
{
    // Here we're pretty sure that an image is connected to 0th input
    rpr_image image = nullptr;
    rpr_int status = rprMaterialNodeGetInputInfo(input_node, 0, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(image), &image, nullptr);
    RETURN_IF_FAILED(status);

    rpr_char image_name[256];
    status = rprImageGetInfo(image, RPR_OBJECT_NAME, sizeof(image_name), image_name, nullptr);
    RETURN_IF_FAILED(status);

    rpr_char input_name[256];
    status = rprMaterialNodeGetInfo(input_node, RPR_OBJECT_NAME, sizeof(input_name), input_name, nullptr);
    RETURN_IF_FAILED(status);

    printer.OpenElement("Input");
    {
        printer.PushAttribute("name", input_name);
        printer.PushAttribute("id", input_id);
        // InputMapType::kSampler = 2
        printer.PushAttribute("type", 2);
        // TODO: tex2name map?
        // TODO: store only relative path to image
        // TODO: save new image if it is not present on a drive
        printer.PushAttribute("value", image_name);
    }
    printer.CloseElement();

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::WriteArithmeticInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id)
{
    // Write child inputs and save their ids, also get type of arithmetic operation
    std::unordered_map<std::string, std::int64_t> child_ids;
    rpr_int op_type;

    std::uint64_t input_count = 0;
    rpr_int status = rprMaterialNodeGetInfo(input_node, RPR_MATERIAL_NODE_INPUT_COUNT, sizeof(input_count), &input_count, nullptr);
    RETURN_IF_FAILED(status);

    for (rpr_int i = 0; i < static_cast<rpr_int>(input_count); ++i)
    {
        rpr_char child_name[64];
        status = rprMaterialNodeGetInputInfo(input_node, i, RPR_MATERIAL_NODE_INPUT_NAME_STRING, sizeof(child_name), child_name, nullptr);

        if (strcmp(child_name, "op") == 0)
        {
            status = rprMaterialNodeGetInputInfo(input_node, i, RPR_MATERIAL_NODE_INPUT_VALUE, sizeof(op_type), &op_type, nullptr);
            RETURN_IF_FAILED(status);
        }
        else
        {
            std::int64_t child_id = -1;
            status = WriteInputMap(printer, input_node, i, &child_id);
            RETURN_IF_FAILED(status);
            child_ids[child_name] = child_id;
        }
    }

    rpr_char input_name[256];
    status = rprMaterialNodeGetInfo(input_node, RPR_OBJECT_NAME, sizeof(input_name), input_name, nullptr);
    RETURN_IF_FAILED(status);

    printer.OpenElement("Input");
    printer.PushAttribute("name", input_name);
    printer.PushAttribute("id", input_id);

    auto it = kArithmeticOpRpr2Baikal.find(op_type);
    if (it == kArithmeticOpRpr2Baikal.end())
    {
        // Unsupported arithmetic operation
        return RPR_ERROR_UNSUPPORTED;
    }

    printer.PushAttribute("type", it->second);

    switch (op_type)
    {
    case RPR_MATERIAL_NODE_OP_ADD:
    case RPR_MATERIAL_NODE_OP_SUB:
    case RPR_MATERIAL_NODE_OP_MUL:
    case RPR_MATERIAL_NODE_OP_DIV:
    case RPR_MATERIAL_NODE_OP_MIN:
    case RPR_MATERIAL_NODE_OP_MAX:
    case RPR_MATERIAL_NODE_OP_DOT3:
    case RPR_MATERIAL_NODE_OP_DOT4:
    case RPR_MATERIAL_NODE_OP_CROSS3:
    // kCross4 - not present in standard RPR
    case RPR_MATERIAL_NODE_OP_POW:
    case RPR_MATERIAL_NODE_OP_MOD:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("input1", child_ids["color1"]);
        break;

    case RPR_MATERIAL_NODE_OP_SIN:
    case RPR_MATERIAL_NODE_OP_COS:
    case RPR_MATERIAL_NODE_OP_TAN:
    case RPR_MATERIAL_NODE_OP_ASIN:
    case RPR_MATERIAL_NODE_OP_ACOS:
    case RPR_MATERIAL_NODE_OP_ATAN:
    case RPR_MATERIAL_NODE_OP_LENGTH3:
    case RPR_MATERIAL_NODE_OP_NORMALIZE3:
    case RPR_MATERIAL_NODE_OP_FLOOR:
    case RPR_MATERIAL_NODE_OP_ABS:
        printer.PushAttribute("input0", child_ids["color0"]);
        break;

    // kLerp - not present in standard RPR
    case RPR_MATERIAL_NODE_OP_SELECT_X:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 0u);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_Y:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 1u);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_Z:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 2u);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_W:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 3u);
        break;

    case RPR_MATERIAL_NODE_OP_AVERAGE_XYZ:
    case RPR_MATERIAL_NODE_OP_AVERAGE:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("input1", child_ids["color1"]);
        printer.PushAttribute("control", 0.5);
        break;

    case RPR_MATERIAL_NODE_OP_COMBINE:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("input1", child_ids["color1"]);
        printer.PushAttribute("mask", "0 4 1 5");
        break;

    case RPR_MATERIAL_NODE_OP_SHUFFLE_YZWX:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("mask", "1 2 3 0");
        break;

    case RPR_MATERIAL_NODE_OP_SHUFFLE_ZWXY:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("mask", "2 3 0 1");
        break;

    case RPR_MATERIAL_NODE_OP_SHUFFLE_WXYZ:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("mask", "3 0 1 2");
        break;

    case RPR_MATERIAL_NODE_OP_MAT_MUL:
        //printer.PushAttribute("input0", child_ids["color0"]);
        return RPR_ERROR_UNIMPLEMENTED;
        break;

    case RPR_MATERIAL_NODE_OP_LOG:
        // Not present in Baikal
        return RPR_ERROR_UNIMPLEMENTED;
        break;

    default:
        assert(!"Incorrect arithmetic operation type!");
        return RPR_ERROR_UNSUPPORTED;

    }
    printer.CloseElement();

    return RPR_SUCCESS;
}

rpr_int MaterialIoXML::LoadTwoArgInput(rpr_context context, rpr_material_system material_system, rpr_uint operation,
    XMLElement* xml_input, std::map<std::int64_t, XMLElement*> const& xml_inputs, rpr_material_node* out_node)
{
    assert(out_node);

    std::int64_t arg1_id = xml_input->Int64Attribute("input0");
    std::int64_t arg2_id = xml_input->Int64Attribute("input1");

    rpr_material_node material_node = nullptr;
    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_ARITHMETIC, &material_node);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputU(material_node, "op", operation);
    RETURN_IF_FAILED(status);

    status = LoadInput(context, material_system, material_node, "color0", arg1_id, xml_inputs);
    RETURN_IF_FAILED(status);

    status = LoadInput(context, material_system, material_node, "color1", arg2_id, xml_inputs);

    *out_node = material_node;
    return status;

}

rpr_int MaterialIoXML::LoadOneArgInput(rpr_context context, rpr_material_system material_system, rpr_uint operation,
    XMLElement* xml_input, std::map<std::int64_t, XMLElement*> const& xml_inputs, rpr_material_node* out_node)
{
    assert(out_node);
    std::int64_t arg1_id = xml_input->Int64Attribute("input0");

    rpr_material_node material_node = nullptr;
    rpr_int status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_ARITHMETIC, &material_node);
    RETURN_IF_FAILED(status);

    status = rprMaterialNodeSetInputU(material_node, "op", operation);
    RETURN_IF_FAILED(status);

    status = LoadInput(context, material_system, material_node, "color0", arg1_id, xml_inputs);
    return status;

}

rpr_int MaterialIoXML::LoadNodeInput(rpr_context context, rpr_material_system material_system,
    XMLElement * xml_input, std::map<std::int64_t, XMLElement*> const & xml_inputs, rpr_material_node * out_node)
{
    assert(out_node);

    std::int64_t id = xml_input->UnsignedAttribute("id");

    auto input = m_loaded_inputs.find(id);
    if (input != m_loaded_inputs.end())
    {
        *out_node = input->second;
        return RPR_SUCCESS;
    }

    rpr_material_node input_node = nullptr;
    rpr_int status = RPR_SUCCESS;
    rpr_uint input_type = xml_input->UnsignedAttribute("type");
    switch (input_type)
    {
    // Leafs
    case Baikal::InputMapType::kSampler:
    case Baikal::InputMapType::kSamplerBumpmap:
        {
            std::string image_filename(xml_input->Attribute("value"));

            auto iter = m_name2tex.find(image_filename);
            rpr_image image = nullptr;

            if (iter != m_name2tex.cend())
            {
                image = iter->second;
            }
            else
            {
                status = rprContextCreateImageFromFile(context, image_filename.c_str(), &image);
                RETURN_IF_FAILED(status);
                status = rprObjectSetName(image, image_filename.c_str());
                RETURN_IF_FAILED(status);
                m_name2tex.insert(std::make_pair(image_filename, image));
            }

            status = rprMaterialSystemCreateNode(material_system, RPR_MATERIAL_NODE_IMAGE_TEXTURE, &input_node);
            RETURN_IF_FAILED(status);

            status = rprMaterialNodeSetInputImageData(input_node, "data", image);
            RETURN_IF_FAILED(status);

            break;
        }
        // Two inputs
        case Baikal::InputMapType::kAdd:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_ADD, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kSub:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_SUB, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kMul:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_MUL, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kDiv:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_DIV, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kMin:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_MIN, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kMax:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_MAX, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kDot3:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_DOT3, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kDot4:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_DOT4, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kCross3:
        case Baikal::InputMapType::kCross4:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_CROSS3, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kPow:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_POW, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kMod:
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_MOD, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        //Single input
        case Baikal::InputMapType::kSin:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_SIN, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kCos:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_COS, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kTan:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_TAN, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kAsin:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_ASIN, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kAcos:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_ACOS, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kAtan:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_ATAN, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kLength3:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_LENGTH3, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kNormalize3:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_NORMALIZE3, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kFloor:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_FLOOR, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        case Baikal::InputMapType::kAbs:
            status = LoadOneArgInput(context, material_system, RPR_MATERIAL_NODE_OP_ABS, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        // Specials
        case Baikal::InputMapType::kLerp:
        {
            // Ignore "control" attribute!!! (Since we cannot set it in RPR layer)
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_AVERAGE, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        }
        case Baikal::InputMapType::kSelect:
        {
            std::uint32_t selection = xml_input->UnsignedAttribute("selection");
            rpr_uint operation;
            switch (selection)
            {
            case 0:
                operation = RPR_MATERIAL_NODE_OP_SELECT_X;
                break;
            case 1:
                operation = RPR_MATERIAL_NODE_OP_SELECT_Y;
                break;
            case 2:
                operation = RPR_MATERIAL_NODE_OP_SELECT_Z;
                break;
            case 3:
                operation = RPR_MATERIAL_NODE_OP_SELECT_W;
                break;
            default:
                // Unsupported RPR operation
                return RPR_ERROR_UNSUPPORTED;
            }

            status = LoadTwoArgInput(context, material_system, operation, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        }
        case Baikal::InputMapType::kShuffle:
        {
            rpr_uint operation;
            char const* mask = xml_input->Attribute("mask");
            if (strcmp(mask, "1 2 3 0") == 0)
            {
                operation = RPR_MATERIAL_NODE_OP_SHUFFLE_YZWX;
            }
            else if (strcmp(mask, "2 3 0 1") == 0)
            {
                operation = RPR_MATERIAL_NODE_OP_SHUFFLE_ZWXY;
            }
            else if (strcmp(mask, "3 0 1 2") == 0)
            {
                operation = RPR_MATERIAL_NODE_OP_SHUFFLE_WXYZ;
            }
            else
            {
                // Unsupported RPR operation
                return RPR_ERROR_UNSUPPORTED;
            }
            status = LoadOneArgInput(context, material_system, operation, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);

            break;
        }
        case Baikal::InputMapType::kShuffle2:
        {
            // Ignore "mask" attribute!!! (Since we cannot set it in RPR layer)
            status = LoadTwoArgInput(context, material_system, RPR_MATERIAL_NODE_OP_COMBINE, xml_input, xml_inputs, &input_node);
            RETURN_IF_FAILED(status);
            break;
        }
//        case Baikal::InputMapType::kMatMul:
//        {
//            uint32_t arg1_id = element->UnsignedAttribute("input0");
//            InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//
//            RadeonRays::matrix mat;
//            std::istringstream iss(element->Attribute("matrix"));
//            for (int i = 0; i < 4; ++i)
//            {
//                for (int j = 0; j < 4; ++j)
//                {
//                    iss >> mat.m[i][j];
//                }
//            }
//
//            result = InputMap_MatMul::Create(arg1, mat);
//            break;
//        }
//        case Baikal::InputMapType::kRemap:
//        {
//            uint32_t src_id = element->UnsignedAttribute("src");
//            InputMap::Ptr src = LoadInputMap(io, input_map_cache.at(src_id), input_map_cache, loaded_inputs);
//            uint32_t dst_id = element->UnsignedAttribute("dst");
//            InputMap::Ptr dst = LoadInputMap(io, input_map_cache.at(dst_id), input_map_cache, loaded_inputs);
//            uint32_t data_id = element->UnsignedAttribute("data");
//            InputMap::Ptr data = LoadInputMap(io, input_map_cache.at(data_id), input_map_cache, loaded_inputs);
//
//            result = InputMap_Remap::Create(src, dst, data);
//            break;
//        }
        default:
        {
            return RPR_ERROR_UNSUPPORTED;
        }
    }

    status = rprObjectSetName(input_node, xml_input->Attribute("name"));
    RETURN_IF_FAILED(status);

    m_loaded_inputs.insert(std::make_pair(id, input_node));

    *out_node = input_node;

    return RPR_SUCCESS;
}
