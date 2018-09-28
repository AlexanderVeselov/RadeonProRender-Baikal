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

#define RETURN_IF_FAILED(status) \
    if ((status) != RPR_SUCCESS) \
    {                            \
        return (status);         \
    }


// XML based material IO implememtation
class MaterialIoXML : public MaterialIo
{
public:
    rpr_int LoadMaterials(rpr_char const* filename, std::vector<rpr_material_node> & materials) override;
    rpr_int SaveMaterials(rpr_char const* filename, std::vector<rpr_material_node> const& materials) override;

private:
    // Write single material
    rpr_int WriteMaterial(XMLPrinter& printer, const rpr_material_node material);

    rpr_int GetInputMapId(const rpr_material_node material, rpr_uint input_index, std::int64_t* out_id);
    // Write single InputMap
    rpr_int WriteInputMap(XMLPrinter& printer, const rpr_material_node material, rpr_uint input_index, std::int64_t* out_input_id = nullptr);
    rpr_int WriteFloatInput(XMLPrinter& printer, const rpr_material_node base_node, rpr_uint input_index, std::int64_t input_id);
    rpr_int WriteTextureInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id);
    rpr_int WriteArithmeticInput(XMLPrinter& printer, const rpr_material_node input_node, std::int64_t input_id);


//    // Load inputs
//    InputMap::Ptr LoadInputMap(XMLElement* element,
//        const std::map<uint32_t, XMLElement*> &input_map_cache,
//        std::map<uint32_t, InputMap::Ptr> &loaded_inputs);
//    // Load single material
//    Material::Ptr LoadMaterial(XMLElement& element, const std::map<uint32_t, InputMap::Ptr> &loaded_inputs);

    // Texture to name map
    std::map<rpr_image, std::string> m_tex2name;

    std::map<std::string, rpr_image> m_name2tex;
    std::map<std::uint64_t, rpr_material_node> m_id2mat;
    std::set<std::int64_t> m_saved_inputs;

    std::string m_base_path;
    std::uint32_t m_current_material_index;
    std::uint32_t m_current_input_index;

//    template<class T>
//    InputMap::Ptr LoadTwoArgInput(ImageIo& io, XMLElement* element,
//        const std::map<uint32_t, XMLElement*> &input_map_cache,
//        std::map<uint32_t, InputMap::Ptr> &loaded_inputs)
//    {
//        uint32_t arg1_id = element->UnsignedAttribute("input0");
//        uint32_t arg2_id = element->UnsignedAttribute("input1");
//        InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//        InputMap::Ptr arg2 = LoadInputMap(io, input_map_cache.at(arg2_id), input_map_cache, loaded_inputs);

//        return T::Create(arg1, arg2);
//    }

//    template<class T>
//    InputMap::Ptr LoadOneArgInput(ImageIo& io, XMLElement* element,
//        const std::map<uint32_t, XMLElement*> &input_map_cache,
//        std::map<uint32_t, InputMap::Ptr> &loaded_inputs)
//    {
//        uint32_t arg1_id = element->UnsignedAttribute("input0");
//        InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);

//        return T::Create(arg1);
//    }

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

rpr_int MaterialIoXML::SaveMaterials(rpr_char const* filename, std::vector<rpr_material_node> const& materials)
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

    //m_tex2name.clear();
    rpr_int status = RPR_SUCCESS;

    m_current_input_index = 0;
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

    m_current_input_index = 0;
    printer.OpenElement("Materials");
    for (rpr_material_node material : materials)
    {
        status = WriteMaterial(printer, material);
        RETURN_IF_FAILED(status);
    }
    printer.CloseElement();

    doc.Parse(printer.CStr());
    doc.Parse(printer.CStr());

    doc.SaveFile(filename);

    return RPR_SUCCESS;
}

//Material::Ptr MaterialIoXML::LoadMaterial(ImageIo& io, XMLElement& element, const std::map<uint32_t, InputMap::Ptr> &loaded_inputs)
//{
//    std::string name(element.Attribute("name"));
//
//    auto attribute_thin = element.Attribute("thin");
//    std::string thin(attribute_thin ? attribute_thin : "");
//    auto id = static_cast<std::uint64_t>(std::atoi(element.Attribute("id")));
//
//    UberV2Material::Ptr material = UberV2Material::Create();
//    material->SetLayers(std::atoi(element.Attribute("layers")));
//    material->SetThin(thin == "true");
//    material->SetDoubleSided(strcmp(element.Attribute("emission_doublesided"), "true") == 0);
//    material->LinkRefractionIOR(strcmp(element.Attribute("refraction_link_ior"), "true") == 0);
//    material->SetMultiscatter(strcmp(element.Attribute("sss_multyscatter"), "true") == 0);
//    material->SetName(name);
//
//    auto num_inputs = material->GetNumInputs();
//    for (std::size_t a = 0u; a < num_inputs; ++a)
//    {
//        auto inputs = material->GetInput(a);
//        uint32_t input_id = element.UnsignedAttribute(inputs.info.name.c_str());
//        material->SetInputValue(inputs.info.name, loaded_inputs.at(input_id));
//    }
//
//    m_id2mat[id] = material;
//
//    return material;
//}

rpr_int MaterialIoXML::LoadMaterials(rpr_char const* filename, std::vector<rpr_material_node> & materials)
{
    return RPR_ERROR_UNIMPLEMENTED;


//    m_id2mat.clear();
//    m_name2tex.clear();
//    m_resolve_requests.clear();
//
//    auto slash = file_name.find_last_of('/');
//    if (slash == std::string::npos) slash = file_name.find_last_of('\\');
//    if (slash != std::string::npos)
//        m_base_path.assign(file_name.cbegin(), file_name.cbegin() + slash + 1);
//    else
//        m_base_path.clear();
//
//    XMLDocument doc;
//    doc.LoadFile(file_name.c_str());
//
//    auto image_io = ImageIo::CreateImageIo();
//
//    std::map<uint32_t, XMLElement*> input_map_cache;
//    auto inputs = doc.FirstChildElement("Inputs");
//    for (auto element = inputs->FirstChildElement(); element; element = element->NextSiblingElement())
//    {
//        uint32_t id = element->UnsignedAttribute("id");
//        input_map_cache.insert(std::make_pair(id, element));
//    }
//
//    std::map<uint32_t, InputMap::Ptr> loaded_elements;
//    for (auto element = inputs->FirstChildElement(); element; element = element->NextSiblingElement())
//    {
//        LoadInputMap(*image_io, element, input_map_cache, loaded_elements);
//    }
//
//    std::set<Material::Ptr> materials;
//    auto materials_node = doc.FirstChildElement("Materials");
//    for (auto element = materials_node->FirstChildElement(); element; element = element->NextSiblingElement())
//    {
//        auto material = LoadMaterial(*image_io, *element, loaded_elements);
//        materials.insert(material);
//    }
//
//    // Fix up non-resolved stuff
//    for (auto& i : m_resolve_requests)
//    {
//        i.material->SetInputValue(i.input, m_id2mat[i.id]);
//    }
//
//    return std::make_unique<ContainerIterator<std::set<Material::Ptr>>>(std::move(materials));
}

rpr_int MaterialIo::SaveMaterialsFromScene(rpr_char const* filename, rpr_scene scene)
{
    std::size_t shape_count = 0;
    rpr_int status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_COUNT, sizeof(shape_count), &shape_count, nullptr);
    RETURN_IF_FAILED(status);

    std::vector<rpr_shape> shapes(shape_count);
    status = rprSceneGetInfo(scene, RPR_SCENE_SHAPE_LIST, sizeof(rpr_shape) * shape_count, shapes.data(), nullptr);
    RETURN_IF_FAILED(status);

    // TODO: Check if materials not unique
    std::vector<rpr_material_node> materials;

    for (rpr_shape shape : shapes)
    {
        rpr_material_node material = nullptr;
        status = rprShapeGetInfo(shape, RPR_SHAPE_MATERIAL, sizeof(material), &material, nullptr);
        RETURN_IF_FAILED(status);
        materials.push_back(material);
    }

    status = SaveMaterials(filename, materials);
    RETURN_IF_FAILED(status);

    return RPR_SUCCESS;
}

//void MaterialIo::ReplaceSceneMaterials(Scene1& scene, Iterator& iterator, MaterialMap const& mapping)
//{
//    std::map<std::string, Material::Ptr> name2mat;
//
//    for (iterator.Reset(); iterator.IsValid(); iterator.Next())
//    {
//        auto material = iterator.ItemAs<Material>();
//        auto name = material->GetName();
//        name2mat[name] = material;
//    }
//
//    auto shape_iter = scene.CreateShapeIterator();
//
//    for (; shape_iter->IsValid(); shape_iter->Next())
//    {
//        auto shape = shape_iter->ItemAs<Shape>();
//        auto material = shape->GetMaterial();
//
//        if (!material)
//            continue;
//
//        auto name = material->GetName();
//        auto citer = mapping.find(name);
//
//        if (citer != mapping.cend())
//        {
//            auto mat_iter = name2mat.find(citer->second);
//
//            if (mat_iter != name2mat.cend())
//            {
//                shape->SetMaterial(mat_iter->second);
//            }
//        }
//    }
//}
//
//MaterialIo::MaterialMap MaterialIo::LoadMaterialMapping(std::string const& filename)
//{
//    MaterialMap map;
//
//    XMLDocument doc;
//    doc.LoadFile(filename.c_str());
//
//    for (auto element = doc.FirstChildElement(); element; element = element->NextSiblingElement())
//    {
//        std::string from(element->Attribute("from"));
//        std::string to(element->Attribute("to"));
//        map.emplace(from, to);
//    }
//
//    return map;
//}

//void MaterialIo::SaveIdentityMapping(std::string const& filename, Scene1 const& scene)
//{
//    XMLDocument doc;
//    XMLPrinter printer;
//
//    auto shape_iter = scene.CreateShapeIterator();
//    std::set<Material::Ptr> serialized_mats;
//
//    for (; shape_iter->IsValid(); shape_iter->Next())
//    {
//        auto material = shape_iter->ItemAs<Shape>()->GetMaterial();
//
//        if (material && serialized_mats.find(material) == serialized_mats.cend())
//        {
//            auto name = material->GetName();
//            printer.OpenElement("Mapping");
//            printer.PushAttribute("from", name.c_str());
//            printer.PushAttribute("to", name.c_str());
//            printer.CloseElement();
//            serialized_mats.emplace(material);
//        }
//    }
//
//    doc.Parse(printer.CStr());
//
//    doc.SaveFile(filename.c_str());
//}

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
        *out_id = static_cast<std::int64_t>(hash_fn(oss.str()));

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

        // Type: 0 � kConstantFloat3 (to save compatibility)
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
        // kSampler: type 2 � kSampler (to save compatibility)
        printer.PushAttribute("type", 2);
        // TODO: tex2name map
        // TODO: make image loaders to set only relative paths to names
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

    static const std::map<rpr_int, unsigned> kArithmeticOpRpr2Baikal =
    {
        { RPR_MATERIAL_NODE_OP_ADD, 3},             // InputMapType::kAdd = 3
        { RPR_MATERIAL_NODE_OP_SUB, 4},             // InputMapType::kSub = 4
        { RPR_MATERIAL_NODE_OP_MUL, 5},             // InputMapType::kMul = 5
        { RPR_MATERIAL_NODE_OP_DIV, 6},             // InputMapType::kDiv = 6
        { RPR_MATERIAL_NODE_OP_SIN, 7},             // InputMapType::kSin = 7
        { RPR_MATERIAL_NODE_OP_COS, 8},             // InputMapType::kCos = 8
        { RPR_MATERIAL_NODE_OP_TAN, 9},             // InputMapType::kTan = 9
        { RPR_MATERIAL_NODE_OP_SELECT_X, 10},       // InputMapType::kSelect = 10
        { RPR_MATERIAL_NODE_OP_SELECT_Y, 10},       // InputMapType::kSelect = 10
        { RPR_MATERIAL_NODE_OP_SELECT_Z, 10},       // InputMapType::kSelect = 10
        { RPR_MATERIAL_NODE_OP_COMBINE, 26},        // InputMapType::kShuffle2 = 26
        { RPR_MATERIAL_NODE_OP_DOT3, 11},           // InputMapType::kDot3 = 11
        { RPR_MATERIAL_NODE_OP_CROSS3, 12},         // InputMapType::kCross3 = 12
        { RPR_MATERIAL_NODE_OP_LENGTH3, 13},        // InputMapType::kLength3 = 13
        { RPR_MATERIAL_NODE_OP_NORMALIZE3, 14},     // InputMapType::kNormalize3 = 14
        { RPR_MATERIAL_NODE_OP_POW, 15},            // InputMapType::kPow = 15
        { RPR_MATERIAL_NODE_OP_ACOS, 16},           // InputMapType::kAcos = 16
        { RPR_MATERIAL_NODE_OP_ASIN, 17},           // InputMapType::kAsin = 17
        { RPR_MATERIAL_NODE_OP_ATAN, 18},           // InputMapType::kAtan = 18
        { RPR_MATERIAL_NODE_OP_AVERAGE_XYZ, 19},    // InputMapType::kLerp = 19
        { RPR_MATERIAL_NODE_OP_AVERAGE, 19},        // InputMapType::kLerp = 19
        { RPR_MATERIAL_NODE_OP_MIN, 20},            // InputMapType::kMin = 20
        { RPR_MATERIAL_NODE_OP_MAX, 21},            // InputMapType::kMax = 21
        { RPR_MATERIAL_NODE_OP_FLOOR, 22},          // InputMapType::kFloor = 22
        { RPR_MATERIAL_NODE_OP_MOD, 23},            // InputMapType::kMod = 23
        { RPR_MATERIAL_NODE_OP_ABS, 24},            // InputMapType::kAbs = 24
        { RPR_MATERIAL_NODE_OP_SHUFFLE_YZWX, 25},   // InputMapType::kShuffle = 25
        { RPR_MATERIAL_NODE_OP_SHUFFLE_ZWXY, 25},   // InputMapType::kShuffle = 25
        { RPR_MATERIAL_NODE_OP_SHUFFLE_WXYZ, 25},   // InputMapType::kShuffle = 25
        { RPR_MATERIAL_NODE_OP_MAT_MUL, 29},        // InputMapType::kMatMul = 29
        { RPR_MATERIAL_NODE_OP_SELECT_W, 10},       // InputMapType::kSelect = 10
        { RPR_MATERIAL_NODE_OP_DOT4, 27},           // InputMapType::kDot4 = 27
        //{ RPR_MATERIAL_NODE_OP_LOG, -1},          // Not present in Baikal
    };
    
    auto it = kArithmeticOpRpr2Baikal.find(op_type);
    if (it == kArithmeticOpRpr2Baikal.end())
    {
        // Unsupported operation
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
        printer.PushAttribute("selection", 0);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_Y:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 1);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_Z:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 2);
        break;

    case RPR_MATERIAL_NODE_OP_SELECT_W:
        printer.PushAttribute("input0", child_ids["color0"]);
        printer.PushAttribute("selection", 3);
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
        return RPR_ERROR_UNIMPLEMENTED;
        break;
    default:
        assert(!"Incorrect arithmetic operation type!");
        return RPR_ERROR_UNSUPPORTED;
    }
    printer.CloseElement();

    return RPR_SUCCESS;
}

//
//InputMap::Ptr MaterialIoXML::LoadInputMap(ImageIo& io, XMLElement* element,
//    const std::map<uint32_t, XMLElement*> &input_map_cache,
//    std::map<uint32_t, InputMap::Ptr> &loaded_inputs)
//{
//    InputMap::InputMapType type = static_cast<InputMap::InputMapType>(element->UnsignedAttribute("type"));
//    std::string name = element->Attribute("name");
//    int id = element->UnsignedAttribute("id");
//
//    auto input = loaded_inputs.find(id);
//    if (input != loaded_inputs.end())
//    {
//        return input->second;
//    }
//
//    InputMap::Ptr result;
//
//    switch (type)
//    {
//        //Leafs
//        case InputMap::InputMapType::kConstantFloat:
//        {
//            result = InputMap_ConstantFloat::Create(element->FloatAttribute("value"));
//            break;
//        }
//        case InputMap::InputMapType::kConstantFloat3:
//        {
//
//            std::istringstream iss(element->Attribute("value"));
//            RadeonRays::float3 value;
//            iss >> value.x >> value.y >> value.z;
//
//            result = InputMap_ConstantFloat3::Create(value);
//            break;
//        }
//        case InputMap::InputMapType::kSampler:
//        {
//            std::string filename(element->Attribute("value"));
//
//            auto iter = m_name2tex.find(filename);
//            Texture::Ptr texture;
//
//            if (iter != m_name2tex.cend())
//            {
//                texture = iter->second;
//            }
//            else
//            {
//                texture = io.LoadImage(m_base_path + filename);
//                m_name2tex[name] = texture;
//            }
//
//            result = InputMap_Sampler::Create(texture);
//
//            break;
//        }
//        case InputMap::InputMapType::kSamplerBumpmap:
//        {
//            std::string filename(element->Attribute("value"));
//
//            auto iter = m_name2tex.find(filename);
//            Texture::Ptr texture;
//
//            if (iter != m_name2tex.cend())
//            {
//                texture = iter->second;
//            }
//            else
//            {
//                texture = io.LoadImage(m_base_path + filename);
//                m_name2tex[name] = texture;
//            }
//
//            result = InputMap_SamplerBumpMap::Create(texture);
//            break;
//        }
//
//        // Two inputs
//        case InputMap::InputMapType::kAdd:
//            result = LoadTwoArgInput<InputMap_Add>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kSub:
//            result = LoadTwoArgInput<InputMap_Sub>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kMul:
//            result = LoadTwoArgInput<InputMap_Mul>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kDiv:
//            result = LoadTwoArgInput<InputMap_Div>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kMin:
//            result = LoadTwoArgInput<InputMap_Min>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kMax:
//            result = LoadTwoArgInput<InputMap_Max>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kDot3:
//            result = LoadTwoArgInput<InputMap_Dot3>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kDot4:
//            result = LoadTwoArgInput<InputMap_Dot4>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kCross3:
//            result = LoadTwoArgInput<InputMap_Cross3>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kCross4:
//            result = LoadTwoArgInput<InputMap_Cross4>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kPow:
//            result = LoadTwoArgInput<InputMap_Pow>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kMod:
//            result = LoadTwoArgInput<InputMap_Mod>(io, element, input_map_cache, loaded_inputs);
//            break;
//        //Single input
//        case InputMap::InputMapType::kSin:
//            result = LoadOneArgInput<InputMap_Sin>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kCos:
//            result = LoadOneArgInput<InputMap_Cos>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kTan:
//            result = LoadOneArgInput<InputMap_Tan>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kAsin:
//            result = LoadOneArgInput<InputMap_Asin>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kAcos:
//            result = LoadOneArgInput<InputMap_Acos>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kAtan:
//            result = LoadOneArgInput<InputMap_Atan>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kLength3:
//            result = LoadOneArgInput<InputMap_Length3>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kNormalize3:
//            result = LoadOneArgInput<InputMap_Normalize3>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kFloor:
//            result = LoadOneArgInput<InputMap_Floor>(io, element, input_map_cache, loaded_inputs);
//            break;
//        case InputMap::InputMapType::kAbs:
//            result = LoadOneArgInput<InputMap_Abs>(io, element, input_map_cache, loaded_inputs);
//            break;
//        // Specials
//        case InputMap::InputMapType::kLerp:
//        {
//            uint32_t arg1_id = element->UnsignedAttribute("input0");
//            uint32_t arg2_id = element->UnsignedAttribute("input1");
//            uint32_t control_id = element->UnsignedAttribute("control");
//            InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//            InputMap::Ptr arg2 = LoadInputMap(io, input_map_cache.at(arg2_id), input_map_cache, loaded_inputs);
//            InputMap::Ptr control = LoadInputMap(io, input_map_cache.at(control_id), input_map_cache, loaded_inputs);
//
//            result = InputMap_Lerp::Create(arg1, arg2, control);
//            break;
//        }
//        case InputMap::InputMapType::kSelect:
//        {
//            uint32_t arg1_id = element->UnsignedAttribute("input0");
//            InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//            InputMap_Select::Selection selection = 
//                static_cast<InputMap_Select::Selection>(element->UnsignedAttribute("selection"));
//
//            result = InputMap_Select::Create(arg1, selection);
//            break;
//        }
//        case InputMap::InputMapType::kShuffle:
//        {
//            uint32_t arg1_id = element->UnsignedAttribute("input0");
//            InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//            std::array<uint32_t, 4> mask;
//            std::istringstream iss(element->Attribute("mask"));
//            iss >> mask[0] >> mask[1] >> mask[2] >> mask[3];
//
//            result = InputMap_Shuffle::Create(arg1, mask);
//            break;
//        }
//        case InputMap::InputMapType::kShuffle2:
//        {
//            uint32_t arg1_id = element->UnsignedAttribute("input0");
//            InputMap::Ptr arg1 = LoadInputMap(io, input_map_cache.at(arg1_id), input_map_cache, loaded_inputs);
//            uint32_t arg2_id = element->UnsignedAttribute("input1");
//            InputMap::Ptr arg2 = LoadInputMap(io, input_map_cache.at(arg2_id), input_map_cache, loaded_inputs);
//
//            std::array<uint32_t, 4> mask;
//            std::istringstream iss(element->Attribute("mask"));
//            iss >> mask[0] >> mask[1] >> mask[2] >> mask[3];
//
//            result = InputMap_Shuffle2::Create(arg1, arg2, mask);
//            break;
//        }
//        case InputMap::InputMapType::kMatMul:
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
//        case InputMap::InputMapType::kRemap:
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
//    }
//
//    result->SetName(name);
//    loaded_inputs.insert(std::make_pair(id, result));
//    return result;
//}
