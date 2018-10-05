if (BAIKAL_NEXT_ENABLE_GLTF)
    set(GLTF_SOURCE_DIR ${BaikalNext_SOURCE_DIR}/3rdparty/RadeonProRender-GLTF/Source)

    set(GLTF_SOURCES
        ${GLTF_SOURCE_DIR}/BBox.cpp
        ${GLTF_SOURCE_DIR}/ImportBase.h
        ${GLTF_SOURCE_DIR}/ProRenderGLTF.h
        ${GLTF_SOURCE_DIR}/gltf2.h
        ${GLTF_SOURCE_DIR}/Import10.cpp
        ${GLTF_SOURCE_DIR}/ProRenderGLTF.cpp
        ${GLTF_SOURCE_DIR}/ImportBase.cpp
        ${GLTF_SOURCE_DIR}/Export.cpp
        ${GLTF_SOURCE_DIR}/Import11.cpp
        ${GLTF_SOURCE_DIR}/base64.h
        ${GLTF_SOURCE_DIR}/Import11.h
        ${GLTF_SOURCE_DIR}/Import10.h
        ${GLTF_SOURCE_DIR}/base64.cpp
        ${GLTF_SOURCE_DIR}/gltf2.cpp
    )
    set(GLTF_EXTENSIONS_SOURCES
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_lights.h
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_camera.cpp
        ${GLTF_SOURCE_DIR}/Extensions/Extensions.h
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_lights.cpp
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_heteroVolumes.h
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_post_effects.h
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_material.h
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_heteroVolumes.cpp
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_post_effects.cpp
        ${GLTF_SOURCE_DIR}/Extensions/KHR_materials_pbrSpecularGlossiness.h
        ${GLTF_SOURCE_DIR}/Extensions/KHR_materials_pbrSpecularGlossiness.cpp
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_material.cpp
        ${GLTF_SOURCE_DIR}/Extensions/AMD_RPR_camera.h
    )

    source_group("GLTF" FILES ${GLTF_SOURCES})
    source_group("GLTF\\Extensions" FILES ${GLTF_EXTENSIONS_SOURCES})

    add_library(RprGltf STATIC ${GLTF_SOURCES} ${GLTF_EXTENSIONS_SOURCES})

    target_include_directories(RprGltf
        PUBLIC ${GLTF_SOURCE_DIR}
        PRIVATE ${BaikalNext_SOURCE_DIR}/3rdparty/RadeonProRender-GLTF/3rdParty/json/include
        )

    find_package(Freeimage)
    target_link_libraries(RprGltf PUBLIC RadeonProRender64 RprSupport64 BaikalNextUtils Freeimage::Freeimage)
    target_compile_features(RprGltf PRIVATE cxx_std_14)

    if (UNIX AND NOT APPLE)
        target_link_libraries(RprGltf PUBLIC stdc++fs)
    elseif (APPLE)
        target_link_libraries(RprGltf PUBLIC c++experimental)
    endif ()

    set_property(TARGET RprGltf PROPERTY POSITION_INDEPENDENT_CODE ON)

endif (BAIKAL_NEXT_ENABLE_GLTF)

set(SCENE_IO_SOURCES
    scene_io.h
    scene_obj_io.cpp
    scene_test_io.cpp
    )

if (BAIKAL_NEXT_ENABLE_GLTF)
    list(APPEND SCENE_IO_SOURCES scene_gltf_io.cpp)
endif (BAIKAL_NEXT_ENABLE_GLTF)

if (BAIKAL_NEXT_ENABLE_FBX)
    list(APPEND SCENE_IO_SOURCES scene_fbx_io.cpp)
endif (BAIKAL_NEXT_ENABLE_FBX)

set(MATERIAL_IO_SOURCES
    material_io.cpp
    material_io.h
    )

set(XML_SOURCES
    XML/tinyxml2.cpp
    XML/tinyxml2.h
    )

set(SOURCES
    ${SCENE_IO_SOURCES}
    ${MATERIAL_IO_SOURCES}
    ${XML_SOURCES}
    RadeonProRenderIO.cpp
    RadeonProRenderIO.h
)

source_group("SceneIO" FILES ${SCENE_IO_SOURCES})
source_group("MaterialIO" FILES ${MATERIAL_IO_SOURCES})
source_group("XML" FILES ${XML_SOURCES})

add_library(RadeonProRenderIO SHARED ${SOURCES})
target_compile_definitions(RadeonProRenderIO PRIVATE RPRIO_EXPORT_API)
target_compile_features(RadeonProRenderIO PRIVATE cxx_std_14)
target_include_directories(RadeonProRenderIO PUBLIC . )
target_include_directories(RadeonProRenderIO PRIVATE ${BaikalNext_SOURCE_DIR}/3rdparty/tinyobjloader)
target_link_libraries(RadeonProRenderIO PRIVATE RadeonProRender64 RprSupport64 OpenImageIO::OpenImageIO)

if (BAIKAL_NEXT_ENABLE_GLTF)
    target_link_libraries(RadeonProRenderIO PRIVATE RprGltf)
endif (BAIKAL_NEXT_ENABLE_GLTF)

if (BAIKAL_NEXT_ENABLE_FBX)
    target_link_libraries(RadeonProRenderIO PRIVATE fbxsdk::fbxsdk)
endif (BAIKAL_NEXT_ENABLE_FBX)

if (UNIX AND NOT APPLE)
    target_link_libraries(RadeonProRenderIO PUBLIC stdc++fs)
elseif (APPLE)
    target_link_libraries(RadeonProRenderIO PUBLIC c++experimental)
endif ()

if (WIN32)
    target_compile_options(RadeonProRenderIO PRIVATE /WX)
elseif (UNIX)
    target_compile_options(RadeonProRenderIO PRIVATE -Wall -Werror)
endif (WIN32)

# Install target
if (WIN32)
    install(TARGETS RadeonProRenderIO RUNTIME DESTINATION bin)
else ()
    install(TARGETS RadeonProRenderIO LIBRARY DESTINATION lib)
endif ()
