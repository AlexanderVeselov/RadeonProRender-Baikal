/**********************************************************************
 Copyright (c) 2018 Advanced Micro Devices, Inc. All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 ********************************************************************/

#pragma once

#include "RadeonProRender.h"
#include "RprSupport.h"

#ifdef WIN32
#ifdef RPRIO_EXPORT_API
#define RPRIO_API_ENTRY __declspec(dllexport)
#else
#define RPRIO_API_ENTRY __declspec(dllimport)
#endif
#else
#define RPRIO_API_ENTRY __attribute__((visibility ("default")))
#endif

/** @brief Load scene from file
*
*  @param  filename         Scene filename.
*  @param  basepath         Path to scene resources.
*  @param  context          The pre-initialized Radeon ProRender context handle to create API objects from.
*  @param  materialSystem   The pre-initialized Radeon ProRender material system handle to create API objects from.
*  @param  uberMatContext   The pre-initialized Radeon ProRender uber material system context handle to create API objects from.
*  @param  scene            The loaded scene and stored in this handle.
*  @return                  RPR_SUCCESS in case of success, error code otherwise
*/

extern RPRIO_API_ENTRY rpr_int rprLoadScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene* scene);

/** @brief Save scene to file
*
*  @param  filename         Scene filename.
*  @param  basepath         Path to scene resources.
*  @param  context          The pre-initialized Radeon ProRender context handle to create API objects from.
*  @param  materialSystem   The pre-initialized Radeon ProRender material system handle to create API objects from.
*  @param  uberMatContext   The pre-initialized Radeon ProRender uber material system context handle to create API objects from.
*  @param  scene            The scene to be written.
*  @return                  RPR_SUCCESS in case of success, error code otherwise
*/

extern RPRIO_API_ENTRY rpr_int rprSaveScene(rpr_char const* filename, rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rprx_context uberMatContext, rpr_scene scene);

/** @brief Replace scene materials based on descriptions provided in XML files
*
*  @param  materials_xml    File with descriptions of new materials.
*  @param  mapping_xml      File that contains mapping from old material names to new material names.
*  @param  basepath         Path to scene resources.
*  @param  context          The pre-initialized Radeon ProRender context handle to create API objects from.
*  @param  materialSystem   The pre-initialized Radeon ProRender material system handle to create API objects from.
*  @param  scene            The scene handle.
*  @return                  RPR_SUCCESS in case of success, error code otherwise
*/

extern RPRIO_API_ENTRY rpr_int rprReplaceSceneMaterials(rpr_char const* materials_xml, rpr_char const* mapping_xml, rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rpr_scene scene);

/** @brief Save new material configuration to XML files
*
*  @param  materials_xml    File with descriptions of new materials.
*  @param  mapping_xml      File that contains mapping from old material names to new material names.
*  @param  basepath         Path to scene resources.
*  @param  context          The pre-initialized Radeon ProRender context handle to create API objects from.
*  @param  materialSystem   The pre-initialized Radeon ProRender material system handle to create API objects from.
*  @param  scene            The scene handle.
*  @return                  RPR_SUCCESS in case of success, error code otherwise
*/

extern RPRIO_API_ENTRY rpr_int rprSaveSceneMaterials(rpr_char const* materials_xml, rpr_char const* mapping_xml, rpr_char const* basepath, rpr_context context, rpr_material_system materialSystem, rpr_scene scene);
