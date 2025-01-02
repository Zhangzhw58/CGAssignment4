/*
* 加载 .obj, .mtl和纹理图图像到模型
* 基于 https://github.com/distancemay5/learn_cg 继续完善
*/
#pragma once
#ifndef LOAD_MODEL_H
#define LOAD_MODEL_H

#define TINYOBJLOADER_IMPLEMENTATION
#include <string>
#include "hittable_list.h"
#include <unordered_map>
#include "external/tiny_obj_loader.h"
#include <cassert>
#include "triangle.h"
#include "material.h"
#include "texture.h"

// 加载模型
void load_model(std::string obj_path, hittable_list& world) {
    // 创建一个 ObjReader 对象
    tinyobj::ObjReader reader;

    // 解析 OBJ 文件
    if (!reader.ParseFromFile(obj_path)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        return;
    }

    if (!reader.Warning().empty()) {
        std::cout << "Warning: " << reader.Warning() << std::endl;
    }

    // 提取模型文件的目录部分
    std::string base_dir = obj_path.substr(0, obj_path.find_last_of('/') + 1);

    // 获取顶点属性、形状和材质
    const tinyobj::attrib_t& attrib = reader.GetAttrib();// 顶点属性，包括位置、法线和纹理坐标
    const std::vector<tinyobj::shape_t>& shapes = reader.GetShapes();//形状
    const std::vector<tinyobj::material_t>& materials = reader.GetMaterials();//材质

    // 材质映射
    // 材质映射(.obj材质名->通用材质类)
    std::unordered_map<std::string, std::shared_ptr<generic_material>> material_map;
    // 读取材质信息并存储到材质映射中
    for (const auto& mat : materials) {
        // 材质信息（可能加入更多）
        vec3 ambient(mat.ambient[0], mat.ambient[1], mat.ambient[2]);
        vec3 diffuse(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);
        vec3 specular(mat.specular[0], mat.specular[1], mat.specular[2]);
        double shininess = mat.shininess;
        double ref_idx = mat.ior; // 折射率

        std::shared_ptr<texture> diffuse_texture = nullptr;
        if (!mat.diffuse_texname.empty()) {
            std::string texture_path = base_dir + mat.diffuse_texname;
            diffuse_texture = std::make_shared<image_texture>(texture_path.c_str());
        }

        // 创建通用材质对象
        auto material = std::make_shared<generic_material>(ambient, diffuse, specular, shininess, ref_idx, diffuse_texture);
        material_map[mat.name] = material; // 插入映射
    }

    // 解析形状，创建三角形对象
    for (const auto& shape : shapes) {
        for (size_t i = 0; i < shape.mesh.indices.size(); i += 3) {
            tinyobj::index_t idx0 = shape.mesh.indices[i];
            tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
            tinyobj::index_t idx2 = shape.mesh.indices[i + 2];
            // 获取顶点位置
            vec3 v0_pos(attrib.vertices[3 * idx0.vertex_index + 0],
                attrib.vertices[3 * idx0.vertex_index + 1],
                attrib.vertices[3 * idx0.vertex_index + 2]);
            vec3 v1_pos(attrib.vertices[3 * idx1.vertex_index + 0],
                attrib.vertices[3 * idx1.vertex_index + 1],
                attrib.vertices[3 * idx1.vertex_index + 2]);
            vec3 v2_pos(attrib.vertices[3 * idx2.vertex_index + 0],
                attrib.vertices[3 * idx2.vertex_index + 1],
                attrib.vertices[3 * idx2.vertex_index + 2]);

            // 获取顶点法线（如果存在）
            glm::vec3 v0_normal(attrib.normals[3 * idx0.normal_index + 0],
                attrib.normals[3 * idx0.normal_index + 1],
                attrib.normals[3 * idx0.normal_index + 2]);
            glm::vec3 v1_normal(attrib.normals[3 * idx1.normal_index + 0],
                attrib.normals[3 * idx1.normal_index + 1],
                attrib.normals[3 * idx1.normal_index + 2]);
            glm::vec3 v2_normal(attrib.normals[3 * idx2.normal_index + 0],
                attrib.normals[3 * idx2.normal_index + 1],
                attrib.normals[3 * idx2.normal_index + 2]);

            // 获取纹理坐标（如果存在）
            glm::vec2 v0_tex(attrib.texcoords[2 * idx0.texcoord_index + 0],
                attrib.texcoords[2 * idx0.texcoord_index + 1]);
            glm::vec2 v1_tex(attrib.texcoords[2 * idx1.texcoord_index + 0],
                attrib.texcoords[2 * idx1.texcoord_index + 1]);
            glm::vec2 v2_tex(attrib.texcoords[2 * idx2.texcoord_index + 0],
                attrib.texcoords[2 * idx2.texcoord_index + 1]);

            // 获取材质
            int mat_id = shape.mesh.material_ids[i / 3];
            std::shared_ptr<material> material = nullptr;
            if (mat_id >= 0 && mat_id < materials.size()) {
                material = material_map[materials[mat_id].name];
            }

            // 颜色

            vertex v0(v0_pos.x(), v0_pos.y(), v0_pos.z(), v0_normal, glm::vec3(1.0, 1.0, 1.0), v0_tex);
            vertex v1(v1_pos.x(), v1_pos.y(), v1_pos.z(), v1_normal, glm::vec3(1.0, 1.0, 1.0), v1_tex);
            vertex v2(v2_pos.x(), v2_pos.y(), v2_pos.z(), v2_normal, glm::vec3(1.0, 1.0, 1.0), v2_tex);

            

            auto triangle_ = std::make_shared<triangle>(v0, v1, v2, material);
            world.add(triangle_);
        }
    }

    std::cout << "模型加载完成" << std::endl;
    
};

#endif