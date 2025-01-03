#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H
#include "utils.h"
#include "hittable.h"
#include "texture.h"
struct hit_record;

// 材料
class material {
public:
	// 散射
	virtual color emitted(double u, double v, const point3& p) const {
		// 没有光线返回黑色
		return color(0, 0, 0);
	}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const = 0;
};

// 点光源
class point_light : public material {
public:
	point_light(const point3& p, const color& c) : position(p), color(c) {}

	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		return false;
	}

	virtual color emitted(double u, double v, const point3& p) const override {
		double d = distance(p);
		return color / (d * d);
	}

	double distance(const point3& p) const {
		return (position - p).length();
	}

	bool in_shadow(const point3& p, const hittable& world) const {
		ray shadow_ray(p, position - p);
		hit_record rec;
		if (world.hit(shadow_ray, 0.001, infinity, rec)) {
			return true;
		}
		return false;
	}

public:
	point3 position;
	color color;
};
// 金属材料，镜面反射
class metal : public material {
public:
	metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {} // 模糊扰动固定<=1
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const override {
		vec3 reflected = reflect(unit_vector(r_in.direction()),
			rec.normal);
		scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere(),
			r_in.time());
		attenuation = albedo;
		return (dot(scattered.direction(), rec.normal) > 0);
	}
public:
	color albedo;
	double fuzz; // 模糊扰动
};

// 电介质材料，完全折射
class dielectric : public material {
public:
	dielectric(double index_of_refraction) : ir(index_of_refraction) {}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const override {
		attenuation = color(1.0, 1.0, 1.0);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;
		vec3 unit_direction = unit_vector(r_in.direction());
		// 判断全反射
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;
		if (cannot_refract || reflectance(cos_theta, refraction_ratio) >
			random_double())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal,
				refraction_ratio);
		scattered = ray(rec.p, direction, r_in.time());
		return true;
	}
public:
	double ir; // Index of Refraction
private:
	static double reflectance(double cosine, double ref_idx) {
		// Use Schlick's approximation for reflectance.
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};
// 朗伯材料, 散射所有光并以反射率 R 衰减
class lambertian : public material {
public:
	lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
	lambertian(shared_ptr<texture> a) : albedo(a) {}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const override {
		auto scatter_direction = rec.normal + random_unit_vector();
		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		scattered = ray(rec.p, scatter_direction, r_in.time());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}
public:
	shared_ptr<texture> albedo;
};
class diffuse_light : public material {
public:
	diffuse_light(shared_ptr<texture> a) : emit(a) {}
	diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const override {
		return false;
	}
	virtual color emitted(double u, double v, const point3& p) const
		override {
		return emit->value(u, v, p);
	}
public:
	shared_ptr<texture> emit;
};

// 各向同性材料(烟雾等)
class isotropic : public material {
public:
	isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
	isotropic(shared_ptr<texture> a) : albedo(a) {}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const override {
		scattered = ray(rec.p, random_in_unit_sphere(), r_in.time());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}
public:
	shared_ptr<texture> albedo;
};

// 通用材质类 (of .obj)
// 包含环境光、漫反射光、镜面反射光、光泽度和可选的漫反射纹理（图片）
// 根据光照模型进行对应的运算
class generic_material : public material {
public:
	// 构造函数
	generic_material(const vec3& ambient, const vec3& diffuse, const vec3& specular, double shininess, double ref_idx, std::shared_ptr<texture> diffuse_texture, const int illum, const vec3& emitted)
		: ambient(ambient), diffuse(diffuse), specular(specular), shininess(shininess), ref_idx(ref_idx), diffuse_texture(diffuse_texture), illum(illum), light_color(emitted) {}

	// 光线与材质相交后的散射
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		// 计算入射光线的单位方向向量
		vec3 unit_direction = unit_vector(r_in.direction());
		// 计算反射方向
		vec3 reflected = reflect(unit_direction, rec.normal);
		// 计算入射角的余弦值
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		// 计算入射角的正弦值
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);


		// 判断是否发生全反射
		bool cannot_refract = ref_idx * sin_theta > 1.0;
		vec3 direction;

		// 根据光照模型选择不同的散射行为
		switch (illum) {
		case 4: // 高光反射和透明度
			if (cannot_refract || reflectance(cos_theta, ref_idx) > random_double()) {
				direction = reflected + fuzz() * random_in_unit_sphere(); // 添加模糊扰动;
			}
			else {
				direction = refract(unit_direction, rec.normal, ref_idx);
			}
			break;
		case 7: // 反射和透明度（即为光源）
			direction = refract(unit_direction, rec.normal, ref_idx);
			return false;
			break;
		case 3: // 镜面反射(添加模糊扰动)
			direction = reflected + fuzz() * random_in_unit_sphere(); // 添加模糊扰动;
			break;
		default:
			// 默认漫反射
			direction = rec.normal + random_unit_vector();
			break;
		}

		// 生成新的散射光线
		scattered = ray(rec.p, direction, r_in.time());

		// 如果有纹理图片
		if (diffuse_texture) {
			// 使用纹理图片的颜色值作为衰减
			attenuation = diffuse_texture->value(rec.u, rec.v, rec.p);
		}
		else {
			// 否则使用漫反射光的颜色
			attenuation = diffuse;
		}

		return true;
	}

	// 自发光颜色
	virtual vec3 emitted(double u, double v, const vec3& p) const override {
		return light_color;
	}

	// 材质属性
	vec3 ambient; // 环境光
	vec3 diffuse; // 漫反射光
	vec3 specular; // 镜面反射光
	double shininess; // 光泽度(镜面反射指数)
	double ref_idx; // 折射率
	std::shared_ptr<texture> diffuse_texture; // 漫反射纹理
	vec3 light_color; // 自发光颜色
	int illum; // 光照模型

private:
	// 使用 Schlick 的近似公式计算反射率
	static double reflectance(double cosine, double ref_idx) {
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}

	// 计算模糊因子
	double fuzz() const {
		return 1.0 / (shininess + 1.0);
	}
};

#endif