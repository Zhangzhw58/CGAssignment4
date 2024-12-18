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
#endif