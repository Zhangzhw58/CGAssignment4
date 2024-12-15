#ifndef MATERIAL_H
#define MATERIAL_H
#include "utils.h"

struct hit_record;

// 材料
class material {
public:
	// 散射
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const = 0;
};

// 朗伯材料, 散射所有光并以反射率 R 衰减
class lambertian : public material {
public:
	lambertian(const color& a) : albedo(a) {}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		auto scatter_direction = rec.normal + random_unit_vector(); // 散射方向

		// 如果在0附近，改成法线方向
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;

		scattered = ray(rec.p, scatter_direction);
		attenuation = albedo;
		return true;
	}
public:
	color albedo;
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
		scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere()); // 反射方向添加模糊扰动
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
		vec3 refracted = refract(unit_direction, rec.normal,
			refraction_ratio);
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
		scattered = ray(rec.p, direction);
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


#endif
