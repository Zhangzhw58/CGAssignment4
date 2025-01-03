#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H
#include "utils.h"
#include "hittable.h"
#include "texture.h"
struct hit_record;

// ����
class material {
public:
	// ɢ��
	virtual color emitted(double u, double v, const point3& p) const {
		// û�й��߷��غ�ɫ
		return color(0, 0, 0);
	}
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray&
		scattered
	) const = 0;
};

// ���Դ
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
// �������ϣ����淴��
class metal : public material {
public:
	metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {} // ģ���Ŷ��̶�<=1
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
	double fuzz; // ģ���Ŷ�
};

// ����ʲ��ϣ���ȫ����
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
		// �ж�ȫ����
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
// �ʲ�����, ɢ�����йⲢ�Է����� R ˥��
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

// ����ͬ�Բ���(�����)
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

// ͨ�ò����� (of .obj)
// ���������⡢������⡢���淴��⡢����ȺͿ�ѡ������������ͼƬ��
// ���ݹ���ģ�ͽ��ж�Ӧ������
class generic_material : public material {
public:
	// ���캯��
	generic_material(const vec3& ambient, const vec3& diffuse, const vec3& specular, double shininess, double ref_idx, std::shared_ptr<texture> diffuse_texture, const int illum, const vec3& emitted)
		: ambient(ambient), diffuse(diffuse), specular(specular), shininess(shininess), ref_idx(ref_idx), diffuse_texture(diffuse_texture), illum(illum), light_color(emitted) {}

	// ����������ཻ���ɢ��
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		// ����������ߵĵ�λ��������
		vec3 unit_direction = unit_vector(r_in.direction());
		// ���㷴�䷽��
		vec3 reflected = reflect(unit_direction, rec.normal);
		// ��������ǵ�����ֵ
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		// ��������ǵ�����ֵ
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);


		// �ж��Ƿ���ȫ����
		bool cannot_refract = ref_idx * sin_theta > 1.0;
		vec3 direction;

		// ���ݹ���ģ��ѡ��ͬ��ɢ����Ϊ
		switch (illum) {
		case 4: // �߹ⷴ���͸����
			if (cannot_refract || reflectance(cos_theta, ref_idx) > random_double()) {
				direction = reflected + fuzz() * random_in_unit_sphere(); // ���ģ���Ŷ�;
			}
			else {
				direction = refract(unit_direction, rec.normal, ref_idx);
			}
			break;
		case 7: // �����͸���ȣ���Ϊ��Դ��
			direction = refract(unit_direction, rec.normal, ref_idx);
			return false;
			break;
		case 3: // ���淴��(���ģ���Ŷ�)
			direction = reflected + fuzz() * random_in_unit_sphere(); // ���ģ���Ŷ�;
			break;
		default:
			// Ĭ��������
			direction = rec.normal + random_unit_vector();
			break;
		}

		// �����µ�ɢ�����
		scattered = ray(rec.p, direction, r_in.time());

		// ���������ͼƬ
		if (diffuse_texture) {
			// ʹ������ͼƬ����ɫֵ��Ϊ˥��
			attenuation = diffuse_texture->value(rec.u, rec.v, rec.p);
		}
		else {
			// ����ʹ������������ɫ
			attenuation = diffuse;
		}

		return true;
	}

	// �Է�����ɫ
	virtual vec3 emitted(double u, double v, const vec3& p) const override {
		return light_color;
	}

	// ��������
	vec3 ambient; // ������
	vec3 diffuse; // �������
	vec3 specular; // ���淴���
	double shininess; // �����(���淴��ָ��)
	double ref_idx; // ������
	std::shared_ptr<texture> diffuse_texture; // ����������
	vec3 light_color; // �Է�����ɫ
	int illum; // ����ģ��

private:
	// ʹ�� Schlick �Ľ��ƹ�ʽ���㷴����
	static double reflectance(double cosine, double ref_idx) {
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}

	// ����ģ������
	double fuzz() const {
		return 1.0 / (shininess + 1.0);
	}
};

#endif