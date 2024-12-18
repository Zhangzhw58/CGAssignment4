#pragma once
#include"color.h"
#include"hittable_list.h"
#include"bvh.h"
#include"material.h"
#define CONTINUE_PROBABILITY 0.95
using namespace std;

color background = color(0.5, 0.7, 1.0); // ���豳����ɫ��ǳ��ɫ
const int max_depth = 50; // ������������

// ��������Ⱦ��ɫ
color ray_color(const ray& r, const hittable& world, int depth) {
	hit_record rec;

	// ����������������ṩ����
	if (depth <= 0)
		return color(0, 0, 0);

	if (world.hit(r, 0.0001, infinity, rec)) {
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(0, 0, 0);
	}
	vec3 unit_direction = unit_vector(r.direction()); // ���߷���
	auto t = 0.5 * (unit_direction.y() + 1.0); // ���ŵ�[0,1]
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0); // ����t��������ɫ���ֵ
}

color ray_color(const ray& r, const color& background, const hittable& world,
	int depth) {
	//RR
	double choice = random_double();
	// ����ﵽ�����ȣ��򷵻غ�ɫ
	if (depth <= 0)
		return color(0, 0, 0);

	hit_record rec;
	// �������û�л����κ����壬���ر���ɫ
	if (!world.hit(r, 0.001, infinity, rec))
		return background;

	ray scattered;  // ɢ���Ĺ���
	color attenuation; // ��ɫ˥��
	color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p); // �Է�����ɫ

	// ������ʲ�ɢ�䣬�����Է�����ɫ
	if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		return emitted;

	// ����ɢ���Ĺ��ױ���
	//double p = fmax(attenuation.x(), fmax(attenuation.y(), attenuation.z()));
	if (depth < max_depth) {
		// ʹ�� Russian Roulette ��ֹ����
		if (choice < CONTINUE_PROBABILITY) {
			return emitted + attenuation * ray_color(scattered, background, world, depth - 1) / CONTINUE_PROBABILITY;
		}
		else {
			return emitted;
		}
	}
	else {
		// ����ﵽ�����ȣ�ֱ�ӷ����ۻ�����ɫ
		return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
	}
}