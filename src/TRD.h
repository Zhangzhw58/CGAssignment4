#pragma once
#include"color.h"
#include"hittable_list.h"
#include"bvh.h"
#include"material.h"
#define CONTINUE_PROBABILITY 0.95
using namespace std;

color background = color(0.5, 0.7, 1.0); // 假设背景颜色是浅蓝色
const int max_depth = 50; // 最大漫反射次数

// 射线上渲染颜色
color ray_color(const ray& r, const hittable& world, int depth) {
	hit_record rec;

	// 如果超过最大深度则不提供贡献
	if (depth <= 0)
		return color(0, 0, 0);

	if (world.hit(r, 0.0001, infinity, rec)) {
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(0, 0, 0);
	}
	vec3 unit_direction = unit_vector(r.direction()); // 射线方向
	auto t = 0.5 * (unit_direction.y() + 1.0); // 缩放到[0,1]
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0); // 根据t在两个颜色间插值
}

color ray_color(const ray& r, const color& background, const hittable& world,
	int depth) {
	//RR
	double choice = random_double();
	// 如果达到最大深度，则返回黑色
	if (depth <= 0)
		return color(0, 0, 0);

	hit_record rec;
	// 如果光线没有击中任何物体，返回背景色
	if (!world.hit(r, 0.001, infinity, rec))
		return background;

	ray scattered;  // 散射后的光线
	color attenuation; // 颜色衰减
	color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p); // 自发光颜色

	// 如果材质不散射，返回自发光颜色
	if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		return emitted;

	// 计算散射光的贡献比例
	//double p = fmax(attenuation.x(), fmax(attenuation.y(), attenuation.z()));
	if (depth < max_depth) {
		// 使用 Russian Roulette 终止策略
		if (choice < CONTINUE_PROBABILITY) {
			return emitted + attenuation * ray_color(scattered, background, world, depth - 1) / CONTINUE_PROBABILITY;
		}
		else {
			return emitted;
		}
	}
	else {
		// 如果达到最大深度，直接返回累积的颜色
		return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
	}
}