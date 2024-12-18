#ifndef COLOR_H
#define COLOR_H

#include"vec3.h"
#include<iostream>

vec3 get_color(color pixel_color, int samples_per_pixel) {
	auto r = pixel_color.x();
	auto g = pixel_color.y();
	auto b = pixel_color.z();
	//根据样本数对颜色取平均值
	auto scale = 1.0 / samples_per_pixel;
	r = sqrt(r * scale);
	g = sqrt(g * scale);
	b = sqrt(b * scale);
	return vec3(r, g, b);
}



#endif // !COLOR_H