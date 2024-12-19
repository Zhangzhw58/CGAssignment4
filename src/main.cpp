/*The MIT License (MIT)

Copyright (c) 2021-Present, Wencong Yang (yangwc3@mail2.sysu.edu.cn).

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.*/

#include <array>
#include <vector>
#include <thread>
#include <iostream>

#include "WindowsApp.h"
#include "utils.h"
#include "hittable_list.h"

#include "material.h"
#include "bvh.h"
#include "aarect.h"
#include "texture.h"
#include "color.h"
#include "TRD.h"
#include "box.h"
#include "constant_media.h"
#include "moving_sphere.h"

static std::vector<std::vector<color>> gCanvas;		//Canvas
const std::string save_path = "output.bmp";			//运行结果保存路径

// The width and height of the screen
const auto aspect_ratio = 3.0 / 2.0;/*3.0 / 2.0;*/
const int gWidth = 800;
const int gHeight = static_cast<int>(gWidth / aspect_ratio);
const int samples_per_pixel = 5000; //500; // 每点采样数

void rendering();

// 光线
color ray_color(const ray& r, const hittable& world, int depth);
color ray_color(const ray& r, const color& background, const hittable& world,
	int depth);

// 射线与球相交
double hit_sphere(const point3& center, double radius, const ray& r) {
	vec3 oc = r.origin() - center; // 球心
	// 解二元一次方程
	auto a = r.direction().length_squared();
	auto half_b = dot(oc, r.direction());
	auto c = oc.length_squared() - radius * radius;
	auto discriminant = half_b * half_b - a * c;

	if (discriminant < 0) {
		return -1.0;
	}
	else {
		return (-half_b - sqrt(discriminant)) / a;
	}
}

// 球体类
class sphere : public hittable {
public:
	sphere() : center(point3(0, 0, 0)), radius(0), mat_ptr(nullptr){};
	sphere(point3 cen, double r, shared_ptr<material> m)
		: center(cen), radius(r), mat_ptr(m) {};
	virtual bool hit(
		const ray& r, double t_min, double t_max, hit_record& rec) const
		override;
	virtual bool bounding_box(double time0, double time1, aabb& output_box)
		const override;
public:
	point3 center; // 球心
	double radius; // 半径
	shared_ptr<material> mat_ptr; // 材料
	
private:
	// 纹理坐标
	static void get_sphere_uv(const point3& p, double& u, double& v) {
		// p: a given point on the sphere of radius one, centered at the origin.
		// u: returned value [0,1] of angle around the Y axis from X=-1.
		// v: returned value [0,1] of angle from Y=-1 to Y=+1.
		// <1 0 0> yields <0.50 0.50> <-1 0 0> yields <0.00 0.50 >
		// <0 1 0> yields <0.50 1.00> < 0 -1 0> yields <0.50 0.00 >
		// <0 0 1> yields <0.25 0.50> < 0 0 -1> yields <0.75 0.50 >
		auto theta = acos(-p.y());
		auto phi = atan2(-p.z(), p.x()) + pi;
		u = phi / (2 * pi);
		v = theta / pi;
	}

};

// 判断光线与球体碰撞
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec)
const {
	vec3 oc = r.origin() - center;
	auto a = r.direction().length_squared();
	auto half_b = dot(oc, r.direction());
	auto c = oc.length_squared() - radius * radius;
	auto discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
	auto sqrtd = sqrt(discriminant);
	// Find the nearest root that lies in the acceptable range.
	auto root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root) {
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
			return false;
	}
	rec.t = root;
	rec.p = r.at(rec.t);
	vec3 outward_normal = (rec.p - center) / radius;
	rec.set_face_normal(r, outward_normal);
	get_sphere_uv(outward_normal, rec.u, rec.v);
	rec.mat_ptr = mat_ptr;

	return true;
}

// 球体边界框
bool sphere::bounding_box(double time0, double time1, aabb& output_box) const {
	output_box = aabb(
		center - vec3(radius, radius, radius),
		center + vec3(radius, radius, radius));
	return true;
}

// 相机类
class camera {
public:
	camera(
		point3 lookfrom,
		point3 lookat,
		vec3 vup,
		double vfov, // vertical field-of-view in degrees
		double aspect_ratio,
		double aperture, // 散焦模糊参数
		double focus_dist,
		double _time0 = 0,
		double _time1 = 0
	) {
		auto theta = degrees_to_radians(vfov);
		auto h = tan(theta / 2);
		auto viewport_height = 2.0 * h;
		auto viewport_width = aspect_ratio * viewport_height;

		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		origin = lookfrom;
		horizontal = focus_dist * viewport_width * u;
		vertical = focus_dist * viewport_height * v;
		lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w; // 散焦模糊中心
		lens_radius = aperture / 2; // 散焦模糊半径

		time0 = _time0;
		time1 = _time1;
	}
	ray get_ray(double s, double t) const {
		vec3 rd = lens_radius * random_in_unit_disk();
		vec3 offset = u * rd.x() + v * rd.y();
		return ray(
			origin + offset,
			lower_left_corner + s * horizontal + t * vertical - origin -
			offset,
			random_double(time0, time1)
		);
	}
private:
	point3 origin;
	point3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 u, v, w;
	double lens_radius;
	double time0, time1; // shutter open/close times
};


// World 一些场景
// 1、小球场景
hittable_list random_scene() {
	hittable_list world;
	// 地板
	auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
	world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));
	//棋盘格纹理地板
	/*auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1),
		color(0.9, 0.9, 0.9));
	world.add(make_shared<sphere>(point3(0, -1000, 0), 1000,
		make_shared<lambertian>(checker)));*/
	//// 大理石纹理
	//auto pertext = make_shared<noise_texture>(4); 
	//world.add(make_shared<sphere>(point3(0, -1000, 0), 1000,
	//	make_shared<lambertian>(pertext)));
	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			auto choose_mat = random_double();
			point3 center(a + 0.9 * random_double(), 0.2, b +
				0.9 * random_double());
			if ((center - point3(4, 0.2, 0)).length() > 0.9) {
				shared_ptr<material> sphere_material;
				if (choose_mat < 0.8) {
					// diffuse
					auto albedo = color::random() * color::random();
					sphere_material = make_shared<lambertian>(albedo);
					world.add(make_shared<sphere>(center, 0.2,
						sphere_material));
				}
				else if (choose_mat < 0.95) {
					// metal
					auto albedo = color::random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<sphere>(center, 0.2,
						sphere_material));
				}
				else {
					// glass
					sphere_material = make_shared<dielectric>(1.5);
					world.add(make_shared<sphere>(center, 0.2,
						sphere_material));
				}
			}
		}
	}
	auto material1 = make_shared<dielectric>(1.5);
	world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));
	 auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
	world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));
	//// 贴图
	//auto star_texture = make_shared<image_texture>("../image/star.jpeg");
	//auto star_surface = make_shared<lambertian>(star_texture);
	//world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, star_surface));

	 auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
	 world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));
	/*auto earth_texture = make_shared<image_texture>("../image/earthmap.jpg");
	auto earth_surface = make_shared<lambertian>(earth_texture);
	world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, earth_surface));*/

	return world;
}

// 两个球场景: 棋盘格纹理
hittable_list two_spheres() {
	hittable_list objects;
	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1),
		color(0.9, 0.9, 0.9));
	objects.add(make_shared<sphere>(point3(0, -10, 0), 10,
		make_shared<lambertian>(checker)));
	objects.add(make_shared<sphere>(point3(0, 10, 0), 10,
		make_shared<lambertian>(checker)));
	return objects;
}

// 球+平面场景: perlin噪声纹理
hittable_list two_perlin_spheres() {
	hittable_list objects;
	auto pertext = make_shared<noise_texture>(4); // 大理石纹理
	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000,
		make_shared<lambertian>(pertext)));
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>
		(pertext)));
	return objects;
}

// 4、图片纹理场景，把图片贴在球上: 地球
hittable_list earth() {
	auto earth_texture = make_shared<image_texture>("../image/earthmap.jpg");
	auto earth_surface = make_shared<lambertian>(earth_texture);
	auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);
	return hittable_list(globe);
}

// 5、点光源场景
hittable_list simple_light() {
	hittable_list objects;

	// 添加一个噪声纹理的地面球
	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));

	// 添加一个噪声纹理的小球体
	objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));
	
	// 添加一个球形光源
	auto sphere_light = make_shared<diffuse_light>(color(4, 4, 4));
	objects.add(make_shared<sphere>(point3(0, 7, 0), 2, sphere_light));

	// 创建一个点光源
	auto light = make_shared<point_light>(point3(2, 5, 2), color(4, 4, 4));
	objects.add(make_shared<sphere>(point3(2, 5, 2), 0.05, light));
	return objects;
}

// 6、康奈尔盒子
hittable_list cornell_box() {
	hittable_list objects;
	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(15, 15, 15));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
	// 加两个立方体
	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330,
		165), white);
	box1 = make_shared<rotate_y>(box1, 15); // 旋转
	box1 = make_shared<translate>(box1, vec3(265, 0, 295)); // 平移
	objects.add(box1);
	shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0),
		point3(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, vec3(130, 0, 65));
	objects.add(box2);


	return objects;
}

// 7、烟雾+康奈尔盒子
hittable_list cornell_smoke() {
	hittable_list objects;
	auto red = make_shared<lambertian>(color(.65, .05, .05));
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	auto green = make_shared<lambertian>(color(.12, .45, .15));
	auto light = make_shared<diffuse_light>(color(7, 7, 7));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
	shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0),
		point3(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, vec3(265, 0, 295));
	shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0),
		point3(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, vec3(130, 0, 65));
	objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
	objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));
	return objects;
}

// 8、最终场景
hittable_list final_scene() {
	hittable_list boxes1;
	auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));
	const int boxes_per_side = 20;
	for (int i = 0; i < boxes_per_side; i++) {
		for (int j = 0; j < boxes_per_side; j++) {
			auto w = 100.0;
			auto x0 = -1000.0 + i * w;
			auto z0 = -1000.0 + j * w;
			auto y0 = 0.0;
			auto x1 = x0 + w;
			auto y1 = random_double(1, 101);
			auto z1 = z0 + w;
			boxes1.add(make_shared<box>(point3(x0, y0, z0), point3(x1, y1, z1),
				ground));
		}
	}
	hittable_list objects;
	objects.add(make_shared<bvh_node>(boxes1, 0, 1));
	auto light = make_shared<diffuse_light>(color(7, 7, 7));
	objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));
	auto center1 = point3(400, 400, 200);
	auto center2 = center1 + vec3(30, 0, 0);
	auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3,
		0.1));
	objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50,
		moving_sphere_material));
	objects.add(make_shared<sphere>(point3(260, 150, 45), 50,
		make_shared<dielectric>(1.5)));
	objects.add(make_shared<sphere>(
		point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)
	));
	auto boundary = make_shared<sphere>(point3(360, 150, 145), 70,
		make_shared<dielectric>(1.5));
	objects.add(boundary);
	objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4,
		0.9)));
	boundary = make_shared<sphere>(point3(0, 0, 0), 5000,
		make_shared<dielectric>(1.5));
	objects.add(make_shared<constant_medium>(boundary, .0001, color(1, 1, 1)));
	auto emat = make_shared<lambertian>(make_shared<image_texture>
		("../image/earthmap.jpg"));
	objects.add(make_shared<sphere>(point3(400, 200, 400), 100, emat));
	auto pertext = make_shared<noise_texture>(0.1);
	objects.add(make_shared<sphere>(point3(220, 280, 300), 80,
		make_shared<lambertian>(pertext)));
	hittable_list boxes2;
	auto white = make_shared<lambertian>(color(.73, .73, .73));
	int ns = 1000;
	for (int j = 0; j < ns; j++) {
		boxes2.add(make_shared<sphere>(point3::random(0, 165), 10, white));
	}
	objects.add(make_shared<translate>(
		make_shared<rotate_y>(
			make_shared<bvh_node>(boxes2, 0.0, 1.0), 15),
		vec3(-100, 270, 395)
	)
	);
	return objects;
}


// main 函数
int main(int argc, char* args[])
{
	// Create window app handle
	WindowsApp::ptr winApp = WindowsApp::getInstance(gWidth, gHeight, "Group 19: Ray Tracing");
	if (winApp == nullptr)
	{
		std::cerr << "Error: failed to create a window handler" << std::endl;
		return -1;
	}

	// Memory allocation for canvas
	gCanvas.resize(gHeight, std::vector<color>(gWidth));

	// Launch the rendering thread
	// Note: we run the rendering task in another thread to avoid GUI blocking
	std::thread renderingThread(rendering);

	// Window app loop
	while (!winApp->shouldWindowClose())
	{
		// Process event
		winApp->processEvent();

		// Display to the screen
		winApp->updateScreenSurface(gCanvas);

	}

	renderingThread.join();

	// Save the rendered image to a file
    winApp->saveCanvasToImage(gCanvas, save_path);

	return 0;
}

// 每个像素点赋值颜色
void write_color(int x, int y, color pixel_color)
{
	// Out-of-range detection
	if (x < 0 || x >= gWidth)
	{
		std::cerr << "Warnning: try to write the pixel out of range: (x,y) -> (" << x << "," << y << ")" << std::endl;
		return;
	}

	if (y < 0 || y >= gHeight)
	{
		std::cerr << "Warnning: try to write the pixel out of range: (x,y) -> (" << x << "," << y << ")" << std::endl;
		return;
	}

	// Note: x -> the column number, y -> the row number
	gCanvas[y][x] = pixel_color;

}

void write_color(int x, int y, color pixel_color, int samples_per_pixel) {
	// Out-of-range detection
	if (x < 0 || x >= gWidth)
	{
		std::cerr << "Warnning: try to write the pixel out of range: (x,y) -> (" << x << "," << y << ")" << std::endl;
		return;
	}

	if (y < 0 || y >= gHeight)
	{
		std::cerr << "Warnning: try to write the pixel out of range: (x,y) -> (" << x << "," << y << ")" << std::endl;
		return;
	}

	/*gCanvas[y][x] = pixel_color / samples_per_pixel;*/

	auto r = pixel_color.x();
	auto g = pixel_color.y();
	auto b = pixel_color.z();
	// Divide the color by the number of samples and gamma-correct for gamma = 2.0.
	auto scale = 1.0 / samples_per_pixel;
	r = sqrt(scale * r);
	g = sqrt(scale * g);
	b = sqrt(scale * b);

	gCanvas[y][x] = color(r, g, b);
}

vec3 normalize(const vec3& v) {
	double len = v.length();
	if (len > 0) {
		return v / len;
	}
	return v; // 如果向量长度为 0，返回原向量
}

// 平行光源类
class directional_light {
public:
	directional_light(const vec3& direction, const color& light_color)
		: direction(normalize(direction)), light_color(light_color) {}

	// 获取光源方向
	vec3 get_direction() const {
		return direction;
	}

	// 获取光源颜色
	color get_color() const {
		return light_color;
	}

private:
	vec3 direction;    // 光源方向
	color light_color; // 光源颜色
};


// 射线上渲染颜色
// 方向性光源
color ray_color(const ray& r, const color& background, const hittable& world, const directional_light& light,
	int depth) {
	hit_record rec;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0)
		return color(0, 0, 0);
	// If the ray hits nothing, return the background color.
	if (!world.hit(r, 0.001, infinity, rec))
		return background;
	ray scattered;
	color attenuation;
	color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

	// 处理平行光源
	vec3 light_direction = light.get_direction();
	color light_color = light.get_color();

	// 发射阴影光线，检查是否被遮挡
	ray shadow_ray(rec.p, -light_direction);
	hit_record shadow_rec;
	if (!world.hit(shadow_ray, 0.001, infinity, shadow_rec)) {
		// 如果没有被遮挡，计算光照
		double light_intensity = dot(rec.normal, -light_direction);
		if (light_intensity > 0) {
			emitted += light_color * light_intensity;
		}
}

	if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		return emitted;
	return emitted + attenuation * ray_color(scattered, background, world, light,
		depth - 1);
}

// 点光源
color ray_color1(const ray& r, const color& background, const hittable& world, const point_light& light, int depth) {
	hit_record rec;
	if (depth <= 0)
		return color(0, 0, 0);
	if (!world.hit(r, 0.001, infinity, rec))
		return background;

	ray scattered;
	color attenuation;
	color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

	if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		return emitted;

	color light_color = color(0, 0, 0);
	if (!light.in_shadow(rec.p, world)) {
		double d = light.distance(rec.p);
		light_color = light.color / (d * d);
	}

	return emitted + attenuation * ray_color1(scattered, background, world, light, depth - 1) + light_color;
}
// 渲染
// Render

// 渲染单行
void render_row(int j, int image_width, int image_height, int samples_per_pixel, const camera& cam, const color& background, const hittable& world, const point_light& light, int max_depth) {
    for (int i = 0; i < image_width; i++) {
        color pixel_color(0, 0, 0);
        for (int s = 0; s < samples_per_pixel; ++s) {
            auto u = (double(i) + random_double()) / (image_width - 1);
            auto v = (double(j) + random_double()) / (image_height - 1);
            ray r = cam.get_ray(u, v);
            pixel_color += ray_color1(r, background, world, light, max_depth); // 点光源
			//pixel_color += ray_color(r, background, world, max_depth); // 蒙特卡洛积分
			/*pixel_color += ray_color(r, world, max_depth);*/
		}
        write_color(i, j, pixel_color, samples_per_pixel);
    }
}

// 渲染全屏幕
void render_image(int image_width, int image_height, int samples_per_pixel, const camera& cam, const color& background, const hittable& world, const point_light& light, int max_depth) {
    std::vector<std::thread> threads;// 多线程

    for (int j = image_height - 1; j >= 0; j--) {
		// 创建线程，每个线程渲染一行
        threads.push_back(std::thread(render_row, j, image_width, image_height, samples_per_pixel, std::ref(cam), std::ref(background), std::ref(world), std::ref(light), max_depth));
    }
	// 等待所有线程结束
    for (auto& t : threads) {
        t.join();
    }
}

// 场景、相机选择与渲染
void rendering()
{
	double startFrame = clock();

	printf("CGAssignment4 (built %s at %s) \n", __DATE__, __TIME__);
	std::cout << "Ray-tracing based rendering launched..." << std::endl;

	// Image

	const int image_width = gWidth;
	const int image_height = gHeight;

	// World
	hittable_list world;
	point3 lookfrom;
	point3 lookat;
	// 添加一个平行光源
	const point_light light = point_light(point3(2, 5, 2), color(4, 4, 4));

	auto vfov = 40.0;
	auto aperture = 0.0;
	// 选择场景
	switch (8) {
	case 1:
		// 很多小球场景
		world = random_scene();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		aperture = 0.1;
		break;
	case 2:
		// 两个球场景
		world = two_spheres();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;
	case 3:
		// perlin噪声场景
		world = two_perlin_spheres();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;
	case 4:
		world = earth();
		lookfrom = point3(13, 2, 3);
		lookat = point3(0, 0, 0);
		vfov = 20.0;
		break;
	case 5: // 点光源
		background = color(0, 0, 0); // 设置背景颜色是黑色
		world = simple_light();
		lookfrom = point3(26, 3, 6);
		lookat = point3(0, 2, 0);
		vfov = 20.0;
		break;
	case 6: // 康奈尔盒子
		world = cornell_box();
		background = color(0, 0, 0);
		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		break;
	case 7: // 烟雾康奈尔盒子
		world = cornell_smoke();
		lookfrom = point3(278, 278, -800);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		break;
	default:
	case 8: // 最终场景
		world = final_scene();
		background = color(0, 0, 0);
		lookfrom = point3(478, 278, -600);
		lookat = point3(278, 278, 0);
		vfov = 40.0;
		break;
	
	}
	
	
	// Camera
	vec3 vup(0, 1, 0);
	auto dist_to_focus = 10.0;

	camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus,
		0.0, 1.0);

	// Render
	render_image(image_width, image_height, samples_per_pixel, cam, background, world, light, max_depth);

	double endFrame = clock();
	double timeConsuming = static_cast<double>(endFrame - startFrame) / CLOCKS_PER_SEC;
	std::cout << "Ray-tracing based rendering over..." << std::endl;
	std::cout << "The rendering task took " << timeConsuming << " seconds" << std::endl;
}