#ifndef HITTABLE_H
#define HITTABLE_H
#include "ray.h"

#include "utils.h"
#include "aabb.h"

class material;

// 击中记录
struct hit_record {
	point3 p;
	vec3 normal;
	shared_ptr<material> mat_ptr;
	double t;
	double u;
	double v;
	bool front_face; // 正面
	inline void set_face_normal(const ray& r, const vec3& outward_normal) {
		// 设置法线
		front_face = dot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

// 可击中类
class hittable {
public:
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
	virtual bool bounding_box(double time0, double time1, aabb& output_box)const = 0;
};

#endif