#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "vec3.h"
#include"glm/glm.hpp"
#include "hittable.h"
#include <array>

class vertex : public point3
{
public:
    vertex(){}
    vertex(double x, double y, double z , glm::vec3 norm , glm::vec3 c , glm::vec2 tx):tex(tx) { 
        e[0] = x; 
        e[1] = y; 
        e[2] = z;
        normal = { norm.x , norm.y , norm.z };
        col = { c.x , c.y , c.z };
    }
    vec3 normal = {0, 0, 0};
    vec3 col = { 0 , 0 , 0 };
    glm::vec2 tex = { 0 , 0 };

    double u = tex.x;
    double v = tex.y;
};

class triangle : public hittable{
public:
    triangle(const vertex &_v0, const vertex &_v1, const vertex &_v2, shared_ptr<material> m) : v0(_v0), v1(_v1), v2(_v2), mat_ptr(m)
    {
        // 计算三角形的法线
        face_normal = unit_vector(cross(v1 - v0, v2 - v0));
        v0v1=v1-v0;
        v0v2=v2-v0;
    }
    // 计算光线与三角形的相交
    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const
    {
        // 计算光线起点到三角形顶点的向量
        vec3 tvec = r.origin() - v0;
        // 计算光线方向与三角形边的叉积
        vec3 pvec = cross(r.direction(), v0v2);
        // 计算行列式
        auto det = dot(v0v1, pvec);
        // 如果行列式接近于零，说明光线与三角形平行，返回 false
        if (fabs(det) < 1e-5) return false;
        auto invDet = 1.0 / det;
        // 计算重心坐标 u
        double u = dot(tvec, pvec) * invDet;
        // 如果 u 不在 [0, 1] 范围内，返回 false
        if (u < 0 || u > 1) return false;
        // 计算光线方向与三角形另一边的叉积
        vec3 qvec = cross(tvec, v0v1);
        // 计算重心坐标 v
        double v = dot(r.direction(), qvec) * invDet;
        // 如果 v 不在 [0, 1] 范围内，或者 u + v > 1，返回 false
        if (v < 0 || u + v > 1) return false;
        // 计算光线参数 t
        double t = dot(v0v2, qvec) * invDet;
        // 如果 t 不在 [t_min, t_max] 范围内，返回 false
        if (t < t_min || t > t_max) return false;
        // 设置击中记录
        rec.t = t;
        rec.p = r.at(t);
        rec.u = u;
        rec.v = v;
        rec.mat_ptr = mat_ptr;
        rec.set_face_normal(r, face_normal);
        return true;
    }

    bool bounding_box(double time0, double time1, aabb& output_box) const
    {
        vec3 min_point(
            fmin(v0.x(), fmin(v1.x(), v2.x())),
            fmin(v0.y(), fmin(v1.y(), v2.y())),
            fmin(v0.z(), fmin(v1.z(), v2.z()))
        );
        vec3 max_point(
            fmax(v0.x(), fmax(v1.x(), v2.x())),
            fmax(v0.y(), fmax(v1.y(), v2.y())),
            fmax(v0.z(), fmax(v1.z(), v2.z()))
        );
        output_box = aabb(min_point, max_point);

        return true;
    }
    // 计算概率密度函数值
    double pdf_value(const point3& origin, const vec3& direction) const
    {
        hit_record rec;
        if (!this->hit(ray(origin, direction), 0.001, infinity, rec))
        return 0;

        vec3 temp = cross(v0v1, v0v2);
        double area = temp.length()/2;
        double distance_squared = rec.t*rec.t*direction.length_squared();
        double cosine = fabs(dot(direction, rec.normal)/direction.length());
        return distance_squared/(cosine*area);
    }
    // 随机采样三角形上的一点，用于光线追踪中的重要性采样
    vec3 random(const point3& origin) const
    {
        double length1 = random_double();
        double length2 = random_double();
        point3 random_point = v0 + v0v1*length1 + v0v2*length2;
        return random_point-origin;
    }
private:
    vertex v0, v1, v2; // 三角形的三个顶点
    shared_ptr<material> mat_ptr; // 材质
    vec3 face_normal; // 三角形的法线
    vec3 v0v1, v0v2; // 三角形两条边
};

#endif