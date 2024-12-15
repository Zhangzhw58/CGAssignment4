#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <cstdlib>
#include "random.h"

// Usings
using std::shared_ptr;
using std::make_shared;
using std::sqrt;
// Constants
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;
// Utility Functions
inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

// 裁剪 x，使超范围的x值变成[min 或 max]
inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

#include "ray.h"
#include "vec3.h"

#endif