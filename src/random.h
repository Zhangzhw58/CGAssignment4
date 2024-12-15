#ifndef RANDOM_H
#define RANDOM_H

// ���������
#include <random>

inline double random_double() {
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;
	return distribution(generator);
}

// ����[min, max)�������
inline double random_double(double min, double max) {
	// Returns a random real in [min,max).
	return min + (max - min) * random_double();
}

// ����[min, max)�������(int����)
inline int random_int(int min, int max) {
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}

#endif