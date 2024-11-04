#include "aabb.h" 

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.
bool AABB::intersect(Ray ray, double t_min, double t_max)  {
    for(int i=0; i<3; i++) {
		if (ray.direction[i] == 0.0) {
            if (ray.origin[i] < min[i] || ray.origin[i] > max[i]) return false;
            continue;
        }

		double t0 = (min[i] - ray.origin[i])/ray.direction[i];
		double t1 = (max[i] - ray.origin[i])/ray.direction[i];

		if(ray.direction[i]<0.0) std::swap(t0, t1);

		t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min) return false;
	}

	return true;
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
	return {
        {aabb.min.x, aabb.min.y, aabb.min.z},
        {aabb.min.x, aabb.min.y, aabb.max.z},
        {aabb.min.x, aabb.max.y, aabb.min.z},
        {aabb.min.x, aabb.max.y, aabb.max.z},
        {aabb.max.x, aabb.min.y, aabb.min.z},
        {aabb.max.x, aabb.min.y, aabb.max.z},
        {aabb.max.x, aabb.max.y, aabb.min.z},
        {aabb.max.x, aabb.max.y, aabb.max.z}
    };
};

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.
AABB construct_aabb(std::vector<double3> points) {
	// Si aucun point n'est fourni, retourner un AABB infini.
	if(points.empty()) {
		return AABB{double3{-DBL_MAX,-DBL_MAX,-DBL_MAX},double3{DBL_MAX,DBL_MAX,DBL_MAX}};
	}

	double3 min_point = {DBL_MAX, DBL_MAX, DBL_MAX};
    double3 max_point = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

	// Parcourir les points pour déterminer les limites de l'AABB.
    for(const auto& point:points) {
        min_point.x = std::min(min_point.x, point.x);
        min_point.y = std::min(min_point.y, point.y);
        min_point.z = std::min(min_point.z, point.z);

        max_point.x = std::max(max_point.x, point.x);
        max_point.y = std::max(max_point.y, point.y);
        max_point.z = std::max(max_point.z, point.z);
    }

    return AABB{min_point, max_point};
	
};

AABB combine(AABB a, AABB b) {
	return AABB{min(a.min,b.min),max(a.max,b.max)};
};

bool compare(AABB a, AABB b, int axis){
	return a.min[axis] < b.min[axis];
};