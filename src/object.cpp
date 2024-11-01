#include "object.h"
#include "linalg/linalg.h"
using namespace linalg::aliases;

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
	return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray, double t_min, double t_max, Intersection *hit) 
{
    //printf("Sphere");
    // Vector L is from the center of the sphere (which is at (0, 0, 0)) to the origin of the ray
    //double3 L = ray.origin - this->center; // If you are at the origin, this could simply be ray.origin.
	double3 L = ray.origin; 
    
    // It is usually 1 if the beam is normalized.
    double a = dot(ray.direction, ray.direction); 
    double b = 2 * dot(L, ray.direction); 
    double c = dot(L, L) - pow(this->radius, 2); 
    
    double discriminant = pow(b, 2) - 4 * a * c;

    if (discriminant < 0) return false;

    double sqrt_discriminant = sqrt(discriminant);
    double t0 = (-b - sqrt_discriminant) / (2 * a); 
    double t1 = (-b + sqrt_discriminant) / (2 * a);

    // Arrange t0 and t1 so that t0 is the smallest
    if (t0 > t1) std::swap(t0, t1);

    // Check the nearest intersection that is within range
    if (t0 < t_min || t0 > t_max) {
        // If t0 is invalid, check t1
        if (t1 < t_min || t1 > t_max) {
            return false; 
        }
        t0 = t1; // Use t1 if it is the only valid intersection
    }

    // Updates the intersection data
    hit->depth = t0;
    hit->position = ray.origin + t0 * ray.direction; // Position of the intersection
	//hit->normal = normalize(hit->position - this->center); // Vector normal (normalizado)
    hit->normal = normalize(hit->position); 

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
	return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray ray, 
							double t_min, double t_max, 
							Intersection *hit)
{
    //printf("Quad");
	// The normal of the plane is in the Z+ direction (0, 0, 1)
    double3 normal(0, 0, 1);
    
    // The distance from the origin of the square to the ray (assuming the square is at Z=0)
    double d = dot(normal, double3(0, 0, 0));

    // Calculate the denominator for the equation of the plane
    double denominator = dot(normal, ray.direction);
    
    // If the ray is parallel to the plane (denominator is zero), there is no intersection
    if (fabs(denominator) < 1e-6) return false;

    // Calculate the value of t at the intersection
    double t = (d - dot(normal, ray.origin)) / denominator;

    // Check if t is within the range [t_min, t_max]
    if (t < t_min || t > t_max) return false;

    // Calculate the position of the intersection
    double3 intersection_point = ray.origin + t * ray.direction;

    // Check if the intersection position is within the boundaries of the square
    if (intersection_point.x < -half_size || intersection_point.x > half_size ||
        intersection_point.y < -half_size || intersection_point.y > half_size) {
        return false; 
    }

    hit->depth = t;
    hit->position = intersection_point; 
    hit->normal = normal; 

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
	return Object::compute_aabb();
	//return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
//
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Cylinder::local_intersect(Ray ray, 
							   double t_min, double t_max, 
							   Intersection *hit)
{
    // Transform the ray to the local space of the cylinder
    double3 origin = ray.origin;
    double3 direction = ray.direction;

    double a = direction.x * direction.x + direction.z * direction.z;
    double b = 2 * (origin.x * direction.x + origin.z * direction.z);
    double c = origin.x * origin.x + origin.z * origin.z - radius * radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return false;

    // Calculate the two possible solutions (t0 and t1)
    double sqrt_discriminant = sqrt(discriminant);
    double t0 = (-b - sqrt_discriminant) / (2 * a);
    double t1 = (-b + sqrt_discriminant) / (2 * a);

    if (t0 > t1) std::swap(t0, t1);

    // Check the nearest intersection that is within range
    double y0 = origin.y + t0 * direction.y;
    double y1 = origin.y + t1 * direction.y;

    if ((y0 < -half_height || y0 > half_height) && (y1 < -half_height || y1 > half_height)) {
        return false; 
    }

    // Find the first value of t that is within the height limits
    if (y0 >= -half_height && y0 <= half_height) {
        // t0 is within the height of the cylinder
        hit->depth = t0;
    } else {
        // t1 is within the height of the cylinder
        hit->depth = t1;
    }

    // Calculates the position of the intersection
    hit->position = ray.origin + hit->depth * ray.direction;

    // Calculate the normal at the intersection point
    // Project the intersection point onto the plane of the cylinder
    hit->normal = normalize(double3(hit->position.x, 0, hit->position.z)); // Normal en la superficie

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
	return Object::compute_aabb();
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
//
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
//
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
//
bool Mesh::local_intersect(Ray ray,  
						   double t_min, double t_max, 
						   Intersection* hit)
{
	bool intersected = false; 

    for (const Triangle& tri : triangles) {
        Intersection tempHit;
        
        // We try to intersect the ray with the current triangle
        if (intersect_triangle(ray, t_min, t_max, tri, &tempHit)) {
            // If the intersection is closer than any other found before
            if (tempHit.depth < hit->depth) {
                *hit = tempHit; // We updated the information of the nearest intersection
                t_max = tempHit.depth;
                intersected = true;
            }
        }
    }

    return intersected; 
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray  ray, 
							  double t_min, double t_max,
							  Triangle const tri,
							  Intersection *hit)
{
	// Extrait chaque position de sommet des données du maillage.
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A (Pour faciliter les explications)
	double3 const &p1 = positions[tri[1].pi]; // ou Sommet B
	double3 const &p2 = positions[tri[2].pi]; // ou Sommet C

	// Triangle en question. Respectez la convention suivante pour vos variables.
	//
	//     A
	//    / \
	//   /   \
	//  B --> C
	//
	// Respectez la règle de la main droite pour la normale.

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Pour plus de d'informations sur la géométrie, référez-vous à la classe dans object.hpp.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.

    // Calculate the vectors of the triangle
    double3 edge1 = p1 - p0; 
    double3 edge2 = p2 - p0; 

    double3 h = cross(ray.direction, edge2);
    double a = dot(edge1, h);

    // Check if the ray is parallel to the triangle
    if (fabs(a) < 1e-6) return false; 

    double f = 1.0 / a;
    double3 s = ray.origin - p0;
    double u = f * dot(s, h);

    // Check if the intersection is outside the triangle
    if (u < 0.0 || u > 1.0) return false;

    double3 q = cross(s, edge1);
    double v = f * dot(ray.direction, q);

    // Check if the intersection is outside the triangle
    if (v < 0.0 || u + v > 1.0) return false;

    //Calculate t to determine if the ray intersects within the t_min and t_max limits
    double t = f * dot(edge2, q);
    
    // Check if the intersection is within the allowed range
    if (t < t_min || t > t_max) return false;

    // If there is intersection, fill the hit structure
    hit->depth = t;
    hit->position = ray.origin + t * ray.direction; 

    // Calculate the normal of the triangle
    hit->normal = normalize(cross(edge1, edge2)); // Normal del triángulo

    // Opcional: calcular coordenadas de textura, si son necesarias
    // Esto puede variar dependiendo de cómo almacenes las coordenadas UV.
    // hit->tex_coords = ...; // Asignar coordenadas de textura si es necesario.

    return true; 
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	return Object::compute_aabb();
}