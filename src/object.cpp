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
    // El vector L es desde el centro de la esfera (que está en (0, 0, 0)) al origen del rayo
    //double3 L = ray.origin - this->center; // Si está en el origen, esto podría ser simplemente ray.origin.
	double3 L = ray.origin; // Si está en el origen, esto podría ser simplemente ray.origin.
    
    double a = dot(ray.direction, ray.direction); // Generalmente es 1 si el rayo está normalizado.
    double b = 2 * dot(L, ray.direction); 
    double c = dot(L, L) - pow(this->radius, 2); 
    
    double discriminant = pow(b, 2) - 4 * a * c;

    // Si el discriminante es negativo, no hay intersección
    if (discriminant < 0) return false;

    // Calcula las dos posibles soluciones (t0 y t1)
    double sqrt_discriminant = sqrt(discriminant);
    double t0 = (-b - sqrt_discriminant) / (2 * a); 
    double t1 = (-b + sqrt_discriminant) / (2 * a);

    // Ordena t0 y t1 para que t0 sea el menor
    if (t0 > t1) std::swap(t0, t1);

    // Chequea la intersección más cercana que esté dentro del rango
    if (t0 < t_min || t0 > t_max) {
        // Si t0 no es válido, verifica t1
        if (t1 < t_min || t1 > t_max) {
            return false; // No hay intersección válida
        }
        t0 = t1; // Usa t1 si es la única intersección válida
    }

    // Actualiza los datos de la intersección
    hit->depth = t0;
    hit->position = ray.origin + t0 * ray.direction; // Posición de la intersección
	//hit->normal = normalize(hit->position - this->center); // Vector normal (normalizado)
    hit->normal = normalize(hit->position); // Vector normal (normalizado)

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
	return false;
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
    return false;
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
	return false;
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

	return false;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
	return Object::compute_aabb();
}