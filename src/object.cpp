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
    // Calcul des coefficients de l'équation quadratique ax^2 + bx + c = 0.
    double a = dot(ray.direction, ray.direction); 
    double b = 2 * dot(ray.origin, ray.direction); 
    double c = dot(ray.origin, ray.origin) - pow(this->radius, 2); 
    
    // Calcul du discriminant pour déterminer l'existence de solutions
    double discriminant = pow(b, 2) - 4 * a * c;
    if (discriminant < 0) return false;

    // Calcul des deux solutions de l'équation quadratique (t0 et t1).
    double sqrt_discriminant = sqrt(discriminant);
    double t0 = (-b - sqrt_discriminant) / (2 * a); 
    double t1 = (-b + sqrt_discriminant) / (2 * a);

    // On s'assure que t0 est la plus petite valeur en permutant les valeurs si nécessaire.
    if (t0 > t1) std::swap(t0, t1);

    // Vérifie si l'intersection la plus proche (t0) est dans les limites [t_min, t_max].
    if (t0 < t_min || t0 > t_max) {
        // Si t0 est hors des limites, vérifie t1.
        if (t1 < t_min || t1 > t_max) {
            return false; // Aucune intersection valide.
        }
        // Utilise t1 si c'est la seule intersection valide.
        t0 = t1; 
    }

    // Met à jour les informations de l'intersection dans la structure hit
    hit->depth = t0;
    hit->position = ray.origin + t0 * ray.direction;
    hit->normal = normalize(hit->position); 

    // Calcul des coordonnées UV en utilisant les coordonnées sphériques.
    double theta = acos(hit->normal.y); 

    // On joute PI/2 à phi, pour décaler l'angle de 90 degrés pour orienter correctement la texture.
    double phi = atan2(hit->normal.z, hit->normal.x) + PI/2; 

    // Normalise phi pour qu'il reste dans la plage [-pi, pi]
    if (phi > PI) {
        phi -= 2 * PI;
    } else if (phi < -PI) {
        phi += 2 * PI;
    }

    // Normalisation des coordonnées UV dans la plage [0, 1].
    // On inverse l'axe U pour que la texture soit à l'endroit.
    hit->uv.x = 1-(phi + PI) / (2 * PI);  
    hit->uv.y = theta / PI; 

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
    // Crée un AABB local autour de la sphère, en utilisant son rayon. La boîte est centrée autour de l'origine avec 
    // des limites de -radius à +radius dans toutes les directions.
    AABB local_aabb{double3{-radius, -radius, -radius}, double3{radius, radius, radius}};
    
    std::vector<double3> global_corners;
    std::vector<double3> corners = retrieve_corners(local_aabb);
    // Récupère les huit coins du AABB local pour les transformer en coordonnées globales.
    for(const auto& corner:corners) {
        double3 global_corner = mul(transform, {corner, 1}).xyz();
        global_corners.push_back(global_corner);
    }

    // Construit un nouveau AABB global englobant tous les coins transformés,
	return construct_aabb(global_corners);
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
    // Normale du plan sur lequel se situe le quad (face orientée vers +Z dans l'espace local).
    double3 normal(0, 0, 1); 
    
    // Calcul le dénominateur de l'équation du plan pour vérifier l'angle entre le rayon et le plan.
    double denominator = dot(normal, ray.direction);
    
    // Si le rayon est parallèle au plan (dénominateur proche de zéro), il n'y a pas d'intersection.
    if (fabs(denominator) < 0) return false;

    // Calcul la valeur de t pour déterminer le point d'intersection sur le plan.
    double t = -dot(normal, ray.origin) / denominator;

    // Vérifie que t est dans la plage [t_min, t_max].
    if (t < t_min || t > t_max) return false;

    // Vérifie dans l'espace local si le point d'intersection se situe dans les limites du quad.
    double3 intersection_point = ray.origin + t * ray.direction;
    if (intersection_point.x < -half_size || intersection_point.x > half_size ||
        intersection_point.y < -half_size || intersection_point.y > half_size) {
        return false; 
    }
    
    // Inverse la normale si elle pointe dans la même direction que le rayon, pour une orientation correcte.
    if (dot(normal, ray.direction) > 0) {
        normal = -normal;
    }

    // Met à jour les informations de l'intersection dans la structure hit.
    hit->depth = t;
    hit->position = intersection_point; 
    hit->normal = normal; 

    // Le quad s'étend de -half_size à +half_size, on convertit donc cette plage en [0,1] pour calculer les 
    // coordonnées UV.
    hit->uv.x = (intersection_point.x + half_size) / (2*half_size);
    // On inverse l'axe V pour que la texture soit à l'endroit
    hit->uv.y = 1-(intersection_point.y + half_size) / (2*half_size);

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Quad::compute_aabb() {
    // Crée un AABB local pour le quad, en utilisant half_size. La boîte est centrée autour de l'origine avec des 
    // limites de -half_size à +half_size en XY, et a une épaisseur nulle (z=0)
    AABB local_aabb{double3{-half_size, -half_size, 0}, double3{half_size, half_size, 0}};
    
    std::vector<double3> global_corners;
    std::vector<double3> corners = retrieve_corners(local_aabb);
    // Récupère les huit coins du AABB local pour les transformer en coordonnées globales.
    for(const auto& corner:corners) {
        double3 global_corner = mul(transform, {corner, 1}).xyz();
        global_corners.push_back(global_corner);
    }

    // Construit un nouveau AABB global englobant tous les coins transformés,
	AABB global_aabb = construct_aabb(global_corners);

    // Ajout d'une petite marge (epsilon) pour éviter les AABB dégénérés.
    const double epsilon = 1e-6;
    if(global_aabb.min.x == global_aabb.max.x) global_aabb.max.x += epsilon;
    if(global_aabb.min.y == global_aabb.max.y) global_aabb.max.y += epsilon;
    if(global_aabb.min.z == global_aabb.max.z) global_aabb.max.z += epsilon;

    return global_aabb;
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
    // Calcul des coefficients de l'équation quadratique ax^2 + bx + c = 0.
    double a = ray.direction.x * ray.direction.x + ray.direction.z * ray.direction.z;
    double b = 2 * (ray.origin.x * ray.direction.x + ray.origin.z * ray.direction.z);
    double c = ray.origin.x * ray.origin.x + ray.origin.z * ray.origin.z - radius * radius;

    // Calcul du discriminant pour déterminer l'existence de solutions.
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;

    // Calcul des deux solutions de l'équation quadratique (t0 et t1).
    double sqrt_discriminant = sqrt(discriminant);
    double t0 = (-b - sqrt_discriminant) / (2 * a);
    double t1 = (-b + sqrt_discriminant) / (2 * a);

    // On s'assure que t0 est la plus petite valeur en permutant les valeurs si nécessaire.
    if (t0 > t1) std::swap(t0, t1);

    // Vérifie si t0 est dans les limites de hauteur du cylindre.
    double y0 = ray.origin.y + t0 * ray.direction.y;
    if (y0 < -half_height || y0 > half_height) {
        // Si t0 est hors des limites, vérifie t1.
        double y1 = ray.origin.y + t1 * ray.direction.y;
        if (y1 < -half_height || y1 > half_height) {
            return false; // Aucune intersection valide.
        }
        // Utilise t1 si c'est la seule intersection valide.
        t0 = t1;
    }

    // Vérifie que t0 est dans l'intervalle [t_min, t_max].
    if (t0 < t_min || t0 > t_max) return false;

    // Met à jour les informations d'intersection.
    hit->depth = t0;
    hit->position = ray.origin + t0 * ray.direction;
    hit->normal = normalize(double3(hit->position.x, 0, hit->position.z));

    // Calcul des coordonnées UV en utilisant les coordonnées cylindriques.
    // On joute PI à phi, pour décaler l'angle de 180 degrés pour orienter correctement la texture.
    double phi = atan2(hit->position.z, hit->position.x) + PI;

    // Normalise phi pour qu'il reste dans la plage [-pi, pi].
    if (phi > PI) {
        phi -= 2 * PI;
    } else if (phi < -PI) {
        phi += 2 * PI;
    }

    // Normalisation des coordonnées UV dans la plage [0, 1].
    // On inverse l'axe U et V pour que la texture soit à l'endroit.
    hit->uv.x = 1-(phi + PI) / (2 * PI);  
    hit->uv.y = 1-(hit->position.y + half_height) / (2 * half_height); 

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
    // Crée un AABB local autour du cylindre, en utilisant sont rayon et sa demi-hauteur. La boîte est centrée autour 
    // de l'origine avec des limites de -half_height à half_height en Y et -radius à +radius dans toutes les autres 
    // directions.
    AABB local_aabb{double3{-radius, -half_height, -radius}, double3{radius, half_height, radius}};

    std::vector<double3> global_corners;
    std::vector<double3> corners = retrieve_corners(local_aabb);
    // Récupère les huit coins du AABB local pour les transformer en coordonnées globales.
    for(const auto& corner:corners) {
        double3 global_corner = mul(transform, {corner, 1}).xyz();
        global_corners.push_back(global_corner);
    }

    // Construit un nouveau AABB global englobant tous les coins transformés,
    return construct_aabb(global_corners);
	// return Object::compute_aabb();
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

    // Parcourt tous les triangles du maillage.     
    for (const Triangle& tri : triangles) {
        Intersection tempHit;
    
        // Teste l'intersection entre le rayon et le triangle courant.
        if (intersect_triangle(ray, t_min, t_max, tri, &tempHit)) {
            // Si l'intersection est plus proche que les intersections précédentes, on met à jour les informations 
            // d'intersection pour garder la plus proche.
            if (tempHit.depth < hit->depth) {
                *hit = tempHit;
                t_max = tempHit.depth;
                intersected = true;
            }
        }
    }

    // Retourne true si une intersection a été trouvée, sinon false.
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
	double3 const &p0 = positions[tri[0].pi]; // ou Sommet A
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

     // Calcul des vecteurs de bord du triangle : de A à B et de A à C.
    double3 edge1 = p1 - p0; 
    double3 edge2 = p2 - p0; 

    // Calcule le produit vectoriel entre la direction du rayon et le deuxième bord du triangle.
    double3 h = cross(ray.direction, edge2);

    // Calcul le dénominateur de l'équation du plan pour vérifier l'angle entre le rayon et le triangle.
    double denominator = dot(edge1, h);

    // Si le rayon est parallèle au triangle (dénominateur proche de zéro), il n'y a pas d'intersection.
    if (fabs(denominator) < 1e-6) return false; 

    // Calcule l'inverse du déterminant pour éviter des divisions répétées.
    double f = 1.0 / denominator;
    double3 s = ray.origin - p0;
    double u = f * dot(s, h);

    // Vérifie si l'intersection est à l'extérieur du triangle en évaluant la coordonnée barycentrique u.
    if (u < 0.0 || u > 1.0) return false;

    double3 q = cross(s, edge1);
    double v = f * dot(ray.direction, q);

    // Vérifie si l'intersection est à l'extérieur du triangle en utilisant les valeurs de u et v.
    if (v < 0.0 || u + v > 1.0) return false;

    //Calculate t to determine if the ray intersects within the t_min and t_max limits
    double t = f * dot(edge2, q);
    
    // Vérifie que t est dans l'intervalle [t_min, t_max].
    if (t < t_min || t > t_max) return false;

    // Remplit les informations d'intersection si un point d'intersection a été trouvé.
    hit->depth = t;
    hit->position = ray.origin + t * ray.direction; 
    hit->normal = normalize(cross(edge1, edge2));

    // Calcule les coordonnées barycentriques du point d'intersection.
    double w = 1 - u - v;

    // Interpolation des coordonnées UV à l'aide des coordonnées barycentriques.
    double2 uv0 = tex_coords[tri[0].ti];
    double2 uv1 = tex_coords[tri[1].ti];
    double2 uv2 = tex_coords[tri[2].ti];
    hit->uv.x = w * uv0.x + u * uv1.x + v * uv2.x;
    hit->uv.y = w * uv0.y + u * uv1.y + v * uv2.y;

    return true; 
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
    std::vector<double3> global_positions;
    // Parcourt chaque position de sommet dans le maillage.
    for(const auto& position:positions) {
        double3 global_position = mul(transform, {position, 1}).xyz();
        global_positions.push_back(global_position);
    }

    // Construit un nouveau AABB global englobant tous les sommets transformés,
    return construct_aabb(global_positions);
	// return Object::compute_aabb();
}