#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	//printf("intersect");
	return true;
}

// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
// 		- Détecter l'intersection avec l'AABB
//			- Si intersection, détecter l'intersection avec la géométrie.
//				- Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	//printf(Naive.IContainer.aabb);

    bool hit_anything = false;
    double closest_so_far = t_max;

    for (auto obj : objects) {
        // Temporales para mantener los datos de intersección actuales
        Intersection temp_hit;

        // Intentamos intersectar el rayo con el objeto actual
        if (obj->intersect(ray, t_min, closest_so_far, &temp_hit)) {
            hit_anything = true;

            // Verificar si esta intersección está más cerca que las anteriores
            if (temp_hit.depth < closest_so_far) {
                closest_so_far = temp_hit.depth;  // Actualizar el z-buffer
                *hit = temp_hit;  // Actualizar la intersección más cercana
            }
        }
    }

    return hit_anything;
	
}
