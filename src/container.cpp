#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud altérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
	printf("intersect");
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
	//printf("intersect");

    // Acceder a la lista de objetos
    bool hit_anything = false;
    double closest_so_far = t_max;

    // Iterar sobre la lista de objetos en la escena
    //std::cout << objects[0]->local_intersect(ray, t_min, t_max, hit) << "\n";
    for (auto obj : objects) {
        hit_anything = obj->intersect(ray, t_min, t_max, hit);
        if(hit_anything){
			//std::cout << hit;
		}
       
    }

    return hit_anything;  
	
}
