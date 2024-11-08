#include "container.h"

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
// 		- S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//		- Sinon, il s'agit d'un noeud intérieur.
//			- Faites l'intersection du rayon avec le AABB gauche et droite. 
//				- S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool hit_anything = false;
    double closest_so_far = t_max; // Initialise la profondeur maximale avec la limite de distance.
    std::vector<BVHNode*> node_stack; // Utilise un vector pour simuler une pile.
    node_stack.push_back(root);

    // Parcours en profondeur de l'arbre BVH.
    while (!node_stack.empty()) {
        BVHNode* node = node_stack.back();
        node_stack.pop_back();

        // Vérifie l'intersection avec l'AABB du nœud courant.
        if (!node->aabb.intersect(ray, t_min, closest_so_far)) {
            continue;  // Si aucune intersection avec cet AABB, passer au prochain nœud
        }

        // Si c'est une feuille, teste l'intersection avec l'objet.
        if (node->left == nullptr && node->right == nullptr) {
            Intersection temp_hit;
            if (objects[node->idx]->intersect(ray, t_min, closest_so_far, &temp_hit)) {
                if (temp_hit.depth < closest_so_far) {
                    closest_so_far = temp_hit.depth;
                    *hit = temp_hit;
                    hit_anything = true;
                }
            }
        } else {
            // Si est un nœud interne, on empile les enfants (gauche et droite) pour les explorer plus tard.
            if (node->left) node_stack.push_back(node->left);
            if (node->right) node_stack.push_back(node->right);
        }
    }

    // Renvoie true si une intersection a été trouvée avec au moins un objet, sinon false.
    return hit_anything;
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
    double closest_so_far = t_max; // Initialise la profondeur maximale avec la limite de distance.

    // Parcourt chaque objet dans la liste.
    for (size_t i = 0; i < objects.size(); i++) {
        // Vérifie d'abord l'intersection avec l'AABB de l'objet.
        if (aabbs[i].intersect(ray, t_min, closest_so_far)) {
            // Effectue un test d'intersection plus précis avec la géométrie de l'objet.
            Intersection temp_hit;
            if (objects[i]->intersect(ray, t_min, closest_so_far, &temp_hit)) {
                hit_anything = true;

                // Met a jour les informations d'intersection si l'objet est plus proche
                if (temp_hit.depth < closest_so_far) {
                    closest_so_far = temp_hit.depth;  // Met à jour la profondeur maximale la plus proche.
                    *hit = temp_hit;  // Met à jour l'intersection la plus proche.
                }
            }
        }
    }

    // Renvoie true si une intersection a été trouvée avec au moins un objet, sinon false.
    return hit_anything;
	
}
