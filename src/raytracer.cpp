#include "raytracer.h"

void Raytracer::render(const Scene& scene, Frame* output)
{       
    // Crée le z_buffer.
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; //Anciennement DBL_MAX. À remplacer avec la valeur de scene.camera.z_far
    }


	//---------------------------------------------------------------------------------------------------------------
	// Nous vous fournissons ci-bas du code pour une caméra orthographique très simple. 
	// Cette caméra peut être utilisée pour tester l’intersection avec la sphère.
	// Vous devez utiliser la scène de test portho.ray pour utiliser cette caméra. 
	// Notez que votre code de caméra ne doit pas être basé sur ce code de caméra. 
	// Ce code n’est là que pour prendre en compte le développement initial du test d’intersection.
	// Pour utiliser cette caméra, vous devez supprimer les commentaires qui rendent inactive cette partie du code, 
	// et mettre en commentaires la boucle d’image originale.

	// CameraOrthographic camOrth;
	// double3 uVec{ 0,1,0 };
	// double3 vVec{ 0,0,1 };
	// double y_shift = 2.0 / scene.resolution[1];
	// double x_shift = 2.0 / scene.resolution[0];
	// double3 uVec2{ 1,0,0 };
	// double3 vVec2{ 0,1,0 };

	// for (int y = 0; y < scene.resolution[1]; y++) {
	// 	if (y % 40) {
	// 		std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
	// 	}

	// 	for (int x = 0; x < scene.resolution[0]; x++) {
	// 		double3 color{ 0,0,0 };

	// 		Intersection hit;
	// 		//double3 rayOrigin = camOrth.minPosition + uVec * x_shift * x + vVec * y_shift * y;
	// 		double3 rayOrigin = camOrth.minPosition + uVec2 * x_shift * x + vVec2 * y_shift * y;
	// 		// std::cout << "Otigin " << rayOrigin.x << "," << rayOrigin.y << "," << rayOrigin.z << "\n"; 
	// 		double3 rayDirection{ 0,0,1 };
	// 		Ray ray = Ray(rayOrigin, rayDirection);
	// 		double itHits = 0;

	// 		double z_depth = scene.camera.z_far;
	// 		if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
	// 			Material& material = ResourceManager::Instance()->materials[hit.key_material];
	// 			color = material.color_albedo;
	// 			itHits = 1.0f;
	// 			//printf("intersect");
	// 		}

	// 		output->set_color_pixel(x, y, color);
	// 		output->set_depth_pixel(x, y, itHits);
	// 		//printf("intersect");
	// 	}
	// }

	//---------------------------------------------------------------------------------------------------------------


	// @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons.
	Camera camera;

	// Vecteurs de direction de la caméra : calcul des axes u, v et w pour orienter les rayons selon l'orientation de la caméra.
	double3 vVec = normalize(scene.camera.up); 								//Direction Y de la caméra
	double3 wVec = normalize(scene.camera.center - scene.camera.position); 	//Direction Z de la caméra
	double3 uVec = normalize(cross(wVec,vVec));								//Direction X de la caméra
	double3 rayOrigin = scene.camera.position;		

	double height = 2 * scene.camera.z_near * tan(deg2rad(scene.camera.fovy/2)); // Hauteur du plan d'image
	double width = height * scene.camera.aspect; // Largeur du plan d'image
	double disk_radius = 0.0;
	bool depth_of_field_enabled = scene.camera.defocus_angle > 0.0;

	// Calcul le rayon du dique de délocalisation si l'angles de défocalisation n'est pas nul.
	if(depth_of_field_enabled) {
		disk_radius = tan(deg2rad(scene.camera.defocus_angle / 2)) * scene.camera.focus_distance;
	}

	//Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
		if (y % 40){
			std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
		}

        for(int x = 0; x < scene.resolution[0]; x++) {

			int avg_z_depth = 0;
			double3 avg_ray_color{0,0,0};

			for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
				// Génère le rayon approprié pour ce pixel.
				Ray ray;
				// Initialise la profondeur de récursivité du rayon.
				int ray_depth = 0;
				// Initialize la couleur du rayon.
				double3 ray_color{0,0,0};

				// @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				// Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius. 
				//Faites la moyenne des différentes couleurs obtenues suite à la récursion.
				//printf("random : %f", rand_double());

				// Position aléatoire dans le pixel pour éviter les effets d'aliasing.
				double xDirection = (x+(rand_double()-0.5)-scene.resolution[0]/2) * (width/scene.resolution[0]);
				double yDirection = (y+(rand_double()-0.5)-scene.resolution[1]/2) * (height/scene.resolution[1]);

 				double3 primaryDirection = normalize((scene.camera.position + uVec*xDirection + vVec*yDirection + wVec*scene.camera.z_near) - scene.camera.position);
				if (depth_of_field_enabled) {
                    // Calcul une position perturbée sur le disque de délocalisation.
                    double2 random_disk_point = random_in_unit_disk() * disk_radius;
                    double3 offset = uVec * random_disk_point.x + vVec * random_disk_point.y;

                    // Déplace l'origine du rayon et réajuste la direction pour pointer vers le point focal.
                    double3 perturbed_origin = scene.camera.position + offset;
                    double3 focus_point = scene.camera.position + primaryDirection * scene.camera.focus_distance;
                    double3 perturbed_direction = normalize(focus_point - perturbed_origin);

                    // Crée un rayon avec l'origine perturbée et la direction ajustée.
                    ray = Ray(perturbed_origin, perturbed_direction);
                } else {
                    // Si la profondeur de champ est désactivée, on utilise l'origine normale de la caméra.
                    ray = Ray(rayOrigin, primaryDirection);
                }
				double z_depth = scene.camera.z_far;

				trace(scene, ray, ray_depth, &ray_color, &z_depth);
				avg_ray_color += ray_color;
				avg_z_depth += z_depth;
			}

			// On fait la moyenne des valeurs de profondeur et de couleur après l'échantillonnage.
			avg_z_depth = avg_z_depth / scene.samples_per_pixel;
			avg_ray_color = avg_ray_color / scene.samples_per_pixel;

			// Test de profondeur
			if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
				avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
				z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

				// Met à jour la couleur de l'image (et sa profondeur).
				output->set_color_pixel(x, y, avg_ray_color);
				output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / 
										(scene.camera.z_far-scene.camera.z_near));
			}
        }
    }

    delete[] z_buffer;
}

double3 refract(const double3& incident, const double3& normal, double eta) {
	// Calcule la direction de réfraction en utilisant l'équation de Snell.
    double dots = dot(normal, incident);
    double first_term = eta * dots;
	double second_term = sqrt(1-(eta*eta)*(1-(dots*dots)));
	double3 refracted_direction = normal * (first_term - second_term) - eta * incident;

    return normalize(refracted_direction); 
}

void Raytracer::trace(const Scene& scene,
                      Ray ray, int ray_depth,
                      double3* out_color, double* out_z_depth)
{
    Intersection hit;
    if (scene.container->intersect(ray, EPSILON, *out_z_depth, &hit)) {
		// Récupère le matériau de l'objet intersecté.
        Material& material = ResourceManager::Instance()->materials[hit.key_material];
        *out_color = shade(scene, hit);

        if (ray_depth < scene.max_ray_depth) {
            // Reflection
            double3 reflected_direction = normalize(ray.direction - 2 * dot(ray.direction, hit.normal) * hit.normal);
			Ray reflected_ray = Ray(hit.position + hit.normal * EPSILON, reflected_direction);
			double3 reflected_color{0, 0, 0};
			double out_z_depth_copy = *out_z_depth;
			trace(scene, reflected_ray, ray_depth + 1, &reflected_color, &out_z_depth_copy);
			*out_color += reflected_color * material.k_reflection;

            // Refraction
			if (material.refractive_index > 2) {
				//if (dot(double3(0,0,-1),ray.direction)>0) hit.depth = DBL_MAX;
				//std::cout << "dot: " << dot(double3(0,0,1),ray.direction) << "\n";
				double3 refracted_color{0, 0, 0};
			    double eta = 1 / material.refractive_index;
			    double3 refracted_direction = refract(ray.direction, hit.normal, eta);
				//double3 refracted_direction = ray.direction;
				Ray refracted_ray = Ray(hit.position - hit.normal * EPSILON, refracted_direction);
				//Ray refracted_ray = Ray(double3(0,0,-2), refracted_direction);
				trace(scene, refracted_ray, ray_depth + 1, &refracted_color, out_z_depth);
				*out_color += refracted_color * material.k_refraction;
			    
			} 
        }

        *out_z_depth = hit.depth;
    }
}

// PETITE OPTIMISATION, on calcule les vecteurs une fois 
// Main function to calculate the vectors in the plane
std::pair<double3, double3> get_vecteurs(const double3& sphere_center, const double3& plane_point) {
    // Normal vector to the plane
    double3 normal = normalize(plane_point - sphere_center);

    // "Right" vector in the plane (perpendicular to the normal and the Z axis)
    double3 arbitrary_direction = {0, 0, 1};
    double3 right_vector = normalize(cross(normal, arbitrary_direction));

    // Vector "up" in the plane (perpendicular to the normal and to the "right")
    double3 up_vector = cross(right_vector, normal);

    // Return only the right and up vectors
    return {right_vector, up_vector};
}

double facteur_lumiere(Scene scene, double3 point, SphericalLight ligth){

	//if (ligth.radius == 0) return 1;

    int misses = 0;
	auto [right_vector, up_vector] = get_vecteurs(ligth.position, point);

    for (int i = 1; i <= RAY_QTE; ++i) {
        // Generate a random point on the unit disk
        double2 random_dir = random_in_unit_disk();

		// Get point position on scene
		// double3 ray_end = move_in_plane(ligth.position, point, ligth.radius, random_dir);
		double3 right_movement = right_vector * (random_dir.x * ligth.radius);
		double3 up_movement = up_vector * (random_dir.y * ligth.radius);
		double3 ray_end = ligth.position + right_movement + up_movement;

		Ray ray_lumiere;

        // Create the ray with random direction
        ray_lumiere = Ray(point, normalize(ray_end - point));

		double ray_depth = length(ray_end - point);

        Intersection hit_info;
        if (scene.container->intersect(ray_lumiere, EPSILON, ray_depth, &hit_info)) {
			Material& material_quad = ResourceManager::Instance()->materials[hit_info.key_material];
			if(material_quad.refractive_index < 2) misses++;
			//std::cout << "miss";
		}
    }

    // Calculate the proportion of rays that miss
    double miss_ratio = static_cast<double>(misses) / RAY_QTE;
	double porcentage_lumiere = 1 - miss_ratio;

	return porcentage_lumiere;
}

int clamp(int value, int min, int max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}


// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		* Calculer la contribution des lumières dans la scène.
//			- Itérer sur toutes les lumières.
//				- Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//	          	- Inclure la contribution diffuse. (Faites attention au produit scalare. >= 0)
//   	  	- Inclure la contribution ambiante
//      * Calculer si le point est dans l'ombre
//			- Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occludé.
//				- Ne pas considérer les points plus loins que la lumière.
//			- Par la suite, intégrer la pénombre dans votre calcul
//		* Déterminer la couleur du point d'intersection.
//        	- Si texture est présente, prende la couleur à la coordonnées uv
//			- Si aucune texture, prendre la couleur associé au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
	Material& material = ResourceManager::Instance()->materials[hit.key_material]; 
	double3 base_color;

	 // Vérifie si une texture est présente en utilisant les dimensions
    if (material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) {
        // Convertit les coordonnées UV en coordonnées de la texture
        int tex_x = static_cast<int>(hit.uv.x * material.texture_albedo.width());
        int tex_y = static_cast<int>(hit.uv.y * material.texture_albedo.height());

        // Assure que les coordonnées restent dans les limites de la texture
        tex_x = clamp(tex_x, 0, material.texture_albedo.width() - 1);
        tex_y = clamp(tex_y, 0, material.texture_albedo.height() - 1);

        // Récupère la couleur du pixel de la texture
        rgb_t tex_color;
        material.texture_albedo.get_pixel(tex_x, tex_y, tex_color);

        // Convertit la couleur de [0..255] à [0..1]
        base_color = double3(tex_color.red / 255.0, tex_color.green / 255.0, tex_color.blue / 255.0);
    } else {
        // Si la texture est absente, utilise la couleur albedo du matériau
        base_color = material.color_albedo;
    }
	
	double3 ambient(0, 0, 0);
    double3 diffuse(0, 0, 0);
    double3 specular(0, 0, 0);

    // Vector pointing towards the eye
    double3 Eye = normalize(scene.camera.position - hit.position);

    // Ambient light calculation: L_aλ * k_aλ * Sλ
    ambient = scene.ambient_light * material.k_ambient * base_color;

	double penumbra;

    // Iterate through each light in the scene
    for (const auto& light : scene.lights) {

		double3 lightDir = normalize(light.position - hit.position);

		// ================= CODE OMBRE SANS PENOMBRE ====================

        // // Light direction and distance to the light source
        
        // double lightDistance = length(light.position - hit.position);

        // // Shadow ray: emitted from the hit position towards the light source
        // // Offset the starting position to avoid "surface acne" issues
        // Ray shadowRay(hit.position + EPSILON * hit.normal, lightDir);

        // // Check if any object obstructs the shadow ray up to the light source
        // bool in_shadow = false;
        // // for (const auto& obj : scene.objects) {
        // //     // If the shadow ray intersects an object within the light distance, the point is in shadow
        // //     if (obj->intersect(shadowRay, 1e-6, lightDistance)) {
        // //         in_shadow = true;
        // //         break;
        // //     }
        // // }
		// if (scene.container->intersect(shadowRay, EPSILON, lightDistance, &hit)) in_shadow = true;

        // // If the point is in shadow, skip diffuse and specular contributions for this light
        // if (in_shadow) {
        //     continue;
        // }

		// ===============================================================================================

		penumbra = facteur_lumiere(scene, hit.position, light);
		//penumbra = 1;

        // Diffuse component: k_dλ * Sλ * (N ⋅ L_i)
		double nDotL = std::max(0.0, dot(hit.normal, lightDir));
		diffuse += material.k_diffuse * base_color * nDotL * penumbra * light.emission;

        // Specular component: k_sλ * [ m * Sλ + (1 - m) ] * (R_i ⋅ E)^shininess
		double3 R_i = normalize(2 * nDotL * hit.normal - lightDir);
		double rDotE = std::max(0.0, dot(R_i, Eye));
		double m = material.metallic;  // Usar el factor de metallicidad
		double shininess = material.shininess;

		rDotE = pow(rDotE, shininess);

		// Calcular el componente especular usando el color del material y metallicidad
		specular += material.k_specular * ((m * base_color) + (1 - m)) * rDotE * penumbra * light.emission;

    }

    // Sum the components to get the final color
    double3 outColor = ambient + diffuse + specular;
    return outColor;
}
