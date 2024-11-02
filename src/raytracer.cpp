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

	double3 vVec = normalize(scene.camera.up); 								//Direction Y de la caméra
	double3 wVec = normalize(scene.camera.center - scene.camera.position); 	//Direction Z de la caméra
	double3 uVec = normalize(cross(vVec,wVec)); 							//Direction X de la caméra
	double3 rayOrigin = scene.camera.position;	
	//toString(wVec);
							
	
	// double y_shift = 2.0 / scene.resolution[1];
	// double x_shift = 2.0 / scene.resolution[0];
	double height = 2 * scene.camera.z_near * tan(deg2rad(scene.camera.fovy/2));
	double width = height * scene.camera.aspect;
	// std::cout << "altura " << height << "\n"; 
	// std::cout << "ancho " << width << "\n"; 
	//std::cout << "aspect " << scene.camera.aspect << "\n"; 


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
				// Initialize la couleur du rayon
				double3 ray_color{0,0,0};

				// @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				// Lancez le rayon de manière uniformément aléatoire à l'intérieur du pixel dans la zone délimité par jitter_radius. 
				//Faites la moyenne des différentes couleurs obtenues suite à la récursion.
				//printf("random : %f", rand_double());
				double xDirection = (x+(rand_double()-0.5)-scene.resolution[0]/2) * (width/scene.resolution[0]);
				double yDirection = (y+(rand_double()-0.5)-scene.resolution[1]/2) * (height/scene.resolution[1]);
				//Intersection hit;
				double3 rayXDirection{xDirection, 0, 0};
				double3 rayYDirection{0, yDirection, 0};
				double3 rayDirection = normalize(wVec + rayXDirection + rayYDirection);
				Intersection hit;
				//std::cout << "Otigin " << rayOrigin.x << "," << rayOrigin.y << "," << rayOrigin.z << "\n"; 
				ray = Ray(scene.camera.position, rayDirection);
				double z_depth = scene.camera.z_far;

				if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
					Material& material = ResourceManager::Instance()->materials[hit.key_material];
					// avg_ray_color += material.color_albedo;
					avg_ray_color += shade(scene, hit);
					avg_z_depth += hit.depth;
				}
			}

			avg_z_depth = avg_z_depth / scene.samples_per_pixel;
			avg_ray_color = avg_ray_color / scene.samples_per_pixel;

			// Test de profondeur
			if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
				avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
				z_buffer[x + y*scene.resolution[0]] = avg_z_depth;

				// Met à jour la couleur de l'image (et sa profondeur)
				output->set_color_pixel(x, y, avg_ray_color);
				output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / 
										(scene.camera.z_far-scene.camera.z_near));
			}
        }
    }

    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
// 		- Détermine si le rayon intersecte la géométrie.
//      	- Calculer la contribution associée à la réflexion.
//			- Calculer la contribution associée à la réfraction.
//			- Mettre à jour la couleur avec le shading + 
//			  Ajouter réflexion selon material.reflection +
//			  Ajouter réfraction selon material.refraction 
//            pour la couleur de sortie.
//          - Mettre à jour la nouvelle profondeure.
void Raytracer::trace(const Scene& scene,
					  Ray ray, int ray_depth,
					  double3* out_color, double* out_z_depth)
{
	Intersection hit;
	// Fait appel à l'un des containers spécifiées.
	if(scene.container->intersect(ray,EPSILON,*out_z_depth,&hit)) {		
		Material& material = ResourceManager::Instance()->materials[hit.key_material];
		//printf("INtersec\n");
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réflection d'un rayon de manière récursive.
		
		// @@@@@@ VOTRE CODE ICI
		// Déterminer la couleur associée à la réfraction d'un rayon de manière récursive.
		// 
		// Assumez que l'extérieur/l'air a un indice de réfraction de 1.
		//
		// Toutes les géométries sont des surfaces et non pas de volumes.

		// *out_color = 
		// *out_z_depth =
	} else {
		//printf("NP-INtersec\n");
	}

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
	// lorsque vous serez rendu à la partie texture.
	Material& material = ResourceManager::Instance()->materials[hit.key_material]; 
	
	double3 ambient(0, 0, 0);
    double3 diffuse(0, 0, 0);
    double3 specular(0, 0, 0);

    // Vector pointing towards the eye
    double3 Eye = normalize(scene.camera.position - hit.position);

    // Ambient light calculation: L_aλ * k_aλ * Sλ
    ambient = scene.ambient_light * material.k_ambient * material.color_albedo;

    // Iterate through each light in the scene
    for (const auto& light : scene.lights) {

		// ================= verifier cette partie du code, cela donne de problems avec le cilindre ====================

        // Light direction and distance to the light source
        double3 lightDir = normalize(light.position - hit.position);
        double lightDistance = length(light.position - hit.position);

        // Shadow ray: emitted from the hit position towards the light source
        // Offset the starting position to avoid "surface acne" issues
        Ray shadowRay(hit.position + EPSILON * hit.normal, lightDir);

        // Check if any object obstructs the shadow ray up to the light source
        bool in_shadow = false;
        // for (const auto& obj : scene.objects) {
        //     // If the shadow ray intersects an object within the light distance, the point is in shadow
        //     if (obj->intersect(shadowRay, 1e-6, lightDistance)) {
        //         in_shadow = true;
        //         break;
        //     }
        // }
		if (scene.container->intersect(shadowRay, EPSILON, lightDistance, &hit)) in_shadow = true;

        // If the point is in shadow, skip diffuse and specular contributions for this light
        if (in_shadow) {
            continue;
        }

		// ===============================================================================================

        // Diffuse component: 2 * k_dλ * Sλ * (N ⋅ L_i)
        double nDotL = std::max(0.0, dot(hit.normal, lightDir));
        diffuse += 2 * material.k_diffuse * material.color_albedo * nDotL;

        // Specular component: k_sλ * [ m * Sλ + (1 - m) ] * (R_i ⋅ E)
        double3 R_i = normalize(2 * nDotL * hit.normal - lightDir);
        double rDotE = std::max(0.0, dot(R_i, Eye));
        double m = material.k_reflection;
        specular += material.k_specular * ((m * material.color_albedo) + (1 - m)) * rDotE;
    }

    // Sum the components to get the final color
    double3 outColor = ambient + diffuse + specular;
    return outColor;
}
