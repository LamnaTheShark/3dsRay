

	#include <stdio.h>

	#include <cstdlib>
	#include <cstdio>
	#include <cmath>
	#include <fstream>
	#include <vector>
	#include <iostream>
	#include <cassert>


	#define M_PI 3.141592653589793
	#define INFINITY 1e12

	unsigned int renderBitmap[256 * 256] = { 0 };

	template<typename T>
	class Vec3
	{
	public:
		T x, y, z;
		Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
		Vec3(T xx) : x(xx), y(xx), z(xx) {}
		Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
		Vec3& normalize()
		{
			T nor2 = length2();
			if (nor2 > 0) {
				T invNor = 1 / sqrt(nor2);
				x *= invNor, y *= invNor, z *= invNor;
			}
			return *this;
		}
		Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
		Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
		T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
		Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
		Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
		Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
		Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
		Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
		T length2() const { return x * x + y * y + z * z; }
		T length() const { return sqrt(length2()); }
		friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
		{
			os << "[" << v.x << " " << v.y << " " << v.z << "]";
			return os;
		}
	};

	typedef Vec3<float> Vec3f;

	class Sphere
	{
	public:
		Vec3f center;                           /// position of the sphere
		float radius, radius2;                  /// sphere radius and radius^2
		Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
		float transparency, reflection;         /// surface transparency and reflectivity
		Sphere(
			const Vec3f &c,
			const float &r,
			const Vec3f &sc,
			const float &refl = 0,
			const float &transp = 0,
			const Vec3f &ec = 0) :
			center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
			transparency(transp), reflection(refl)
		{ /* empty */
		}
		//[comment]
		// Compute a ray-sphere intersection using the geometric solution
		//[/comment]
		bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
		{
			Vec3f l = center - rayorig;
			float tca = l.dot(raydir);
			if (tca < 0) return false;
			float d2 = l.dot(l) - tca * tca;
			if (d2 > radius2) return false;
			float thc = sqrt(radius2 - d2);
			t0 = tca - thc;
			t1 = tca + thc;

			return true;
		}
	};



	//[comment]
	// This variable controls the maximum recursion depth
	//[/comment]
	#define MAX_RAY_DEPTH 5

	float mix(const float &a, const float &b, const float &mix)
	{
		return b * mix + a * (1 - mix);
	}

	//[comment]
	// This is the main trace function. It takes a ray as argument (defined by its origin
	// and direction). We test if this ray intersects any of the geometry in the scene.
	// If the ray intersects an object, we compute the intersection point, the normal
	// at the intersection point, and shade this point using this information.
	// Shading depends on the surface property (is it transparent, reflective, diffuse).
	// The function returns a color for the ray. If the ray intersects an object that
	// is the color of the object at the intersection point, otherwise it returns
	// the background color.
	//[/comment]
	Vec3f trace(
		const Vec3f &rayorig,
		const Vec3f &raydir,
		const std::vector<Sphere> &spheres,
		const int &depth)
	{

		int pixelCounterS = 0;
		//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
		float tnear = INFINITY;
		const Sphere* sphere = NULL;
		// find intersection of this ray with the sphere in the scene
		for (unsigned i = 0; i < spheres.size(); ++i) {
			float t0 = INFINITY, t1 = INFINITY;
			if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
				if (t0 < 0) t0 = t1;
				if (t0 < tnear) {
					tnear = t0;
					sphere = &spheres[i];
					//iprintf(".");
				}
				//iprintf(";");
			}
			//iprintf(":");
		}
		// if there's no intersection return black or background color
		if (!sphere) return Vec3f(0);
		Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
		Vec3f phit = rayorig + raydir * tnear; // point of intersection
		Vec3f nhit = phit - sphere->center; // normal at the intersection point
		nhit.normalize(); // normalize normal direction
		// If the normal and the view direction are not opposite to each other
		// reverse the normal direction. That also means we are inside the sphere so set
		// the inside bool to true. Finally reverse the sign of IdotN which we want
		// positive.
		float bias = 1e-4; // add some bias to the point from which we will be tracing
		bool inside = false;
		if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
		if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
			float facingratio = -raydir.dot(nhit);
			// change the mix value to tweak the effect
			float fresneleffect;
			if(sphere->transparency > 0 && sphere->reflection > 0){
				fresneleffect = mix(pow(1 - facingratio, 3), 1.5, 0.2);
			}
			else{
				fresneleffect = mix(pow(1 - facingratio, 3), 1.5, 0.5);
			}
			
			// compute reflection direction (not need to normalize because all vectors
			// are already normalized)
			Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
			refldir.normalize();
			Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
			Vec3f refraction = 0;
			//iprintf("/");
			// if the sphere is also transparent compute refraction ray (transmission)
			if (sphere->transparency) {
				float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
				float cosi = -nhit.dot(raydir);
				float k = 1 - eta * eta * (1 - cosi * cosi);
				Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
				refrdir.normalize();
				refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
				//iprintf("?");
			}
			// the result is a mix of reflection and refraction (if the sphere is transparent)
			surfaceColor = (
				reflection * fresneleffect +
				refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
		}
		else {
			// it's a diffuse object, maybe a light
			
			for (unsigned i = 0; i < spheres.size(); ++i) {
				
				if (spheres[i].emissionColor.x > 0) {
					// this is a light
					Vec3f transmission = 1;
					Vec3f lightDirection = spheres[i].center - phit;
					lightDirection.normalize();
					for (unsigned j = 0; j < spheres.size(); ++j) {
						if (i != j) {
							float t0, t1;
							if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
								
								transmission = 1;//Was 0 before
								break;
								//iprintf("!");
							}
							//iprintf("§");
						}
						//iprintf("*");
					}
					//iprintf("µ");
					//surfaceColor = lightDirection;
					
						//DISABLE EMISSION OF DIFFUSE SURFACES
						surfaceColor = sphere->surfaceColor * transmission *
						 std::max(float(0), nhit.dot(lightDirection)) 
						* spheres[i].emissionColor;

					//iprintf("ù");
				}
				else{ //It's a diffuse object, must calculate caustics

					for (unsigned k = 0; k < spheres.size(); ++k) { //Parsing all spheres
				
						if (spheres[k].emissionColor.x > 0) {
							// this is a light, getting path from light to diffuse sphere
							Vec3f transmission = 1;
							Vec3f lightDirection = spheres[k].center - phit; //Vecteur du hitPoint vers la lumière
							lightDirection.normalize();
							for (unsigned j = 0; j < spheres.size(); ++j) {
								if (k != j) {
									float t0, t1;
									if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
									
										//transmission = 1;//We found a sphere that intersect ray between hitPoint & light
										if(sphere[j].transparency ){
											//That sphere is transparent, we have a caustic !

											float ior = 1.1, eta = (inside) ? ior : 1 / ior;
											float cosi = -phit.dot(lightDirection);
											float ki = 1 - eta * eta * (1 - cosi * cosi);
											Vec3f refrdir = raydir * ior + phit * (ior *  cosi - sqrt(ki));
											refrdir.normalize();
											//transmission = trace(phit - phit * bias, refrdir, spheres, depth + 1);
											
											transmission =2;
											//surfaceColor +=  0.02 ;
											
											//printf("*");
										}
										//break;
									
									//iprintf("!");
									}
								//iprintf("§");
								}
							//iprintf("*");
							}
						//iprintf("µ");
						}
				}
			}
				
				
				
				//iprintf("%");
			}
			//printf(".");
			pixelCounterS++;
			//iprintf("pCount : %d", pixelCounterS);
		}
		//iprintf("£");

		return surfaceColor + sphere->emissionColor;
	}


	//[comment]
	// Main rendering function. We compute a camera ray for each pixel of the image
	// trace it and return a color. If the ray hits a sphere, we return the color of the
	// sphere at the intersection point, else we return the background color.
	//[/comment]
	void render(const std::vector<Sphere> &spheres)
	{
		printf("\tInit");
		unsigned width = 256, height = 256;
		Vec3f *image = new Vec3f[width * height], *pixel = image;
		float invWidth = 1 / float(width), invHeight = 1 / float(height);
		float fov = 30, aspectratio = width / float(height);
		float angle = tan(M_PI * 0.5 * fov / 180.);
		printf("\tCam Setup finished");
		// Trace rays
		for (unsigned y = 0; y < height; ++y) {
			
			for (unsigned x = 0; x < width; ++x, ++pixel) {
				float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
				float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
				Vec3f raydir(xx, yy, -1);
				raydir.normalize();
				*pixel = trace(Vec3f(0), raydir, spheres, 0);
				
			}
			
		}
		printf("\tRender written to vector");
		// Save result to a PPM image (keep these flags if you compile under Windows)
		std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
		ofs << "P6\n" << width << " " << height << "\n255\n";
		
		int totalPixelSize = width * height;
		int totalIt = 0;
		int ib = 0;
		for (unsigned i = 0; i < width * height; ++i) {
			//iprintf("[%d]",i);
			//renderBitmap[i] = (u16)(((image[i].x + image[i].y + image[i].z) * 20)/3);
			//renderBitmap[i + 1] = (u8)image[i].y;
			//renderBitmap[i + 2] = (u8)image[i].z;
			/*
			renderBitmap[ib] = (u8)image[i].x;
			renderBitmap[ib + 1] = (u8)image[i].y;
			renderBitmap[ib + 2] = (u8)image[i].z;
			*/
			//ib += 3;
			//totalIt++;
			//iprintf("[%d]=[%d]",i,totalPixelSize);
			ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
				(unsigned char)(std::min(float(1), image[i].y) * 255) <<
				(unsigned char)(std::min(float(1), image[i].z) * 255);
		}
		printf("Written[%d]=[%d]", totalIt, totalPixelSize);
		printf("\tRender written to table");
		ofs.close();
		
		delete[] image;
	}




	int main(void)
	{
	    // set the mode for 2 text layers and two extended background layers
		/*videoSetMode(MODE_5_2D);
		vramSetPrimaryBanks(VRAM_A_MAIN_BG_0x06000000, VRAM_B_MAIN_BG_0x06020000,
			VRAM_C_SUB_BG, VRAM_D_LCD);

		consoleDemoInit();

		iprintf("\n\n\tSmall Raytracing DS engine\n");
		iprintf("\tStarting Render...\n");
	*/
		srand48(13);
		std::vector<Sphere> spheres;
		// position, radius, surface color, reflectivity, transparency, emission color

		spheres.push_back(Sphere(Vec3f(0.0, 6.0, -30.0), 1, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(1.2)));// light


		spheres.push_back(Sphere(Vec3f(0.0, -1e4 - 3, 0), 1e4, Vec3f(0.7, 0.7, 0.7), 0, 0.0)); //Floor
		spheres.push_back(Sphere(Vec3f(0.0, 1e4 + 8, 0), 1e4, Vec3f(0.7, 0.7, 0.7), 0, 0.0)); //Roof
		spheres.push_back(Sphere(Vec3f(0.0, 0.0, -1e4 - 45), 1e4, Vec3f(0.5, 0.7, 0), 0, 0.0)); //Back
		spheres.push_back(Sphere(Vec3f(-1e4 - 8, 0.0, 0.0), 1e4, Vec3f(0, 0.5, 0.7), 0, 0.0)); //Left
		spheres.push_back(Sphere(Vec3f(1e4 + 8, 0.0, 0.0), 1e4, Vec3f(0.7, 0, 0.0), 0, 0.0)); //Right


		spheres.push_back(Sphere(Vec3f(-2.0, 0, -30), 3, Vec3f(1.00, 1.00, 1.00), 1, 1.2)); //Glass
		spheres.push_back(Sphere(Vec3f(2.0, 0, -35), 3, Vec3f(1.00, 1.00, 1.00), 1, 0.0)); //Mirror
		spheres.push_back(Sphere(Vec3f(3.0, -1, -20), 1, Vec3f(1.00, 1.00, 1.00), 1, 0.0)); //Mirror
		spheres.push_back(Sphere(Vec3f(-6.0, -1, -38), 1, Vec3f(1.00, 1.00, 1.00), 1, 0.0)); //Mirror
		
		
		render(spheres);

		printf("\tRender finished");
	/*
		int bg3 = bgInit(3, BgType_Bmp16, BgSize_B16_256x256, 0,0);

		//dmaCopy((u16*)renderBitmap, bgGetGfxPtr(bg3), 256 * 256);
		//dmaCopy(drunkenlogoPal, BG_PALETTE, 256*2);

		u16 colorMask = 0xFF;

		u16* backBuffer = (u16*)bgGetGfxPtr(bg3) + 256 * 256;

		int xb = 0;
		int yb = 0;
		for (int iy = 0; iy < 256; iy++){
			for (int ix = 0; ix < 256; ix++){
				backBuffer[iy * 256 + ix] = ((renderBitmap[iy * 256 + ix]) & colorMask) | BIT(15);
				xb += 3;
			}
			xb = 0;
			yb += 1;
		}

		backBuffer = (u16*)bgGetGfxPtr(bg3);
		bgSetMapBase(bg3, 8);

		while(1) {
			swiWaitForVBlank();
			scanKeys();
			if (keysDown()&KEY_START) break;
		}
	*/
		return 0;
	}







