// rt: un lanzador de rayos minimalista
// g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
//#include <random> //para generar numeros aleatorios para el metodo de monte carlo

class Vector 
{
public:        
	double x, y, z; // coordenadas x,y,z 
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }
  
	// operator % para producto cruz
	Vector operator%(Vector&b){return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
	
	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
	
	// normalizar vector 
	Vector& normalize(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;
typedef Vector Color;

/*
Un rayo está formado por un origen o y una dirección d. El constructor recibe ambos valores. Matemáticamente,
cualquier punto del rayo se expresa como r(t) = O + tD, donde t representa la distancia recorrida
*/
class Ray 
{ 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

class Sphere 
{
public:
	double r;	// radio de la esfera
	Point p;	// posicion
	Color c;	// color  

	Sphere(double r_, Point p_, Color c_): r(r_), p(p_), c(c_) {}
  
	// PROYECTO 1
	// determina si el rayo intersecta a esta esfera
	double intersect(const Ray &ray) const {
		// regresar distancia si hay intersección
		// regresar 0.0 si no hay interseccion
		Vector op = p-ray.o; //Calculamos el vector desde el origen del rayo hasta el centro de la esfera: op = P-O
		double b = op.dot(ray.d); //se calcula el producto punto entre op y la dirección del rayo: b = op . D
		double det = b*b - op.dot(op) + r*r; //es el discriminante de la ecuacion det = d*d - op*op + r*r
		if(det < 0.0) // si el discriminante es negativo, no hay intersección, por lo que se retorna 0.0
			return 0.0; //El rayo no toca la esfera.
		
		det = sqrt(det); //si existe la interseccion, entonces se calcula las posibles intersecciones t1 y t2:
		/*tambien puede ocurrir:
		t1 < 0
		t2 > 0
		por eso comprobamos el t2 por que Esto sucede cuando el origen del rayo está dentro de la esfera.
		y se comprueba el t1 y t2
		*/
		double t1 = b - det;//t1 es la interseccion mas cercana
    	double t2 = b + det;// Si t1 no sirve, probar la t2
		if(t1 > 1e-4) //usamos 1e-4 para evitar intersecciones muy cercanas a 0
			return t1;
		if(t2 > 1e-4)
			return t2;
	

		return 0.0; //no hay interseccion valida
	}
};

/*
el arreglo spheres[]:contiene las superficies de la escena. Las cinco primeras esferas tienen un radio enorme y
funcionan como paredes: izquierda, derecha, fondo, suelo y techo. Esto crea una habitación. Después se agregan tres
esferas pequeñas ubicadas dentro de ella. Sus colores se expresan con tres valores RGB en el rango aproximado de 0 a 1.
*/
Sphere spheres[] = {
	//Escena: radio, posicion, color 
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),   Color(0.75, 0.25, 0.25)), // pared izq
	Sphere(1e5,  Point(1e5 + 49, 0, 0),    Color(0.25, 0.25, 0.75)), // pared der
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6), Color(.50, .50, .50)), // pared detras
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0), Color(0.50, 0.50, 0.50)), // suelo
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(0.50, 0.50, 0.50)), // techo
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(.999, .999, .999)), // esfera abajo-izq
	Sphere(16.5, Point(23, -24.3, -3.6),   Color(.999, .999, .999) ), // esfera abajo-der
	Sphere(10.5, Point(0, 24.3, 0),        Color(.999, .999, .999)), // esfera arriba
	//Sphere(10.5, Point(0,-10,-2), 		   Color(1, 1, 0)) 
};

//Point lightPosition(0, 30, 150);
//Color lightColor(1, 1, 1);

//==============================================================================================================================================
/*
la funcion clamp() limita cada componente de color al intervalo [0,1]. Si un valor es menor que cero se convierte en 0 y si es
mayor que uno se convierte en 1. 
la funcion toDisplayValue() transforma posteriormente el valor a un entero entre 0 y 255 y
aplica una corrección gamma aproximada de 2.2, adecuada para visualizar los colores en la imagen.
*/
inline double clamp(const double x) { 
	if(x < 0.0)
		return 0.0;
	else if(x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
	return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); 
}
//==============================================================================================================================================

// PROYECTO 1
// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
//int ignoreId = -1
inline bool intersect(const Ray &r, double &t, int &id,int ignoreId = -1) {
	double distance;
    double closest = 1e20;

    bool hit = false;

    int numSpheres = sizeof(spheres) / sizeof(Sphere); // recorremos todas las esferas en la escena

    for (int i = 0; i < numSpheres; i++) { //considerando 3 esferas.

		if(i==ignoreId) // ignorar la esfera con el id especificado
			continue;

        distance = spheres[i].intersect(r);//probamos cada esfera

        if (distance > 0.0 && distance < closest) {//eligo la esfera mas cercana a la camara

            closest = distance;
            id = i;
            hit = true;
        }
    }

    if (hit) {
        t = closest;
        return true;
    }

    return false;
}

// Calcula el valor de color para el rayo dado
/*
Primero se busca una intersección mediante intersect(). Si no se encuentra ningún objeto, se devuelve el color
negro. Si existe una intersección, se obtiene la esfera mediante su índice id.
1.- El punto de impacto se calcula con x = r.o + r.d * t. 
2.- Después se obtiene la normal de la esfera en ese punto mediante (x - obj.p).normalize(). 
3.- En la versión actual, la normal se calcula pero todavía no participa en el color final. Finalmente se devuelve obj.c, por lo que el objeto aparece con su color base sin iluminación.
*/
Color shade(const Ray &r) {

    double t;
    int id = 0;

	// determinar que esfera (id) y a que distancia (t) el rayo intersecta
	// ¿El rayo golpeó alguna esfera?
    if (!intersect(r, t, id))
        return Color();

    // Obtener la esfera golpeada
    const Sphere &obj = spheres[id];

	// PROYECTO 1
	// determinar coordenadas del punto de interseccion
    // Punto de intersección r(t) = O + t*D
    Point x = r.o + r.d * t;

    // determinar la dirección normal en el punto de interseccion
    Vector n = (x - obj.p).normalize();

    // Color final
    Color colorValue = obj.c;

    return colorValue;
}

int main(int argc, char *argv[]) {
	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );
	

	// parametros de la camara
	//son vectores auxiliares que definen los ejes horizontal y vertical del plano de imagen. El producto cruz ayuda a obtener una dirección perpendicular.
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;

	// auxiliar para valor de pixel y matriz para almacenar la imagen
	// en pocas palabras, se reserva memoria dinamica para almacenar los colores por cada pixel de la imagen.
	Color *pixelColors = new Color[w * h];

	// PROYECTO 1
	// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
	#pragma omp parallel for
	for(int y = 0; y < h; y++) 
	{ 
		// recorre todos los pixeles de la imagen
		fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
		for(int x = 0; x < w; x++ ) {
			//convierte las coordenadas 2D en una posición del arreglo lineal y además invierte verticalmente la imagen para el formato de salida.
			int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
			Color pixelValue = Color(); // pixelValue en negro por ahora
			// para el pixel actual, computar la dirección que un rayo debe tener
			// se construye cameraRayDir combinando la posición horizontal, vertical y dirección de la cámara. Luego se normaliza.
			Vector cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
			
			// computar el color del pixel para el punto que intersectó el rayo desde la camara
			pixelValue = shade( Ray(camera.o, cameraRayDir.normalize()) );

			// limitar los tres valores de color del pixel a [0,1]
			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
	}

	fprintf(stderr,"\n");

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("image.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	return 0;
}


