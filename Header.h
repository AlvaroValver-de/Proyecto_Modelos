// Incluye una vez para evitar inclusiones múltiples del mismo encabezado.
#pragma once

// Previene redefiniciones en caso de que otro archivo incluya este encabezado.
#ifndef HeaderH
#define HeaderH

// Librerías estándar necesarias para el código.
#include <stdio.h>
#include <math.h>
#include <glm/glm.hpp>


// Establece el espacio de nombres a utilizar, en este caso 'std'.
using namespace std;

//---------- Definiciones Globales ----------

// Estos son identificadores para facilitar la identificación de ciertos valores.
#define None                    -1
#define xy                      0
#define xz                      1
#define yz                      2

// Constantes físicas y matemáticas.
#define V_SON                   340.0   // Velocidad del sonido en m/s.
#define PI                      3.1415926535897932384626433832795  // Aproximación del valor de PI.
#define MaxNPoints              200     // Número máximo de puntos permitidos.

//---------------------------------------------------------------------------
class vector1 {
public:
    double x, y, z;  // Coordenadas del vector.

    // Suma de vectores.
    vector1 operator+(vector1 v2) {
        vector1 v1;
        v1.x = x + v2.x;
        v1.y = y + v2.y;
        v1.z = z + v2.z;
        return v1;
    }

    // Resta de vectores.
    vector1 operator-(vector1 v2) {
        vector1 v1;
        v1.x = x - v2.x;
        v1.y = y - v2.y;
        v1.z = z - v2.z;
        return v1;
    }

    // Multiplicación por escalar.
    vector1 operator*(double f) {
        vector1 v;
        v.x = x * f;
        v.y = y * f;
        v.z = z * f;
        return v;
    }

    // División por escalar.
    vector1 operator/(double f) {
        vector1 v;
        v.x = x / f;
        v.y = y / f;
        v.z = z / f;
        return v;
    }

    // Producto punto entre vectores.
    double operator*(vector1 v) {
        return x * v.x + y * v.y + z * v.z;
    }
    vector1 Normal() { // Calcula el vector unitario
        double m = Module();
        vector1 v1;
        vector1 v2;
        v1.x = x;
        v1.y = y;
        v1.z = z;
        if (m == 0)
            v2 = 0;
        else
            v2 = v1 / m;
        return v2;
    }
    // Producto cruz.
    vector1 operator/(vector1 v2) {
        vector1 v1;
        v1.x = y * v2.z - z * v2.y;
        v1.y = -x * v2.z + z * v2.x;
        v1.z = x * v2.y - y * v2.x;
        return v1;
    }

    // Asignar mismo valor a x, y y z.
    void operator=(double f) {
        x = y = z = f;
    }

    // Valor absoluto de las coordenadas.
    vector1 Abs() {
        vector1 v;
        v.x = fabs(x);
        v.y = fabs(y);
        v.z = fabs(z);
        return v;
    }

    // Módulo o magnitud del vector.
    double Module() {
        return sqrt(x * x + y * y + z * z);
    }

    // Retorna el vector unitario.
    vector1 Unitary() {
        double m;
        vector1 u;
        u.x = x;
        u.y = y;
        u.z = z;
        m = u.Module();
        if (m == 0) {
            u.x = 0.0;
            u.y = 0.0;
            u.z = 0.0;
        }
        else {
            u = u / m;
        }
        return u;
    }
};
//---------------------------------------------------------------------------
class point {
public:
    double x, y, z;  // Coordenadas del punto.

    // Constructor predeterminado que inicializa el punto en el origen.
    point() : x(0), y(0), z(0) {}

    // Traslación por un vector.
    point operator+(vector1 v) {
        point p;
        p.x = x + v.x;
        p.y = y + v.y;
        p.z = z + v.z;
        return p;
    }

    // Suma de dos puntos.
    point operator+(point p2) {
        point p1;
        p1.x = x + p2.x;
        p1.y = y + p2.y;
        p1.z = z + p2.z;
        return p1;
    }

    // Resta de dos puntos que da lugar a un vector.
    vector1 operator-(point p) {
        vector1 v;
        v.x = x - p.x;
        v.y = y - p.y;
        v.z = z - p.z;
        return v;
    }

    // Multiplicación por escalar.
    point operator*(double f) {
        point p;
        p.x = x * f;
        p.y = y * f;
        p.z = z * f;
        return p;
    }

    // División por escalar.
    point operator/(double f) {
        point p;
        p.x = x / f;
        p.y = y / f;
        p.z = z / f;
        return p;
    }

    // Asignar mismo valor a x, y y z.
    void operator=(double f) {
        x = y = z = f;
    }

    // Comprueba si dos puntos son iguales.
    bool operator==(point p) {
        return (x == p.x && y == p.y && z == p.z);
    }

    // Comprueba si dos puntos son diferentes.
    bool operator!=(point p) {
        return !(*this == p);  // Reutiliza el operador de igualdad.
    }

    // Resetea el punto a su valor predeterminado.
    void Clear() {
        x = 0;
        y = 0;
        z = 0;
    }

    // Valor absoluto de las coordenadas.
    point Abs() {
        point v;
        v.x = fabs(x);
        v.y = fabs(y);
        v.z = fabs(z);
        return v;
    }

    // Calcula la distancia entre dos puntos.
    double distancia(point p2) {
        return sqrt((p2.x - x) * (p2.x - x) + (p2.y - y) * (p2.y - y) + (p2.z - z) * (p2.z - z));
    }
};
//---------------------------------------------------------------------------
class triangle {
public:
    point p0, p1, p2, bc;  // Vértices del triángulo y baricentro
    int Projection;        // Proyección
    double a0;             // Constante a0 para cálculos futuros
    int ID;                // Identificador único

    // Constructor por defecto
    triangle() {
        p0 = 0;
        p1 = 0;
        p2 = 0;
        bc = 0;
        Projection = 0;
        a0 = 0;
        ID = 0;
    };

    // Sobrecarga del operador =
    void operator=(const triangle& t) {
        p0 = t.p0;
        p1 = t.p1;
        p2 = t.p2;
        bc = t.bc;
        Projection = t.Projection;
        a0 = t.a0;
        ID = t.ID;
    }

    // Función para resetear el triángulo
    void Clear() {
        p0 = 0;
        p1 = 0;
        p2 = 0;
        bc = 0;
        Projection = 0;
        a0 = 0;
        ID = 0;
    }

    // Función para calcular el baricentro del triángulo
    void Centroid() {
        bc = (p0 + p1 + p2) / 3;
    }

    // Función para calcular el área del triángulo usando producto cruzado
    double TriangleArea() {
        vector1 v = (p1 - p0) / (p2 - p0);  // Suponiendo que / es producto cruzado
        return 0.5 * v.Module();
    }

    // Función para calcular la proyección del triángulo
    void CalculateProjection() {
        vector1 n = (p1 - p0) / (p2 - p0);  // Suponiendo que / es producto cruzado
        n.x = n.x * n.x;
        n.y = n.y * n.y;
        n.z = n.z * n.z;

        // Elegir proyección basada en la componente más grande del vector normal
        if ((n.x >= n.y) && (n.x >= n.z)) {
            Projection = yz;
            a0 = 1 / (-p1.y * p0.z + p2.y * p0.z + p0.y * p1.z - p2.y * p1.z - p0.y * p2.z + p1.y * p2.z + 0.000001);
        }
        else if ((n.y >= n.x) && (n.y >= n.z)) {
            Projection = xz;
            a0 = 1 / (-p1.x * p0.z + p2.x * p0.z + p0.x * p1.z - p2.x * p1.z - p0.x * p2.z + p1.x * p2.z + 0.000001);
        }
        else {
            Projection = xy;
            a0 = 1 / (-p1.x * p0.y + p2.x * p0.y + p0.x * p1.y - p2.x * p1.y - p0.x * p2.y + p1.x * p2.y + 0.000001);
        }
    }

    // Función para calcular el ángulo sólido con respecto a un punto
    double AnguloSolido(point bc) {
        double d = 0.2;
        triangle t;
        vector1 v1, v2, v3;
        v1 = p0 - bc;
        v2 = p1 - bc;
        v3 = p2 - bc;
        v1 = v1 / v1.Module();
        v2 = v2 / v2.Module();
        v3 = v3 / v3.Module();
        t.p0 = bc + (v1 * d);
        t.p1 = bc + (v2 * d);
        t.p2 = bc + (v3 * d);
        return t.TriangleArea();
    }
};
//---------------------------------------------------------------------------
class plane {
public:
    int         NP;      // Número de puntos en el plano
    point* p;            // Puntero a los puntos del plano
    int         NT;      // Número de triángulos en el plano
    triangle* t;         // Puntero a los triángulos del plano
    vector1      n;      // Vector normal al plano

    // Constructor por defecto
    plane() {
        NP = 0;
        p = NULL;
        NT = 0;
        t = NULL;
        n = vector1(); // Inicializar el vector normal
    }

    // Función para agregar nuevos puntos al plano
    void NewPoints(int N) {
        point* tp;
        tp = new point[NP + N];

        if (NP > 0) {
            delete[] p;
            p = NULL;
        }
        p = tp;
        NP += N;
    }

    // Función para eliminar un punto específico del plano
    void DeletePoint(int IP) {
        int j = 0;
        if (IP >= 0 && IP < NP) {
            point* tp;
            tp = new point[NP - 1];
            for (int P = 0; P < NP; P++) {
                if (P != IP) {
                    tp[j] = p[P];
                    j++;
                }
            }
            delete[] p;
            p = tp;
            NP -= 1;
        }
    }

    // Función para agregar nuevos triángulos al plano
    void NewTriangle(int N) {
        triangle* tt;
        tt = new triangle[NT + N];
        for (int T = 0; T < NT; T++) {
            tt[T] = t[T];
        }
        for (int T = NT; T < NT + N; T++) {
            tt[T].Clear();
        }
        if (NP > 0) {
            delete[] t;
            t = NULL;
        }
        t = tt;
        NT += N;
    }

    // Función para eliminar un triángulo específico del plano
    void DeleteTriangle(int IT) {
        int j = 0;
        if (IT >= 0 && IT < NT) {
            triangle* tt;
            tt = new triangle[NT - 1];
            for (int T = 0; T < NT; T++) {
                if (T != IT) {
                    tt[j] = t[T];
                    j++;
                }
            }
            delete[] t;
            t = tt;
            NT -= 1;
        }
    }

    // Función para generar dos triángulos a partir de los vértices de un cuadrado en el plano
    void PointGenTriangle() {
        NewTriangle(NP - 2);
        int i = 1;
        for (int T = 0; T < NT; T++) {
            i--;
            t[T].p0 = p[i];
            i++;
            if (i == NP) i = 0;
            t[T].p1 = p[i];
            i++;
            if (i == NP) i = 0;
            t[T].p2 = p[i];
            i++;
        }
    }

    // Función para calcular el módulo (longitud) de un vector
    double Module(vector1 v) {
        return sqrt(v * v); // Suponiendo que * es un operador sobrecargado para el producto punto
    }

    // Función para obtener un vector unitario a partir de un vector dado
    vector1 VectorUnitario(vector1 v1) {
        double m = Module(v1);
        if (m == 0) {
            return vector1(); // Devolver un vector nulo si el módulo es 0
        }
        return v1 / m; // Suponiendo que / es una sobrecarga de operador para dividir un vector por su magnitud
    }

    // Función para calcular el vector normal al plano
    vector1 NormalPlano() {
        return VectorUnitario((p[1] - p[0]) / (p[2] - p[0])); // Suponiendo que / es un operador sobrecargado para el producto cruzado
    }

    // Función para resetear la información del plano
    void Clear() {
        NP = 0;
        delete[] p;
        p = NULL;
        NT = 0;
        delete[] t;
        t = NULL;
        n = vector1(); // Inicializar el vector normal
    }
};
//---------------------------------------------------------------------------
// Estructura para almacenar información sobre reflexiones
struct reflection {
    bool lost;        // Bandera para rayos perdidos
    point r[50];      // Posiciones de reflexión (máximo 50 reflexiones)
    double d[50];     // Distancias para las reflexiones
    int idTriangle[51];
    int Plane[51];
    int Triangle[51];
    int N;            // Número de reflexiones
};
//---------------------------------------------------------------------------
class room {
public:
    int		NP;		// Número de Planos
    plane* p;		// Puntero a Planos
    double		distMax;		// Distancia máxima entre dos puntos en la sala.
    int     NR;  // Número de receptores

    room() {
        NP = 0;
        p = NULL;
        NR = 0;
        distMax = 0.0;
    }

   

    double IntersectionDistance(vector1 n, point pnt, vector1 u, point o) {
        double m = n * u;
        if (fabs(m) < 1e-6) return -1;
        return (n * (pnt - o)) / m;
    }

    void MaxDistance() {
        distMax = 0;
        float tmpd = 0;
        for (int i1 = 0; i1 < NP; i1++) {
            for (int j1 = 0; j1 < p[i1].NP; j1++) {
                for (int i2 = 0; i2 < NP; i2++) {
                    for (int j2 = 0; j2 < p[i2].NP; j2++) {
                        tmpd = p[i1].p[j1].distancia(p[i2].p[j2]);
                        if (distMax < tmpd)
                            distMax = tmpd;
                    }
                }
            }
        }
    };

    void NewPlanes(int N) {
        plane* tp = new plane[NP + N];
        for (int P = 0; P < NP; P++) {
            tp[P] = std::move(p[P]);
        }
        for (int P = NP; P < NP + N; P++) {
            tp[P].Clear();
        }
        delete[] p;
        p = tp;
        NP += N;
    }

    bool Inner(point p, triangle t) {
        double a1, a2, x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2;

        x = p.x;
        y = p.y;
        z = p.z;

        x0 = t.p0.x;
        y0 = t.p0.y;
        z0 = t.p0.z;
        x1 = t.p1.x;
        y1 = t.p1.y;
        z1 = t.p1.z;
        x2 = t.p2.x;
        y2 = t.p2.y;
        z2 = t.p2.z;

        if (t.Projection == yz) {                              //Proje  o yz
            a1 = -t.a0 * (-y0 * z + y2 * z + y * z0 - y2 * z0 - y * z2 + y0 * z2);
            a2 = -t.a0 * (y0 * z - y1 * z - y * z0 + y1 * z0 + y * z1 - y0 * z1);
        }
        if (t.Projection == xz) {                              //Proje  o xz
            a1 = -t.a0 * (-x0 * z + x2 * z + x * z0 - x2 * z0 - x * z2 + x0 * z2);
            a2 = -t.a0 * (x0 * z - x1 * z - x * z0 + x1 * z0 + x * z1 - x0 * z1);
        }
        if (t.Projection == xy) {                              //Proje  o xy
            a1 = -t.a0 * (-x2 * y0 + x * y0 + x0 * y2 - x * y2 - x0 * y + x2 * y);
            a2 = t.a0 * (-x1 * y0 + x * y0 + x0 * y1 - x * y1 - x0 * y + x1 * y);
        }

        if ((a1 + a2 <= 1) && (a1 >= 0) && (a2 >= 0))
            return true;
        else
            return false;

    };

    // Función para calcular la reflexión de un vector incidente respecto a una normal
    vector1 Reflect(vector1 incident, vector1 normal) {
        // Calcular el producto punto entre el vector incidente y la normal
        double dotProduct = incident.x * normal.x + incident.y * normal.y + incident.z * normal.z;

        // Calcular el vector de reflexión usando la fórmula R = I - 2 * (I · N) * N
        vector1 reflection;
        reflection.x = incident.x - 2 * dotProduct * normal.x;
        reflection.y = incident.y - 2 * dotProduct * normal.y;
        reflection.z = incident.z - 2 * dotProduct * normal.z;

        return reflection;
    }


    // Función que realiza el Ray Tracing en una escena y devuelve un array de reflexiones
    reflection* RayTracing(point origin, vector1* rays, int numRays) {
        reflection* reflections = new reflection[numRays];

        // Iterar a través de cada rayo en el conjunto de rayos
        for (int rayIdx = 0; rayIdx < numRays; rayIdx++) {
            vector1 ray = rays[rayIdx];

            // Inicializar la estructura de reflexión para el rayo actual
            reflections[rayIdx] = { false };

            point currentPoint = origin;
            int numReflections = 0;

            // Realizar el trazado de rayos mientras el número de reflexiones sea menor a 50
            while (numReflections < 50) {
                double minDistance = std::numeric_limits<double>::max();
                int intersectedPlane = -1;
                int intersectedTriangle = -1;
                int intersectedTriangleId = -1;

                // Iterar a través de cada plano en la escena
                for (int planeIdx = 0; planeIdx < NP; planeIdx++) {
                    plane plane = p[planeIdx];

                    // Comprobar si el rayo se encuentra en el lado correcto del plano
                    if (plane.n.x * ray.x + plane.n.y * ray.y + plane.n.z * ray.z < 0) {
                        double distance = IntersectionDistance(plane.n, plane.p[0], ray, currentPoint);

                        // Verificar si hay una intersección válida
                        if (distance > 0.000001 && distance < minDistance) {
                            point intersectionPoint;
                            intersectionPoint.x = currentPoint.x + ray.x * distance;
                            intersectionPoint.y = currentPoint.y + ray.y * distance;
                            intersectionPoint.z = currentPoint.z + ray.z * distance;

                            // Verificar si el punto de intersección está dentro de un triángulo
                            for (int triIdx = 0; triIdx < plane.NT; triIdx++) {
                                if (Inner(intersectionPoint, plane.t[triIdx])) {
                                    minDistance = distance;
                                    currentPoint = intersectionPoint;
                                    intersectedPlane = planeIdx;
                                    intersectedTriangle = triIdx;
                                    intersectedTriangleId = plane.t[triIdx].ID;
                                    break;
                                }
                            }
                        }
                    }
                }

                // Verificar si se encontró una intersección válida
                if (minDistance < std::numeric_limits<double>::max()) {
                    // Almacenar información sobre la reflexión
                    reflections[rayIdx].r[numReflections] = currentPoint;
                    reflections[rayIdx].d[numReflections] = minDistance;
                    reflections[rayIdx].idTriangle[numReflections] = intersectedTriangleId;
                    reflections[rayIdx].Plane[numReflections] = intersectedPlane;
                    reflections[rayIdx].Triangle[numReflections] = intersectedTriangle;
                    reflections[rayIdx].N = numReflections + 1;

                    // Calcular la dirección reflejada del rayo y actualizar el número de reflexiones
                    ray = Reflect(ray, p[intersectedPlane].n).Normal();
                    numReflections++;

                    // Verificar si se alcanzó el límite máximo de reflexiones
                    if (numReflections > 50) {
                        break;
                    }
                }
                else {
                    // Si no se encontró una intersección, marcar el rayo como perdido
                    reflections[rayIdx].lost = true;
                    break;
                }
            }
        }

        // Devolver el array de reflexiones
        return reflections;
    }

};
//---------------------------------------------------------------------------
class source {
public:
    point p;                //Posici n
    vector1* Rays;           //Direcciones de partida de los rayos
    int NRAYS;              //N mero de rayos
    double eF;              //Energ a inicial de la fuente
    triangle IcoFace[20];   //Representaci n gr fica de la fuente
    double VisualRadius;    //Tama o de la fuente (radio visual)

    source() {   //Inicializo las variables de la clase.
        p = 0.0;
        eF = 0.0;
        NRAYS = 0;
        Rays = NULL;

        //create icosaedron
        double S, R;
        point IcoVertex[12];

        //create vertexes
        S = 2 / sqrt(5);
        R = (5 - sqrt(5)) / 5;
        IcoVertex[0].x = 0;
        IcoVertex[0].y = 0;
        IcoVertex[0].z = 1;
        for (int i = 1; i < 6; i++) {
            IcoVertex[i].x = S * cos((PI * i * 72) / 180);
            IcoVertex[i].y = S * sin((PI * i * 72) / 180);
            IcoVertex[i].z = 1 - R;
            IcoVertex[i + 5].x = S * cos((72 * PI * i) / 180 + (36 * PI) / 180);
            IcoVertex[i + 5].y = S * sin((72 * PI * i) / 180 + (36 * PI) / 180);
            IcoVertex[i + 5].z = R - 1;
        }
        IcoVertex[11].x = 0;
        IcoVertex[11].y = 0;
        IcoVertex[11].z = -1;

        //create faces
        IcoFace[0].p0 = IcoVertex[0];   IcoFace[0].p1 = IcoVertex[1];   IcoFace[0].p2 = IcoVertex[2];
        IcoFace[1].p0 = IcoVertex[0];   IcoFace[1].p1 = IcoVertex[2];   IcoFace[1].p2 = IcoVertex[3];
        IcoFace[2].p0 = IcoVertex[0];   IcoFace[2].p1 = IcoVertex[3];   IcoFace[2].p2 = IcoVertex[4];
        IcoFace[3].p0 = IcoVertex[0];   IcoFace[3].p1 = IcoVertex[4];   IcoFace[3].p2 = IcoVertex[5];
        IcoFace[4].p0 = IcoVertex[0];   IcoFace[4].p1 = IcoVertex[5];   IcoFace[4].p2 = IcoVertex[1];
        IcoFace[5].p0 = IcoVertex[1];   IcoFace[5].p1 = IcoVertex[6];   IcoFace[5].p2 = IcoVertex[2];
        IcoFace[6].p0 = IcoVertex[2];   IcoFace[6].p1 = IcoVertex[6];   IcoFace[6].p2 = IcoVertex[7];
        IcoFace[7].p0 = IcoVertex[2];   IcoFace[7].p1 = IcoVertex[7];   IcoFace[7].p2 = IcoVertex[3];
        IcoFace[8].p0 = IcoVertex[3];   IcoFace[8].p1 = IcoVertex[7];   IcoFace[8].p2 = IcoVertex[8];
        IcoFace[9].p0 = IcoVertex[3];   IcoFace[9].p1 = IcoVertex[8];   IcoFace[9].p2 = IcoVertex[4];
        IcoFace[10].p0 = IcoVertex[4];  IcoFace[10].p1 = IcoVertex[8];  IcoFace[10].p2 = IcoVertex[9];
        IcoFace[11].p0 = IcoVertex[4];  IcoFace[11].p1 = IcoVertex[9];  IcoFace[11].p2 = IcoVertex[5];
        IcoFace[12].p0 = IcoVertex[5];  IcoFace[12].p1 = IcoVertex[9];  IcoFace[12].p2 = IcoVertex[10];
        IcoFace[13].p0 = IcoVertex[5];  IcoFace[13].p1 = IcoVertex[10]; IcoFace[13].p2 = IcoVertex[1];
        IcoFace[14].p0 = IcoVertex[1];  IcoFace[14].p1 = IcoVertex[10]; IcoFace[14].p2 = IcoVertex[6];
        IcoFace[15].p0 = IcoVertex[6];  IcoFace[15].p1 = IcoVertex[11]; IcoFace[15].p2 = IcoVertex[7];
        IcoFace[16].p0 = IcoVertex[7];  IcoFace[16].p1 = IcoVertex[11]; IcoFace[16].p2 = IcoVertex[8];
        IcoFace[17].p0 = IcoVertex[8];  IcoFace[17].p1 = IcoVertex[11]; IcoFace[17].p2 = IcoVertex[9];
        IcoFace[18].p0 = IcoVertex[9];  IcoFace[18].p1 = IcoVertex[11]; IcoFace[18].p2 = IcoVertex[10];
        IcoFace[19].p0 = IcoVertex[10]; IcoFace[19].p1 = IcoVertex[11]; IcoFace[19].p2 = IcoVertex[6];
    };

    void createRays(double NumberOfRays) {
        //matriz das Arestas {1o ponto da aresta, 2o ponto da aresta, Posi  o dos pontos da aresta na matriz Rays}
        int A[30][3] = { {0,1,0}, {0,2,0}, {0,3,0}, {0,4,0}, {0,5,0},
                        {1,6,0}, {2,6,0}, {2,7,0}, {3,7,0}, {3,8,0},
                        {4,8,0}, {4,9,0}, {5,9,0}, {5,10,0},{1,10,0},
                        {6,11,0},{7,11,0},{8,11,0},{9,11,0},{10,11,0},
                        {1,2,0}, {2,3,0}, {3,4,0}, {4,5,0}, {5,1,0},
                        {6,7,0}, {7,8,0}, {8,9,0}, {9,10,0},{10,6,0}
        };
        //matriz dos triangulos {1a aresta, 2a aresta, [0] V em p  [-1] V de cabe a pra baixo}
        int T[20][3] = { {0,1,0},   {1,2,0},   {2,3,0},   {3,4,0},   {4,0,0},
                        {5,6,-1},  {6,7,0},   {7,8,-1},  {8,9,0},   {9,10,-1},
                        {10,11,0}, {11,12,-1},{12,13,0}, {13,14,-1},{14,5,0},
                        {15,16,-1},{16,17,-1},{17,18,-1},{18,19,-1},{19,15,-1}
        };
        int i, j, k, n, m, RAY;
        double S, R, xB, yB, zB, xC, yC, zC, c[8];
        //create Rays matrix
        if (NRAYS > 0)
            delete[] Rays;
        n = int(floor(sqrt((NumberOfRays - 2) / 10) + 0.5));
        NRAYS = int(2 + 10 * pow(n, 2));
        Rays = new vector1[NRAYS];
        //calculating the icosaedron vertives
        S = 2 / sqrt(5);
        R = (5 - sqrt(5)) / 5;
        Rays[0].x = 0;
        Rays[0].y = 0;
        Rays[0].z = 1;
        for (i = 1; i < 6; i++) {
            Rays[i].x = S * cos((PI * i * 72) / 180);
            Rays[i].y = S * sin((PI * i * 72) / 180);
            Rays[i].z = 1 - R;
            Rays[i + 5].x = S * cos((72 * PI * i) / 180 + (36 * PI) / 180);
            Rays[i + 5].y = S * sin((72 * PI * i) / 180 + (36 * PI) / 180);
            Rays[i + 5].z = R - 1;
        }
        Rays[11].x = 0;
        Rays[11].y = 0;
        Rays[11].z = -1;
        RAY = 12;
        //calculating the rays on the icosaedron edges
        for (j = 0; j < 30; j++) {
            A[j][2] = RAY;
            xB = Rays[A[j][0]].x;
            yB = Rays[A[j][0]].y;
            zB = Rays[A[j][0]].z;
            xC = Rays[A[j][1]].x;
            yC = Rays[A[j][1]].y;
            zC = Rays[A[j][1]].z;
            c[0] = pow(xC, 2) * (pow(yB, 2) + pow(zB, 2)) + pow(yC * zB - yB * zC, 2) - 2 * xB * xC * (yB * yC + zB * zC) + pow(xB, 2) * (pow(yC, 2) + pow(zC, 2));
            c[1] = acos(xB * xC + yB * yC + zB * zC);
            c[2] = -xC * (yB * yC + zB * zC) + xB * (pow(yC, 2) + pow(zC, 2));
            c[3] = xC * (pow(yB, 2) + pow(zB, 2)) - xB * (yB * yC + zB * zC);
            c[4] = pow(xC, 2) * yB - xB * xC * yC + zC * (-yC * zB + yB * zC);
            c[5] = -xB * xC * yB + pow(xB, 2) * yC + zB * (yC * zB - yB * zC);
            c[6] = pow(xC, 2) * zB - xB * xC * zC + yC * (yC * zB - yB * zC);
            c[7] = -xB * xC * zB + pow(xB, 2) * zC + yB * (-yC * zB + yB * zC);
            for (i = 1; i < n; i++) {
                Rays[RAY].x = (c[2] * cos(i * c[1] / n) + c[3] * cos((n - i) * c[1] / n)) / c[0];
                Rays[RAY].y = (c[4] * cos(i * c[1] / n) + c[5] * cos((n - i) * c[1] / n)) / c[0];
                Rays[RAY].z = (c[6] * cos(i * c[1] / n) + c[7] * cos((n - i) * c[1] / n)) / c[0];
                RAY++;
            }
        }
        //calculating the rays on the icosaedron faces
        for (k = 0; k < 20; k++)
            for (j = 1; j < n; j++) {
                xB = Rays[A[T[k][0]][2] + j - 1].x;
                yB = Rays[A[T[k][0]][2] + j - 1].y;
                zB = Rays[A[T[k][0]][2] + j - 1].z;
                xC = Rays[A[T[k][1]][2] + j - 1].x;
                yC = Rays[A[T[k][1]][2] + j - 1].y;
                zC = Rays[A[T[k][1]][2] + j - 1].z;
                c[0] = pow(xC, 2) * (pow(yB, 2) + pow(zB, 2)) + pow(yC * zB - yB * zC, 2) - 2 * xB * xC * (yB * yC + zB * zC) + pow(xB, 2) * (pow(yC, 2) + pow(zC, 2));
                c[1] = acos(xB * xC + yB * yC + zB * zC);
                c[2] = -xC * (yB * yC + zB * zC) + xB * (pow(yC, 2) + pow(zC, 2));
                c[3] = xC * (pow(yB, 2) + pow(zB, 2)) - xB * (yB * yC + zB * zC);
                c[4] = pow(xC, 2) * yB - xB * xC * yC + zC * (-yC * zB + yB * zC);
                c[5] = -xB * xC * yB + pow(xB, 2) * yC + zB * (yC * zB - yB * zC);
                c[6] = pow(xC, 2) * zB - xB * xC * zC + yC * (yC * zB - yB * zC);
                c[7] = -xB * xC * zB + pow(xB, 2) * zC + yB * (-yC * zB + yB * zC);
                if (T[k][2] == 0)m = j;
                else m = n - j;
                for (i = 1; i < m; i++) {
                    Rays[RAY].x = (c[2] * cos(i * c[1] / m) + c[3] * cos((m - i) * c[1] / m)) / c[0];
                    Rays[RAY].y = (c[4] * cos(i * c[1] / m) + c[5] * cos((m - i) * c[1] / m)) / c[0];
                    Rays[RAY].z = (c[6] * cos(i * c[1] / m) + c[7] * cos((m - i) * c[1] / m)) / c[0];
                    RAY++;
                }
            }
    };
};
//---------------------------------------------------------------------------
class MatDouble {
public:
    double** A; // Puntero a un array 2D que contiene los elementos de la matriz.
    int n;// Número de filas de la matriz
    int m; // Número de columnas de la matriz.

    // Constructor
    MatDouble() {
        A = NULL; // Inicializa el puntero de la matriz a NULL.
        // Establece el número inicial de filas y columnas a 0.
        n = 0;
        m = 0;
    }

    // Método para inicializar la matriz con un tamaño específico.
    void Init(int i, int j) {
        // Establece el número de filas y columnas según los argumentos proporcionados.
        n = i;
        m = j;
        A = new double* [n]; // Reserva memoria para las filas de la matriz.
        for (int k = 0; k < n; k++) {
            A[k] = new double[m];// Reserva memoria para las columnas de la fila actual.
            for (int l = 0; l < m; l++) {
                A[k][l] = 0.0; // Inicializa el elemento actual de la matriz a 0.0.
            }
        }
    }
};


//---------------------------------------------------------------------------
class receptor {
public:
    point p;                // Posición del receptor en el espacio.
    double Radius;          // Radio del disco receptor.
    double eR;             // Energía recibida.

    // Constructor
    receptor() {
        p = 0.0;
        eR = NULL;
        Radius = 0.0;
    }
};

#endif
