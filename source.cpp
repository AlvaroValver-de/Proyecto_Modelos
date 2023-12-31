#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <learnopengl/shader.h>
#include <learnopengl/camera.h>
#define STB_IMAGE_IMPLEMENTATION 
#include <learnopengl/stb_image.h>
#include <iostream>
#include <iomanip> 
//==================================
#include "Header.h"
//VariablesGlobales
bool loadedRoom; //si la habitación ha sido cargada .
room r;          // Instancia para representar la habitación. 
MatDouble mD; //Matriz que almacene la distancia que existe entre los centros de los diferentes triángulos de la sala.
//Está matriz tiene una dimensión de 𝑛 × 𝑛 donde 𝑛 es el número de triángulos (NumTri) de la sala.
MatDouble mTV;//matriz que almacena los tiempos de vuelo entre baricentros en milisegundos que demorarían las reflexiones difusas
//para llegar de un centroide a otro(entre todos los triángulos que apliquen, es decir que sean visibles)
// La matriz tendrá una dimensión de 𝑛 × n
MatDouble mAS; //matriz angulos solidos
MatDouble mPE; //matriz porcentage de energia. La matriz tendrá una dimensión de 𝑛 × n
MatDouble mE; //matriz porcentage de energia. (espacio/tiempo) 
int NumTri = 0; //Almacena el número total de triángulos.
double N_RAYOS = 20; //Almacena el numero de rayos
int energiaFuente = 100; //Energia de la fuente
source s;           // Instancia de la fuente.
float alfa = 0.2; //Defincion del coeficiente alfa
float delta = 0.15; //Defincion del coeficiente alfa
int tM = 1000; //tiempo de simulacion en unidad discreta
double v_son = 340.0; //velocidad de transmisión de 340 m/s (constante V_SON)
reflection* arrayreflecciones = NULL; //Refleciones
point o; //Punto de origen
point** arrayRec;  //Array Recorrido
float elapsedTime = 0.0f;
const int NumReceptores = 27; // Número de receptores
MatDouble mDR; //Matriz Distancia Receptor
MatDouble mTVR; //Matriz Tiempo de Vuelo  Receptor
MatDouble mASR; //matriz angulos solidos Receptor
MatDouble mPER; //matriz porcentage de energia Receptor
MatDouble mER;

void laodRoom();
void energytransition();

//General Settings
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;
// camera
Camera camera(glm::vec3(0.0f, 0.0f, 10.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;
// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
//lighting
glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
int main()
{
    // glfw: initialize and configure
    glfwInit();
    //Version de Open gl
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef APPLE
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Proyecto Modelos y Simulacion - 2023 A", NULL, NULL);
    // Check si falla al crear
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    //Introduce th wwindow into the current context
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    // configure global opengl state
    // -----------------------------
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK); // Cull las caras traseras

    // build and compile shaders
    // -------------------------
    Shader room("shaders/room.vs", "shaders/room.fs");
    Shader rayo("shaders/lightcube.vs", "shaders/lightcube.fs");
    Shader receptorShader("shaders/room.vs", "shaders/room.fs");

    //Se carga la sala
    laodRoom();
    float vertices1[108];
    int contradork = 0;

    //El vertex imput para el cubo y la particula
    for (int i = 0; i < r.NP; i++) {
        for (int j = 0; j < r.p[i].NT; j++) {
            vertices1[contradork] = r.p[i].t[j].p0.x;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p0.y;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p0.z;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p1.x;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p1.y;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p1.z;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p2.x;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p2.y;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p2.z;
            contradork++;{}

        }
    }
    const int numSegments = 100;
    float vertices[numSegments * 3];

    for (int i = 0; i < numSegments; i++) {
        float angle = 2.0f * glm::pi<float>() * static_cast<float>(i) / static_cast<float>(numSegments);
        vertices[i * 3] = std::cos(angle);
        vertices[i * 3 + 1] = std::sin(angle);
        vertices[i * 3 + 2] = 0.0f;
    }

    energytransition();

    //Configuracion del Ambiente para el cubo y la particula
    // first, configure the cube's VAO (and VBO)
    unsigned int VBO, cubeVAO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices1), vertices1, GL_STATIC_DRAW);
    glBindVertexArray(cubeVAO);
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //Configuracion del cubo de luz
    unsigned int lightCubeVBO, lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glGenBuffers(1, &lightCubeVBO);
    glBindVertexArray(lightCubeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lightCubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices1), vertices1, GL_STATIC_DRAW);
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    unsigned int discoVBO, discoVAO;
    glGenVertexArrays(1, &discoVAO);
    glGenBuffers(1, &discoVBO);

    // Enlazar el VAO
    glBindVertexArray(discoVAO);

    // Enlazar y configurar el VBO
    glBindBuffer(GL_ARRAY_BUFFER, discoVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // Configurar atributos de vértice
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Desenlazar el VAO
    glBindVertexArray(0);

    int tSimulacion = 0;

    // Declaración de la matriz para registrar el número de impactos por triángulo
    int** impactosPorTriangulo = new int* [NumTri];
    for (int t = 0; t < NumTri; t++) {
        impactosPorTriangulo[t] = new int[tM];
        for (int ti = 0; ti < tM; ti++) {
            impactosPorTriangulo[t][ti] = 0; // Inicializa el contador de impactos a cero
        }
    }
    int* sumImpactosPorTriangulo = new int[NumTri]();
    // Difusión de la energía en los triángulos
    for (int r = 0; r < s.NRAYS; r++) {
        double distAcum = 0;
        int tiempo = 0;
        for (int re = 0; re < 50; re++) {
            int tri = arrayreflecciones[r].idTriangle[re];
            distAcum += arrayreflecciones[r].d[re];
            int tim = int(1000 * distAcum / v_son);
            impactosPorTriangulo[tri][tim] += 1;
            tiempo = tim;
        }
    }

    int* maxImpTri = new int[NumTri]();

    for (int i = 0; i < NumTri; i++) {
        for (int j = 1; j < tM; j++) {
            impactosPorTriangulo[i][j] = impactosPorTriangulo[i][j] + impactosPorTriangulo[i][j - 1];
                maxImpTri[i] = impactosPorTriangulo[i][j];
        }
    }
    // render loop
    while (!glfwWindowShouldClose(window))
    {
        // Actualiza el tiempo transcurrido
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // Actualiza el tiempo total transcurrido
        elapsedTime += deltaTime;

        if (elapsedTime >= 0.5) {
            if (tSimulacion < tM) {
                tSimulacion++;
            }
            elapsedTime = 0.0f;
        }

        // input
        processInput(window);
        // render
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Graficar el room       
       // view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

        glm::mat4 view = camera.GetViewMatrix();
        // world transformation
        glm::mat4 model = glm::mat4(1.0f);

        room.use();
        room.setVec3("lightColor", 0.0f, 1.0f, 0.0f);
        room.setVec3("lightPos", lightPos);
        room.setVec3("viewPos", camera.Position);
        
        int cont_t = 0; //contador del numero de triangulos
        for (int i = 0; i < r.NP; i++) {
            for (int j = 0; j < r.p[i].NT; j++) {

                if (mE.A[cont_t][tSimulacion] != 0.0) {
                    float energyValue = float(impactosPorTriangulo[cont_t][tSimulacion]) / 50.0;
                    // Configura el color calculado en el shader basado en el valor de energía
                    r.p[i].t[j].Color.heatMapColor(float(energyValue));
                }
                // Configura el color calculado en el shader
                room.setVec3("ourColor", glm::vec3(r.p[i].t[j].Color.red, r.p[i].t[j].Color.green, r.p[i].t[j].Color.blue));
                room.setMat4("projection", projection);
                room.setMat4("view", view);
                room.setMat4("model", model);
                // render the cube
                glBindVertexArray(cubeVAO);
                glDrawArrays(GL_TRIANGLES, cont_t * 3, 3);
                cont_t++;
            }
        }

        //Dibujar la fuente
        rayo.use();
        rayo.setMat4("projection", projection);
        rayo.setMat4("view", view);
        rayo.setVec3("ourColor", glm::vec3(0.75f, 0.75f, 0.0f));
        model = glm::mat4(1.0f);

        for (int idRayo = 0; idRayo < s.NRAYS; idRayo++) {
            // Calculate the model matrix for each ray instance
            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, glm::vec3(arrayRec[idRayo][tSimulacion].x,
                arrayRec[idRayo][tSimulacion].y,
                arrayRec[idRayo][tSimulacion].z));

            model = glm::scale(model, glm::vec3(0.025f)); // Tamaño del cubo pequeño
            rayo.setMat4("model", model);

            // Render the cube for the current ray
            glBindVertexArray(lightCubeVAO);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }
        // Renderizar el disco de los receptores
        rayo.use();
        rayo.setMat4("projection", projection);
        rayo.setMat4("view", view);
        for (int i = 0; i < NumReceptores; i++)
        {
            rayo.setVec3("ourColor", glm::vec3(r.r[i].Color.red, r.r[i].Color.green, r.r[i].Color.blue));
            model = glm::mat4(1.0f); // Restablecer la matriz de modelo
            model = glm::translate(model, glm::vec3(float(r.r[i].p.x), float(r.r[i].p.y), float(r.r[i].p.z)));
            model = glm::scale(model, glm::vec3(0.25f));
            rayo.setMat4("model", model);
            glBindVertexArray(discoVAO); // Enlazar el VAO del disco
            glDrawArrays(GL_TRIANGLE_FAN, 0, numSegments); // Dibujar el disco
            glBindVertexArray(0); // Desenlazar el VAO del disco
           
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    // optional: de-allocate all resources once they've outlived their purpose:
// ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteVertexArrays(1, &lightCubeVAO);
    glDeleteBuffers(1, &VBO);
    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;

}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);

}
// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
    glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
}
// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }
    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;
    camera.ProcessMouseMovement(xoffset, yoffset);
}
// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}
//================================================================================================================================================
// Función para para calculos al cargar la habitación:
void laodRoom() {
    if (!loadedRoom) {   // Verifica si la habitación ya ha sido cargada anteriormente.
        r.NewPlanes(6);// Genearra 6 planos
        //square back
        r.p[0].NewPoints(4);// Gnererar los 4 puntos
        r.p[0].p[0].x = -2.0f;
        r.p[0].p[0].y = 2.0f;
        r.p[0].p[0].z = 2.0f;
        r.p[0].p[1].x = -2.0f;
        r.p[0].p[1].y = -2.0f;
        r.p[0].p[1].z = 2.0f;
        r.p[0].p[2].x = -2.0f;
        r.p[0].p[2].y = -2.0f;
        r.p[0].p[2].z = -2.0f;
        r.p[0].p[3].x = -2.0f;
        r.p[0].p[3].y = 2.0f;
        r.p[0].p[3].z = -2.0f;
        r.p[0].PointGenTriangle();
        //square front
        r.p[1].NewPoints(4);// Gnererar los 4 puntos
        r.p[1].p[0].x = 2.0f;
        r.p[1].p[0].y = 2.0f;
        r.p[1].p[0].z = 2.0f;
        r.p[1].p[1].x = 2.0f;
        r.p[1].p[1].y = 2.0f;
        r.p[1].p[1].z = -2.0f;
        r.p[1].p[2].x = 2.0f;
        r.p[1].p[2].y = -2.0f;
        r.p[1].p[2].z = -2.0f;
        r.p[1].p[3].x = 2.0f;
        r.p[1].p[3].y = -2.0f;
        r.p[1].p[3].z = 2.0f;
        r.p[1].PointGenTriangle();
        //square left
        r.p[2].NewPoints(4);
        r.p[2].p[0].x = -2.0f;
        r.p[2].p[0].y = -2.0f;
        r.p[2].p[0].z = 2.0f;
        r.p[2].p[1].x = 2.0f;
        r.p[2].p[1].y = -2.0f;
        r.p[2].p[1].z = 2.0f;
        r.p[2].p[2].x = 2.0f;
        r.p[2].p[2].y = -2.0f;
        r.p[2].p[2].z = -2.0f;
        r.p[2].p[3].x = -2.0f;
        r.p[2].p[3].y = -2.0f;
        r.p[2].p[3].z = -2.0f;
        r.p[2].PointGenTriangle();
        //square right
        r.p[3].NewPoints(4);// Gnererar los 4 puntos
        r.p[3].p[0].x = 2.0f;
        r.p[3].p[0].y = 2.0f;
        r.p[3].p[0].z = 2.0f;
        r.p[3].p[1].x = -2.0f;
        r.p[3].p[1].y = 2.0f;
        r.p[3].p[1].z = 2.0f;
        r.p[3].p[2].x = -2.0f;
        r.p[3].p[2].y = 2.0f;
        r.p[3].p[2].z = -2.0f;
        r.p[3].p[3].x = 2.0f;
        r.p[3].p[3].y = 2.0f;
        r.p[3].p[3].z = -2.0f;
        r.p[3].PointGenTriangle();
        //square top
        r.p[4].NewPoints(4);
        r.p[4].p[0].x = -2.0f;
        r.p[4].p[0].y = -2.0f;
        r.p[4].p[0].z = 2.0f;
        r.p[4].p[1].x = -2.0f;
        r.p[4].p[1].y = 2.0f;
        r.p[4].p[1].z = 2.0f;
        r.p[4].p[2].x = 2.0f;
        r.p[4].p[2].y = 2.0f;
        r.p[4].p[2].z = 2.0f;
        r.p[4].p[3].x = 2.0f;
        r.p[4].p[3].y = -2.0f;
        r.p[4].p[3].z = 2.0f;
        r.p[4].PointGenTriangle();
        //square bottom
        r.p[5].NewPoints(4);
        r.p[5].p[0].x = -2.0f;
        r.p[5].p[0].y = 2.0f;
        r.p[5].p[0].z = -2.0f;
        r.p[5].p[1].x = -2.0f;
        r.p[5].p[1].y = -2.0f;
        r.p[5].p[1].z = -2.0f;
        r.p[5].p[2].x = 2.0f;
        r.p[5].p[2].y = -2.0f;
        r.p[5].p[2].z = -2.0f;
        r.p[5].p[3].x = 2.0f;
        r.p[5].p[3].y = 2.0f;
        r.p[5].p[3].z = -2.0f;
        r.p[5].PointGenTriangle();

        int cont_t = 0; //contador del numero de triangulos
        // Loop para calcular el centroide de todos los triángulos de la sala
        for (int i = 0; i < r.NP; i++) {                 // Recorre los planos de la sala.
            r.p[i].n = r.p[i].NormalPlano();       // Calcula la normal del plano.
            for (int j = 0; j < r.p[i].NT; j++) {        // Recorre los triángulos del plano.
                r.p[i].t[j].CalculateProjection();       // Calcula la proyección del triángulo.
                r.p[i].t[j].Centroid();                  // Calcula el centroide baricentro del triángulo.
                r.p[i].t[j].ID = cont_t;                 // Asigna un ID al triángulo.
                cont_t++;                                // Incrementa el contador de triángulos.
            }
        }
        NumTri = cont_t;  // Asigna el número total de triángulos.

        //Inicializacion de las matrices
        mD.Init(NumTri, NumTri);
        mTV.Init(NumTri, NumTri);;
        mAS.Init(NumTri, NumTri);
        mPE.Init(NumTri, NumTri);

        int cont = 0;
        double* suma_angulos_solidos = new double[NumTri](); //arreglo para optener la sunma de las areas de los angulos solidos
        //Ciclo para el calculo de las matriz de distacia, tiempo de vuelo y angulos solidos
        for (int i = 0; i < r.NP; i++) {                      // Recorre los planos de la sala.
            for (int j = 0; j < r.p[i].NT; j++) {             // Recorre los triángulos del plano.
                int idTri1 = r.p[i].t[j].ID;                  //Obtiene el id del triangulo 1
                for (int k = 0; k < r.NP; k++) {              //Se reccorre de nuevo los planos de la sala.
                    for (int l = 0; l < r.p[k].NT; l++) {     //Recorre los triángulos del plano.
                        int idTri2 = r.p[k].t[l].ID;          //Obtiene el id del triangulo 2
                        if (i != k) {                    //verifica que ambos triángulos no pertenezcan al mismo plano.
                            mD.A[idTri1][idTri2] = r.p[i].t[j].bc.distancia(r.p[k].t[l].bc); //Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mTV.A[idTri1][idTri2] = int(1000 * mD.A[idTri1][idTri2] / V_SON); // Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mAS.A[idTri1][idTri2] = r.p[k].t[l].AnguloSolido(r.p[i].t[j].bc); // Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            suma_angulos_solidos[cont] += mAS.A[idTri1][idTri2];
                        }
                    }
                }
                cont++;
            }
        }
        //Calculo de la matriz porcentaje de energia
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                mPE.A[i][j] = mAS.A[i][j] / suma_angulos_solidos[i];
            }
        }

        std::cout << "Matriz distancia entre Baricentros de los triangulos" << std::endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << std::fixed << std::setprecision(2) << mD.A[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Matriz tiempo de vuelo entre Baricentros de los triangulos" << std::endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << mTV.A[i][j] << " "; // Use setw for formatting integers
            }
            std::cout << std::endl;
        }

        std::cout << "Matriz porcentaje de energia de los triangulos" << std::endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << std::fixed << std::setprecision(3) << mPE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        //libera la memoria
        delete[] suma_angulos_solidos;

        // CREACIÓN DE RECEPTORES 
        r.NewReceptor(NumReceptores);
        int cont_rec = 0;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {
                    r.r[cont_rec].p.x = float(i);
                    r.r[cont_rec].p.y = float(j);
                    r.r[cont_rec].p.z = float(k);
                    cont_rec++;
                }
            }
        }
        // Inicialización de las matrices
        mDR.Init(r.NR, NumTri);
        mTVR.Init(r.NR, NumTri);
        mASR.Init(r.NR, NumTri);
        mPER.Init(r.NR, NumTri);

        double* suma_angulos_solidos_receptor = new double[r.NR]();

        // Ciclo para el cálculo de las matrices de distancia, tiempo de vuelo y ángulos sólidos de los receptores
        for (int i = 0; i < r.NP; i++) {                      // Recorre los planos de la sala.
            for (int j = 0; j < r.p[i].NT; j++) {             // Recorre los triángulos del plano.
                int idTri = r.p[i].t[j].ID;
                // Cálculo de las distancias, tiempo de vuelo y ángulos sólidos de los receptores
                for (int m = 0; m < r.NR; m++) {
                    mDR.A[m][idTri] = r.r[m].p.distancia(r.p[i].t[j].bc); // Cálculo de la distancia de la posición del receptor al baricentro de los triángulos de la sala
                    mTVR.A[m][idTri] = int(1000 * mDR.A[m][idTri] / v_son); // Cálculo del tiempo de vuelo entre el receptor y el baricentro de los triángulos de la sala
                    mASR.A[m][idTri] = r.r[m].CalcularAreaDiscoProyectado(r.p[i].t[j].bc); // Cálculo del ángulo sólido entre el receptor y el baricentro de los triángulos de la sala
                    suma_angulos_solidos_receptor[m] += mASR.A[m][idTri];
                }
            }
        }

        // Cálculo de la matriz porcentaje de energía del receptor
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                mPER.A[i][j] = mASR.A[i][j] / suma_angulos_solidos_receptor[i];
            }
        }

        std::cout << "Matriz distancia entre Baricentros de los receptores" << std::endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << std::fixed << std::setprecision(2) << mDR.A[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Matriz tiempo de vuelo entre Baricentros de los receptores" << std::endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << mTVR.A[i][j] << " "; // Use setw for formatting integers
            }
            std::cout << std::endl;
        }

        std::cout << "Matriz porcentaje de energia de los receptores" << std::endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                std::cout << std::fixed << std::setprecision(3) << mPER.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        loadedRoom = true; // Indica que la habitación ha sido cargada exitosamente.
    }
}
void energytransition() {
    if (loadedRoom) {
        int t_vuelo = 0; // Tiempo de vuelo del rayo

        o.x = 1.0;
        o.y = 1.0;
        o.z = 1.0;

        s.eF = energiaFuente; //Energia de la fuente

        s.createRays(N_RAYOS);  //Numero de rayos

        double eneRayo = s.eF / s.NRAYS; // Energía Inicial

        arrayreflecciones = r.RayTracing(o, s.Rays, s.NRAYS);

        mE.Init(NumTri, tM); //inicializacion de la matriz de energia

        mER.Init(r.NR, tM); //inicializacion de la matriz de energia del receptor

        arrayRec = new point * [s.NRAYS];

        // Difusión de la energía difusa en los triángulos
        for (int rayo = 0; rayo < s.NRAYS; rayo++) {
            double eneResidual = eneRayo;
            double distAcum = 0;
            arrayRec[rayo] = new point[tM];
            arrayRec[rayo][0] = o;
            int tiempo = 0;
            for (int re = 0; re < arrayreflecciones[rayo].N; re++) {
                int tri = arrayreflecciones[rayo].idTriangle[re];
                point pun = arrayreflecciones[rayo].r[re];
                distAcum += arrayreflecciones[rayo].d[re];
                int tim = int(1000 * distAcum / v_son);
                arrayRec[rayo][tim] = pun;
                for (int t = tiempo + 1; t < tim; t++) {
                    arrayRec[rayo][t].x = arrayRec[rayo][t - 1].x + (arrayRec[rayo][tim].x - arrayRec[rayo][tiempo].x) / (tim - tiempo);
                    arrayRec[rayo][t].y = arrayRec[rayo][t - 1].y + (arrayRec[rayo][tim].y - arrayRec[rayo][tiempo].y) / (tim - tiempo);
                    arrayRec[rayo][t].z = arrayRec[rayo][t - 1].z + (arrayRec[rayo][tim].z - arrayRec[rayo][tiempo].z) / (tim - tiempo);
                }
                for (int j = 0; j < r.NR; j++) {
                    mER.A[tri][tim] = eneResidual; //Energia  del rayo al receptor
                }
                tiempo = tim;
                mE.A[tri][tim] += (eneResidual * (1 - alfa) * delta); //Energia difusa en los triangulos
                eneResidual = eneResidual * (1 - alfa) * (1 - delta); //Energia incidente de los rayos
            }
        }

        //Transicion de energia en la matriz esapcio tiempo mE
        for (int t = 0; t < tM; t++) {
            for (int e = 0; e < NumTri; e++) {// Triángulo 1
                for (int ed = 0; ed < NumTri; ed++) { // Triángulo 2
                    //Energia de la los triangulos de la sala
                    if (e != ed) { //No estan en el mismo plano
                        //Es el instante de tiempo de la simulacion + el instate de  tiempo en que ocurre la transmision de energia en ese triangulo a hacia el siguiente triangulo
                        t_vuelo = mTV.A[e][ed] + t;
                        if (t_vuelo <= tM) {
                            mE.A[ed][t_vuelo] += (mE.A[e][t] * mPE.A[e][ed]) * (1 - alfa);
                        }
                    }
                }
                //Energía del Receptor
                for (int k = 0; k < r.NR; k++) {
                    t_vuelo = mTVR.A[k][e] + t;
                    if (t_vuelo < tM) {
                        mER.A[k][t_vuelo] += (mE.A[e][t] * mPER.A[k][e]);
                    }
                };
            }
        }
        std::cout << "Matriz de energia" << std::endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < 15; j++) {
                std::cout << mE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Matriz de energia del receptor" << std::endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < 15; j++) {
                std::cout << mER.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
}
