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
//==================================
#include "Header.h"

//VariablesGlobales

bool loadedRoom;
room r;
MatDouble mD; //matriz distancia entre baricentros
MatDouble mTV;//matriz tiempo de vuelo entre baricentros
MatDouble mAS; //matriz angulos solidos
MatDouble mPE; //matriz porcentage de energia
MatDouble mET; //matriz porcentage de energia

int NumTri = 0;
source s;
float speed = 0.3; //Velocidad de movimiento de la animacion
float alfa = 0.5;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void laodRoom();


// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
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
    // ------------------------------
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
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);


    // build and compile shaders
    // -------------------------
    Shader room("shaders/room.vs", "shaders/room.fs");
    Shader ico("shaders/lightcube.vs", "shaders/lightcube.fs");
    Shader rayo("shaders/lightcube.vs", "shaders/lightcube.fs");


    //Se carga la sala
    laodRoom();


    int numeroTriangulos = NumTri;
    float vertices1[108];
    int contradork = 0;
    float vertices2[180];
    int contadorIco = 0;


    //El vertex imput para el cubo y la particula
    for (int i = 0; i < r.NP; i++) {
        for (int j = 0; j < r.p[i].NT; j++) {
            vertices1[contradork] = r.p[i].t[j].p0.x * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p0.y * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p0.z * 0.1;
            contradork++;


            vertices1[contradork] = r.p[i].t[j].p1.x * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p1.y * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p1.z * 0.1;
            contradork++;


            vertices1[contradork] = r.p[i].t[j].p2.x * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p2.y * 0.1;
            contradork++;
            vertices1[contradork] = r.p[i].t[j].p2.z * 0.1;
            contradork++;

        }
    }

    //el verrtex input para el icocaedro
    contadorIco = 0;
    for (int i = 0; i < 20; i++) {
        vertices2[contadorIco] = s.IcoFace[i].p0.x * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p0.y * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p0.z * 0.1;
        contadorIco++;


        vertices2[contadorIco] = s.IcoFace[i].p1.x * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p1.y * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p1.z * 0.1;
        contadorIco++;


        vertices2[contadorIco] = s.IcoFace[i].p2.x * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p2.y * 0.1;
        contadorIco++;
        vertices2[contadorIco] = s.IcoFace[i].p2.z * 0.1;
        contadorIco++;


    }

    point o;
    o.x = -1.0;
    o.y = 1.5;
    o.z = 1.0;


    s.createRays(20);


    for (int i = 0; i < 20; i++) {
        printf("Rayos: x: %f, y: %f, z: %f\n", s.Rays[i].x, s.Rays[i].y, s.Rays[i].z);
    }



    reflection* arrayreflecciones;

    arrayreflecciones = r.RayTracing(o, s.Rays, s.NRAYS);

    for (int i = 0; i < 50; i++) {
        printf(" %d - Punto de Incidencia: x: %f, y: %f, z: %f\n", i, arrayreflecciones[1].r[i].x, arrayreflecciones[1].r[i].y, arrayreflecciones[1].r[i].z);
    }

    reflection arrayDePuntosDeIncidencia = arrayreflecciones[1];

    mET.Init(NumTri, 1000);

    //Transicion de energia
    for (int t = 0; t < 100; t++) {
        for (int e = 0; e < NumTri; e++) {
            for (int ed = 0; ed < NumTri; ed++) {
                if (e != ed && mPE.A[e][ed] != 0) {
                    int tv = mTV.A[e][ed];
                    if ((t + tv) <= 10) {
                        mET.A[t + tv][ed] += mET.A[t][e] * mPE.A[e][ed]*(1- alfa);
                    }
                }
            }
        }
    }


        for (int i = 0; i < mD.m; i++) {
            for (int j = 0; j < mD.n; j++) {
                std::cout << mET.A[i][j] << " ";
            }
            std::cout << std::endl;
        }


    int nPunto = 0;


    int idRayo = 8;




    //Punto de Origen
    point puntoDeOrigen;

    puntoDeOrigen.x = arrayreflecciones[idRayo].r[nPunto].x * 0.1;
    puntoDeOrigen.y = arrayreflecciones[idRayo].r[nPunto].y * 0.1;
    puntoDeOrigen.z = arrayreflecciones[idRayo].r[nPunto].z * 0.1;

    //Punto de Origen
    point puntoDeIncidencia;
    nPunto++;
    puntoDeIncidencia.x = arrayreflecciones[idRayo].r[nPunto].x * 0.1;
    puntoDeIncidencia.y = arrayreflecciones[idRayo].r[nPunto].y * 0.1;
    puntoDeIncidencia.z = arrayreflecciones[idRayo].r[nPunto].z * 0.1;






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

    //Configuracion del ico
    unsigned int VBO2, cubeVAO2;
    glGenVertexArrays(1, &cubeVAO2);
    glGenBuffers(1, &VBO2);

    glBindBuffer(GL_ARRAY_BUFFER, VBO2);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices2), vertices2, GL_STATIC_DRAW);


    glBindVertexArray(cubeVAO2);

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



    double tiempo1 = 0;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
       // --------------------
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
// ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Graficar el room
          // be sure to activate shader when setting uniforms/drawing objects
        room.use();
        room.setVec3("ourColor", 0.0f, 0.0f, 1.0f);
        room.setVec3("lightColor", 0.0f, 1.0f, 0.0f);
        room.setVec3("lightPos", lightPos);
        room.setVec3("viewPos", camera.Position);
        // view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();


        room.setMat4("projection", projection);
        room.setMat4("view", view);

        // world transformation
        glm::mat4 model = glm::mat4(1.0f);
        room.setMat4("model", model);

        // render the cube
        glBindVertexArray(cubeVAO);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_TRIANGLES, 0, 36);

        //Render del icosaedro
        ico.use();
        ico.setVec3("objectColor", 1.0f, 0.0f, 0.0f);
        ico.setVec3("lightColor", 1.0f, 1.0f, 1.0f);


        ico.setVec3("lightPos", lightPos);


        ico.setVec3("viewPos", camera.Position);

        // view/projection transformations
        glm::mat4 projection2 = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view2 = camera.GetViewMatrix();


        ico.setMat4("projection", projection2);
        ico.setMat4("view", view2);

        // world transformation
        glm::mat4 model2 = glm::mat4(1.0f);
        ico.setMat4("model", model2);


        // render the cube
        glBindVertexArray(cubeVAO2);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_TRIANGLES, 0, 60);


        //Dibujar la fuente
        rayo.use();
        rayo.setMat4("projection", projection);
        rayo.setMat4("view", view);
        rayo.setVec3("ourColor", 0.5f, 1.0f, 0.5f);

        model = glm::mat4(1.0f);

        float distancia = puntoDeOrigen.distancia(puntoDeIncidencia);

        if ((puntoDeOrigen.distancia(puntoDeIncidencia) * (glfwGetTime() - tiempo1) * speed) >= distancia) {

            tiempo1 = glfwGetTime();
            puntoDeOrigen = puntoDeIncidencia;
            nPunto++;
            puntoDeIncidencia.x = arrayreflecciones[idRayo].r[nPunto].x * 0.1;
            puntoDeIncidencia.y = arrayreflecciones[idRayo].r[nPunto].y * 0.1;
            puntoDeIncidencia.z = arrayreflecciones[idRayo].r[nPunto].z * 0.1;
        };

        vector1 mov = ((puntoDeIncidencia - puntoDeOrigen)) * (glfwGetTime() - tiempo1) * speed;

        model = glm::translate(model, glm::vec3(puntoDeOrigen.x + mov.x, puntoDeOrigen.y + mov.y, puntoDeOrigen.z + mov.z));

        model = glm::scale(model, glm::vec3(0.01f)); // tamaño del cubo pequeño


        rayo.setMat4("model", model);
        glBindVertexArray(lightCubeVAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);




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

    //If I want to stay in ground level (xz plane)
    //camera.Position.y = 0.0f;

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

void laodRoom() {
    if (!loadedRoom) {
        int nDivTri;
        nDivTri = 2;
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


        //Calcular los normales del plano
        int cont_t = 0;
        for (int i = 0; i < r.NP; i++) {
            r.p[i].n = r.p[i].NormalPlano(r.p[i]);
            for (int j = 0; j < r.p[i].NT; j++) {
                r.p[i].t[j].CalculateProjection();
                r.p[i].t[j].Centroid();
                r.p[i].t[j].ID = cont_t;
                cont_t++;
            }
        }

        NumTri = cont_t;



        mD.Init(NumTri, NumTri);
        //mD.A[0][0] = 0.0;

        mTV.Init(NumTri, NumTri);
        //mTV.A[0][0] = 0.0;


        mAS.Init(NumTri, NumTri);
        mPE.Init(NumTri, NumTri);
        //mAS.A[0][0] = 0.0;
        
        
        int cont = 0;
        double* suma_angulos_solidos = new double[NumTri]();

        
        for (int i = 0; i < r.NP; i++) {

            for (int j = 0; j < r.p[i].NT; j++) {
                // TRIANGULO 1 filas
                int idTri1 = r.p[i].t[j].ID;
                
               
                for (int k = 0; k < r.NP; k++) {

                    for (int l = 0; l < r.p[k].NT; l++) {
                        // TRIANGULO 2
                        int idTri2 = r.p[k].t[l].ID;

                        if (i != k) { // PARTE 6
                            mD.A[idTri1][idTri2] = r.p[i].t[j].bc.distancia(r.p[k].t[l].bc); //Cálculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mTV.A[idTri1][idTri2] = int(1000 * mD.A[idTri1][idTri2] / V_SON); // Cálculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mAS.A[idTri1][idTri2] = r.p[k].t[l].AnguloSolido(r.p[i].t[j].bc); // Cálculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            suma_angulos_solidos[cont] = abs(mAS.A[idTri1][idTri2]);

                            // areaT[idTri1] += matRoomAngles.d[idTri1][idTri2];


                        }

                    }
                }
                cont++;

            }
        }


        /*
        for (int i =0; i<)

        for (int i = 0; i < mD.m; i++) {
            for (int j = 0; j < mD.n; j++) {
                std::cout <<mPE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
*/
   
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                mPE.A[i][j] = mAS.A[i][j] / suma_angulos_solidos[i];
            }
        }

        
        
        for (int i = 0; i < mD.m; i++) {
            for (int j = 0; j < mD.n; j++) {
                std::cout << mPE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        
    }


}