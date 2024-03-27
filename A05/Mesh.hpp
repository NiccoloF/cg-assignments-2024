

/**************
 Creae the meshes, as described below

 WARNING!
 Since it is a C program, you can use for loops and functions if you think they can be helpful in your solution.
 However, please include all your code in this file, since it will be put in an automatic correction process
 for the final evaluation. Please also be cautious when using standard libraries and symbols, since they
 might not be available in all the development environments (especially, they might not be available
 in the final evaluation environment, preventing your code from compiling).
 This WARNING will be valid far ALL THE ASSIGNMENTs, but it will not be repeated in the following texts,
 so please remember these advices carefully!

***************/

void MakeSquare(float size, std::vector<glm::vec3> &vertices, std::vector<uint32_t> &indices) {
// Creates a square, on the xz-plane, aligned with the axis, and centered in the origin.
// The length of the four sides is in parameter >size<.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a glm::vec3 element.
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: this procedure has already been implemented. You can keep it as is
	vertices = {
				   {-size/2.0f,0.0f,-size/2.0f},
				   {-size/2.0f,0.0f, size/2.0f},
				   { size/2.0f,0.0f,-size/2.0f},
				   { size/2.0f,0.0f, size/2.0f}};
	indices = {0, 1, 2,    1, 3, 2};

}

void MakeCube(float size, std::vector<glm::vec3> &vertices, std::vector<uint32_t> &indices) {
// Creates a cube, with the faces perpendicular to the axis, and centered in the origin.
// The length of one edge of the cube is >size<.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a glm::vec3 element.
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a square. You can use it as a side of the cube (please remember
// to change the value of the y component, otherwise the result will be wrong
  vertices.clear();
  indices.clear();
	vertices = {
				   {-size/2.0f,-size/2.0f,-size/2.0f},
				   {-size/2.0f,-size/2.0f,size/2.0f},
				   {-size/2.0f,size/2.0f,size/2.0f},
				   {-size/2.0f,size/2.0f,-size/2.0f},
           {size/2.0f,-size/2.0f,-size/2.0f},
           {size/2.0f,-size/2.0f,size/2.0f},
           {size/2.0f,size/2.0f,size/2.0f},
           {size/2.0f,size/2.0f,-size/2.0f}};

	indices = {0,1,2,  2,3,0,  0,3,4,  4,3,7,
             6,5,7,  4,7,5,  6,1,5,  2,1,6,
             6,7,2,  3,2,7,  4,5,0,  1,0,5};


}

void MakeCylinder(float radius, float height, int slices, std::vector<glm::vec3> &vertices, std::vector<uint32_t> &indices) {

  vertices.clear();
  indices.clear();
  vertices.resize(2 * slices + 2); // Aggiungiamo 2 vertici per i cerchi superiori e inferiori
  indices.resize(12*slices); // Ogni cerchio richiede slices * 3 indici

  // Aggiungi vertici per i cerchi superiori e inferiori
  vertices[2 * slices] = glm::vec3(0.0f, height/2, 0.0f); // Vertice centrale del cerchio superiore
  vertices[2 * slices + 1] = glm::vec3(0.0f, -height/2, 0.0f); // Vertice centrale del cerchio inferiore

  for (int i = 0; i < slices; ++i) {
      // Create vertices for the side faces
      float ang = 2 * M_PI * (float)i / (float)slices;
      vertices[i] = glm::vec3(radius * cos(ang), height/2, radius * sin(ang)); // Vertici del cerchio superiore
      vertices[i + slices] = glm::vec3(radius * cos(ang), -height/2, radius * sin(ang)); // Vertici del cerchio inferiore


      // top circle
      indices[3*i] = (i+1)%slices;
      indices[3*i+1] = i;
      indices[3*i+2] = 2*slices;

      //bottom circle
      indices[3*(i + slices)] = i + slices;
      indices[3*(i+slices)+1] = (i+1)%slices + slices;
      indices[3*(i+slices)+2] = 2*slices + 1;

      // Vertici del lato superiore del cilindro
      indices[6 * slices + 3 * i] = i;
      indices[6 * slices + 3 * i+1] = (i + 1) % slices;
      indices[6 * slices + 3 * i+2] = i + slices;

      // Vertici del lato inferiore del cilindro
      indices[9 * slices + 3 * i] = (i + 1) % slices;
      indices[9 * slices + 3 * i+1] = (i + 1) % slices + slices;
      indices[9 * slices + 3 * i+2] = i + slices;

  }

}

void MakeCone(float radius, float height, int slices, std::vector<glm::vec3> &vertices, std::vector<uint32_t> &indices) {
// Creates a cone, approximated by a pyramid with a base composed by a regular polygon with >slices< sides.
// The radius of the cone is >radius<, and it height is >height<. The cone has its centere
// in the origin, and spans half above and half below the plane that passes thorugh the origin.
// The bottom cap is aligned with xz-plane and perpendicular to the y-axis.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a glm::vec3 element.
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a triangle. You have to change it, or you will obtain a wrong result
// You should use a for loop, and you should start from the procedure to create a circle seen during the lesson
// Generate vertices for the base of the cone

    vertices.clear();
    indices.clear();
    vertices.resize(slices+2);
    indices.resize(6*slices);

    vertices[slices]= {0.0f,-height/2.0f,0.0f};
    vertices[slices+1]= {0.0f,height/2.0f,0.0f};
    for(int i = 0; i < slices; i++) {
      float ang = 2*M_PI * (float)i / (float)slices;
      vertices[i] = glm::vec3(radius * cos(ang),-height/2.0f,radius * sin(ang));
      indices[3*i  ] = slices;
      indices[3*i+1] = i;
      indices[3*i+2] = (i+1) % slices;

      indices[3*(i+slices)  ] = i;
      indices[3*(i+slices)+1] = slices+1;
      indices[3*(i+slices)+2] = (i+1) % slices;
    }
}

void MakeSphere(float radius, int rings, int slices, std::vector<glm::vec3> &vertices, std::vector<uint32_t> &indices)
{
// Creates a sphere, approximated by a poliedron composed by >slices<, and >rings< rings.
// The radius of the sphere is >radius<, and it is centered in the origin.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a glm::vec3 element.
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a circle. You have to change it, or you will obtain a wrong result
// You should use two nested for loops, one used to span across the rings, and the other that spans along
// the rings.

    // Clear any existing data
    vertices.clear();
    indices.clear();
    vertices.resize(slices*rings+2);
    indices.resize(6*slices*rings);

    // Generate vertices
    vertices[0] = glm::vec3(0.0f, 1.0f, 0.0f) * radius;
    for (int i = 1; i <= rings; ++i) {
        float theta = M_PI * (float)i / ((float)rings + 1.0f);
        float sinTheta = glm::sin(theta);
        float cosTheta = glm::cos(theta);

        for (int j = 0; j < slices; ++j) {
            float phi = 2*M_PI * (float)j / (float)slices;
            float sinPhi = glm::sin(phi);
            float cosPhi = glm::cos(phi);

            float x = cosPhi * sinTheta;
            float y = cosTheta;
            float z = sinPhi * sinTheta;

            vertices[(i-1)*slices+j+1] = glm::vec3(x, y, z) * radius;
        }
    }
    vertices[slices*rings+1] = glm::vec3(0.0f, -1.0f, 0.0f) * radius;

    // Generate polo nord
    for (int j = 0; j < slices; ++j) {
        int firstPole = j+1;
        int secondPole = (j + 1)%slices+1;

        indices[6*j] = 0; // Polo nord
        indices[6*j+1] = secondPole;
        indices[6*j+2] = firstPole;

        firstPole = (rings-1)*slices + j+1;
        secondPole = (j + 1)%slices + (rings-1)*slices+1;

        indices[6*j+3] = slices*rings + 1;
        indices[6*j+4] = firstPole;
        indices[6*j+5] = secondPole;

    }

    // Generate surface triangles
for (int i = 0; i < rings - 1; ++i) {
    for (int j = 0; j < slices; ++j) {
        int currentRingStart = i * slices + 1;
        int nextRingStart = (i + 1) * slices + 1;

        int first = currentRingStart + j;
        int second = currentRingStart + (j + 1) % slices;
        int third = nextRingStart + j;
        int fourth = nextRingStart + (j + 1) % slices;

        // First triangle
        indices[6*slices + 6 * (i * slices + j)] = first;
        indices[6*slices + 6*(i * slices + j)+1] = second;
        indices[6*slices + 6*(i * slices + j)+2] = third;

        // Second triangle
        indices[6*slices + 6*(i * slices + j)+3] = third;
        indices[6*slices + 6*(i * slices + j)+4] = second;
        indices[6*slices + 6*(i * slices + j)+5] = fourth;
    }
}

}
