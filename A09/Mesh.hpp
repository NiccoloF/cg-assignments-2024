

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

void MakeSquare(float size, std::vector<std::array<float,6>> &vertices, std::vector<uint32_t> &indices) {
// Creates a square, on the xz-plane, aligned with the axis, and centered in the origin.
// The length of the four sides is in parameter >size<.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a std::array<float,6> element.
// In particular, the first three elements (index 0,1,2) encode the position of the vertex (x in index 0,
// y in index 1, and z in index 2). The second set of three elements (index 3,4,5) encode the direction
// of the normal vector for the considerd vertex (N.x in index 3, N.y in index 4, and N.z in index 5).
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: this procedure has already been implemented. You can keep it as is
	vertices = {
				   {-size/2.0f,0.0f,-size/2.0f,0.0f,1.0f,0.0f},
				   {-size/2.0f,0.0f, size/2.0f,0.0f,1.0f,0.0f},
				   { size/2.0f,0.0f,-size/2.0f,0.0f,1.0f,0.0f},
				   { size/2.0f,0.0f, size/2.0f,0.0f,1.0f,0.0f}};
	indices = {0, 1, 2,    1, 3, 2};
}

void MakeCube(float size, std::vector<std::array<float,6>> &vertices, std::vector<uint32_t> &indices) {
// Creates a cube, with the faces perpendicular to the axis, and centered in the origin.
// The length of one edge of the cube is >size<.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a std::array<float,6> element.
// In particular, the first three elements (index 0,1,2) encode the position of the vertex (x in index 0,
// y in index 1, and z in index 2). The second set of three elements (index 3,4,5) encode the direction
// of the normal vector for the considerd vertex (N.x in index 3, N.y in index 4, and N.z in index 5).
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a square. You can use it as a side of the cube (please remember
// to change the value of the y component, otherwise the result will be wrong
    // Vertices of a cube
    vertices.clear();
    indices.clear();

    vertices = {
            // Bottom face
            {-size / 2.0f, -size / 2.0f, -size / 2.0f, 0.0f, -1.0f, 0.0f}, // Bottom left
            {-size / 2.0f, -size / 2.0f,  size / 2.0f, 0.0f, -1.0f, 0.0f}, // Top left
            { size / 2.0f, -size / 2.0f,  size / 2.0f, 0.0f, -1.0f, 0.0f}, // Top right
            { size / 2.0f, -size / 2.0f, -size / 2.0f, 0.0f, -1.0f, 0.0f}, // Bottom right

            // Top face
            {-size / 2.0f, size / 2.0f, -size / 2.0f, 0.0f, 1.0f, 0.0f}, // Bottom left
            {-size / 2.0f, size / 2.0f,  size / 2.0f, 0.0f, 1.0f, 0.0f}, // Top left
            { size / 2.0f, size / 2.0f,  size / 2.0f, 0.0f, 1.0f, 0.0f}, // Top right
            { size / 2.0f, size / 2.0f, -size / 2.0f, 0.0f, 1.0f, 0.0f}, // Bottom right

            // Front face
            {-size / 2.0f, -size / 2.0f,  size / 2.0f, 0.0f, 0.0f, 1.0f}, // Bottom left
            {-size / 2.0f,  size / 2.0f,  size / 2.0f, 0.0f, 0.0f, 1.0f}, // Top left
            { size / 2.0f,  size / 2.0f,  size / 2.0f, 0.0f, 0.0f, 1.0f}, // Top right
            { size / 2.0f, -size / 2.0f,  size / 2.0f, 0.0f, 0.0f, 1.0f}, // Bottom right

            // Back face
            {-size / 2.0f, -size / 2.0f, -size / 2.0f, 0.0f, 0.0f, -1.0f}, // Bottom left
            {-size / 2.0f,  size / 2.0f, -size / 2.0f, 0.0f, 0.0f, -1.0f}, // Top left
            { size / 2.0f,  size / 2.0f, -size / 2.0f, 0.0f, 0.0f, -1.0f}, // Top right
            { size / 2.0f, -size / 2.0f, -size / 2.0f, 0.0f, 0.0f, -1.0f}, // Bottom right

            // Left face
            {-size / 2.0f, -size / 2.0f, -size / 2.0f, -1.0f, 0.0f, 0.0f}, // Bottom back
            {-size / 2.0f,  size / 2.0f, -size / 2.0f, -1.0f, 0.0f, 0.0f}, // Top back
            {-size / 2.0f,  size / 2.0f,  size / 2.0f, -1.0f, 0.0f, 0.0f}, // Top front
            {-size / 2.0f, -size / 2.0f,  size / 2.0f, -1.0f, 0.0f, 0.0f}, // Bottom front

            // Right face
            { size / 2.0f, -size / 2.0f, -size / 2.0f, 1.0f, 0.0f, 0.0f}, // Bottom back
            { size / 2.0f,  size / 2.0f, -size / 2.0f, 1.0f, 0.0f, 0.0f}, // Top back
            { size / 2.0f,  size / 2.0f,  size / 2.0f, 1.0f, 0.0f, 0.0f}, // Top front
            { size / 2.0f, -size / 2.0f,  size / 2.0f, 1.0f, 0.0f, 0.0f}  // Bottom front
    };

    // Indices of the cube triangles
    indices = {
            // Bottom face
            1, 0, 2,  2, 0, 3,
            // Top face
            4, 5, 6,  4, 6, 7,
            // Front face
            9, 8, 10,  10, 8, 11,
            // Back face
            12, 13, 14,  12, 14, 15,
            // Left face
            17, 16, 18,  18, 16, 19,
            // Right face
            20, 21, 22,  20, 22, 23
    };

}

void MakeCylinder(float radius, float height, int slices, std::vector<std::array<float,6>> &vertices, std::vector<uint32_t> &indices) {
// Creates a cylinder, approximated by a prism with a base composed by a regular polygon with >slices< sides.
// The radius of the cylinder is >radius<, and it height is >height<. The cylinder has its centere
// in the origin, and spans half above and half below the plane that passes thorugh the origin.
// The top and bottom caps are are aligned with xz-plane and perpendicular to the y-axis.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a std::array<float,6> element.
// In particular, the first three elements (index 0,1,2) encode the position of the vertex (x in index 0,
// y in index 1, and z in index 2). The second set of three elements (index 3,4,5) encode the direction
// of the normal vector for the considerd vertex (N.x in index 3, N.y in index 4, and N.z in index 5).
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a rectangle. You have to change it, or you will obtain a wrong result
// You should use a for loop, and you should start from the procedure to create a circle seen during the lesson
    // Clear existing data
    vertices.clear();
    indices.clear();
    vertices.resize(4 * slices + 2); // Aggiungiamo 2 vertici per i cerchi superiori e inferiori
    indices.resize(12*slices); // Ogni cerchio richiede slices * 3 indici

    // Aggiungi vertici per i cerchi superiori e inferiori
    vertices[4 * slices] = {0.0f, height/2, 0.0f,0.0f,1.0f,0.0f}; // Vertice centrale del cerchio superiore
    vertices[4 * slices + 1] = {0.0f, -height/2, 0.0f,0.0f,-1.0f,0.0f}; // Vertice centrale del cerchio inferiore

    for (int i = 0; i < slices; ++i) {
        // Create vertices for the side faces
        float ang = 2 * M_PI * (float)i / (float)slices;
        float nx = radius * cos(ang);
        float nz = radius * sin(ang);
        vertices[i] = {nx, height/2, nz,0.0f,1.0f,0.0f}; // Vertici del cerchio superiore y-up
        vertices[i + slices] = {nx, -height/2, nz,0.0f,-1.0f,0.0f}; // Vertici del cerchio inferiore-normale y-up
        vertices[2*slices + i] = {nx, height/2, nz,nx,0.0f,nz}; // Vertici del cerchio superiore
        vertices[3*slices + i] = {nx, -height/2, nz,nx,0.0f,nz};// Vertici del cerchio inferiore


        // top circle
        indices[3*i] = (i+1)%slices;
        indices[3*i+1] = i;
        indices[3*i+2] = 4*slices;

        //bottom circle
        indices[3*(i + slices)] = i + slices;
        indices[3*(i+slices)+1] = (i+1)%slices + slices;
        indices[3*(i+slices)+2] = 4*slices + 1;

        // Vertici del lato superiore del cilindro
        indices[6 * slices + 3 * i] = i + 2*slices;
        indices[6 * slices + 3 * i+1] = (i + 1) % slices + 2*slices;
        indices[6 * slices + 3 * i+2] = i + 3*slices;

        // Vertici del lato inferiore del cilindro
        indices[9 * slices + 3 * i] =  (i+1)%slices + 3*slices;
        indices[9 * slices + 3 * i+1] =  i + 3*slices;
        indices[9 * slices + 3 * i+2] = (i + 1) % slices + 2*slices;

    }

}

void MakeCone(float radius, float height, int slices, std::vector<std::array<float,6>> &vertices, std::vector<uint32_t> &indices) {
// Creates a cone, approximated by a pyramid with a base composed by a regular polygon with >slices< sides.
// The radius of the cone is >radius<, and it height is >height<. The cone has its centere
// in the origin, and spans half above and half below the plane that passes thorugh the origin.
// The bottom cap is aligned with xz-plane and perpendicular to the y-axis.
// The procedure should fill the array contained in the >vertices< parameter, with the positions of
// the vertices of the primitive, expressed with their local coordinates in a std::array<float,6> element.
// In particular, the first three elements (index 0,1,2) encode the position of the vertex (x in index 0,
// y in index 1, and z in index 2). The second set of three elements (index 3,4,5) encode the direction
// of the normal vector for the considerd vertex (N.x in index 3, N.y in index 4, and N.z in index 5).
// Indices should be returned in the >indices< array, starting from 0, and up to vertices.size()-1.
// The primitive is encoded as an indexed triangle list, so the size of the >indices< array, should
// be a multiple of 3: each group of three indices, defines a different triangle.
//
// HINT: the procedure below creates a triangle. You have to change it, or you will obtain a wrong result
// You should use a for loop, and you should start from the procedure to create a circle seen during the lesson
    vertices.clear();
    indices.clear();
    vertices.resize(2*slices+2);
    indices.resize(6*slices);

    vertices[2*slices]= {0.0f,-height/2.0f,0.0f,0.0f,-1.0f,0.0f}; // bottom
    vertices[2*slices+1]= {0.0f,height/2.0f,0.0f,0.0f,1.0f,0.0f}; // top
    for(int i = 0; i < slices; i++) {
        float ang = 2*M_PI * (float)i / (float)slices;
        float nx = radius * cos(ang);
        float nz = radius * sin(ang);
        vertices[i] = {nx,-height/2.0f,nz,0.0f,-1.0f,0.0f};
        glm::vec3 approx = glm::normalize(glm::vec3(height*glm::cos(ang),radius,height*glm::sin(ang)));
        vertices[i+slices] = {nx,-height/2.0f,nz,approx[0],approx[1],approx[2]};

        // bottom circle
        indices[3*i] = 2*slices;
        indices[3*i+1] = i;
        indices[3*i+2] = (i+1) % slices;

        // side faces
        indices[3*(i+slices)  ] = i+slices;
        indices[3*(i+slices)+1] = 2*slices+1;
        indices[3*(i+slices)+2] = (i+1) % slices + slices;
    }
}

void MakeSphere(float radius, int rings, int slices, std::vector<std::array<float,6>> &vertices, std::vector<uint32_t> &indices)
{
    vertices.clear();
    indices.clear();

    // Generate vertices
    for (int i = 0; i <= rings; ++i) {
        float phi = float(i) * float(M_PI) / float(rings);
        for (int j = 0; j <= slices; ++j) {
            float theta = float(j) * 2.0f * float(M_PI) / float(slices);
            float x = radius * std::sin(phi) * std::cos(theta);
            float y = radius * std::cos(phi);
            float z = radius * std::sin(phi) * std::sin(theta);
            float nx = x / radius;
            float ny = y / radius;
            float nz = z / radius;
            vertices.push_back({x, y, z, nx, ny, nz});
        }
    }

    // Generate indices
    for (int i = 0; i < rings; ++i) {
        for (int j = 0; j < slices; ++j) {
            int first = (i * (slices + 1)) + j;
            int second = first + slices + 1;
            indices.push_back(second);
            indices.push_back(first);
            indices.push_back(first + 1);

            indices.push_back(second + 1);
            indices.push_back(second);
            indices.push_back(first + 1);
        }
    }
}
