#version 450
#extension GL_ARB_separate_shader_objects : enable

// The attributes associated with each vertex.
// Their type and location must match the definition given in the
// corresponding Vertex Descriptor, and in turn, with the CPP data structure
layout(location = 1) in vec2 inUV;

// this defines the variable passed to the Fragment Shader
// the locations must match the one of its in variables
layout(location = 0) out vec2 fragUV;

// Here the Uniform buffers are defined. In this case, the Transform matrices (Set 1, binding 0)
// are used. Note that the definition must match the one used in the CPP code
layout(set = 0, binding = 0) uniform UniformBufferObject {
	mat4 mvpMat;
} ubo;

// Here the shader simply computes clipping coordinates, and passes to the Fragment Shader
// the position of the point in World Space, the transformed direction of the normal vector,
// and the untouched (but interpolated) UV coordinates
	// Clipping coordinates must be returned in global variable gl_Posision
	gl_Position = ubo.mvpMat * vec4(inPosition, 1.0);
	// Here the value of the out variables passed to the Fragment shader are computed
	fragUV = inUV;
}