#version 450

layout(location = 0) in float real;
layout(location = 0) out vec4 outColor;
layout(set = 0, binding = 0) uniform sampler2D tex;
void main() {
	float m_real = 0.0f, m_img = 0.0f, temp;