#version 450

layout(location = 0) in vec2 fragUV;

layout(location = 0) out vec4 outColor;
	vec3 Emit = texture(tex, fragUV).rgb;
	outColor = vec4(Emit, 1.0f);