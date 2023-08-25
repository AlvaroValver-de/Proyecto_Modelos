#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

// texture samplers
uniform sampler2D texture1;
uniform sampler2D texture2;
uniform sampler2D texture3;

uniform float figura;

void main()
{
	if (figura == 0.0f){
        FragColor = texture(texture1, TexCoord);
    }
    if (figura == 1.0f) {
        FragColor = texture(texture2, TexCoord);
    }
    if (figura == 2.0f) {
        FragColor = texture(texture3, TexCoord);
    }
}