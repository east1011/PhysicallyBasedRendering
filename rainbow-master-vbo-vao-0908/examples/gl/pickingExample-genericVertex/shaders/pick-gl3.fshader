#version 130

uniform vec4 uMaterialColor;

out vec4 fragColor;

void main() {
fragColor = uMaterialColor;
//fragColor = vec4(0, 0, 200/255.0, 1.0);


}
