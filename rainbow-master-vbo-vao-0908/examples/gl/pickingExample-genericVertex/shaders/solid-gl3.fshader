#version 130

uniform vec4 uMaterialColor;
in vec3 vPosition;

out vec4 fragColor;
//out vec4 spect1;

void main() {
  fragColor = uMaterialColor;
 // spect1 = vec4( vPosition, 1.0);
      
}
