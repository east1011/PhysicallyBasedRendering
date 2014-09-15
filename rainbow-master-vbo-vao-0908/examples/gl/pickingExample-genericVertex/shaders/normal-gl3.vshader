#version 130

uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

in vec3 aPosition;
in vec3 aNormal;   // N
in vec3 aTangent;  // S
in vec3 aBinormal;  // B
in vec2 aTexCoord;

out vec2 vTexCoord;
out mat3 vTBNMat;  // normal matrix * tangent frame matrix
out vec3 vPosition; // position in eye space
out vec3 vNormal;

void main() {

  vNormal = vec3( uNormalMatrix * vec4(aNormal, 0.0) );
  vTexCoord = aTexCoord;

  vTBNMat = mat3(uNormalMatrix) * mat3(aTangent, aBinormal, aNormal);
  vec4 posE = uModelViewMatrix * vec4(aPosition, 1.0);
  vPosition  = posE.xyz;
  gl_Position = uProjectionMatrix * posE;
}