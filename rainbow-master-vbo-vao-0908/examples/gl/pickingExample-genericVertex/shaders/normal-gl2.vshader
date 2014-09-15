#version 120
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

attribute vec3 aPosition;
attribute vec3 aNormal;
attribute vec3 aTangent;
attribute vec3 aBinormal;
attribute vec2 aTexCoord;

varying vec2 vTexCoord;
//varying mat3 vNTMat;  // normal matrix * tangent frame matrix
varying mat3 vTBNMat;  // normal matrix * tangent frame matrix
varying vec3 vPosition; // position in eye space
varying vec3 vNormal;

void main() {

  vNormal = vec3( uNormalMatrix * vec4(aNormal, 0.0) );
  vTexCoord = aTexCoord;
  vTBNMat = mat3(uNormalMatrix) * mat3(aTangent, aBinormal, aNormal);
  //vNTMat = mat3(uNormalMatrix[0].xyz, uNormalMatrix[1].xyz, uNormalMatrix[2].xyz) * mat3(aTangent, aBinormal, aNormal);
  vec4 posE = uModelViewMatrix * vec4(aPosition, 1.0);
  vPosition = posE.xyz;
  gl_Position = uProjectionMatrix * posE;
}