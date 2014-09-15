#version 120

precision highp float;

uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

attribute  vec3 aPosition;
attribute vec3 aNormal;
attribute  vec2 aTexCoord;

varying  vec3 vNormal;
varying  vec3 vPosition;
varying  vec2 vTexCoord;

void main() {
  //vNormal = vec3( uNormalMatrix * vec4( aNormal, 0.0)  );
  //vNormal = aNormal;

  // send position (eye coordinates) to fragment shader
  //vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);

  //vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);

 // vPosition = vec3(tPosition);

  //gl_Position = uProjectionMatrix * tPosition;
  
  gl_Position = vec4( aPosition, 1); // no projection, aPosition itself the  position in viewspace
                                     // because the w coordinate is 1

  vTexCoord = aTexCoord;
 
 /*
 All executions of a well-formed vertex shader executable must write a
value into this variable. It can be written at any time during shader execution. It may also be read back
by a vertex shader after being written. This value will be used by primitive assembly, clipping, culling,
and other fixed functionality operations that operate on primitives after vertex processing has occurred.
*/

}