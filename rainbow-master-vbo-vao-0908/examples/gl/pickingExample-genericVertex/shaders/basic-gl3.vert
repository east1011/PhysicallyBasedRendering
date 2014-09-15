#version 130

uniform mat4 uProjectionMatrix;
uniform mat4 uModelMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

in vec3 aPosition;
in vec3 aNormal;
in vec2 aTexCoord;

out vec3 vNormal;
out vec3 vPosition;
out vec2 vTexCoord;

void main() {
  vNormal = vec3(uNormalMatrix * vec4(aNormal, 0.0));

  // send position (eye coordinates) to fragment shader
  vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);

 // vec4 mPosition = uModelMatrix * vec4(aPosition, 1.0);
  vPosition = vec3(tPosition);
  gl_Position = uProjectionMatrix * tPosition;

   vTexCoord = aTexCoord;

  // gl_Position is the screen position; In screen space, the coordinates (?1, ?1) and (1, 1) correspond respectively to the lower-left
  // and upper-right corners of the framebuffer; gl_Position's other two vector components are used in depth testing and perspective projection; 
  
}