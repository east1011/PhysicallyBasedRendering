#version 130

uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

in vec3 aPosition;
in vec3 aNormal;
in vec2 aTexCoord;

out vec3 vNormal;
out vec3 vPosition;
out vec2 vTexCoord;

void main() {

  vNormal = vec3( uNormalMatrix * vec4(aNormal, 0.0) );

  // send position (eye coordinates) to fragment shader
  
  vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);
  //vec4 tPosition =  vec4(aPosition, 1.0);
  vPosition = vec3(tPosition);

  // compute the clip coordinates ( not window coordinates)
  // After clipping, the surviving coordinates are divided by the w component
  // to get normalized device coordinates in (-1, 1). A transformation will then be applied to move from NDC space 
  // to window coordinates, where the X and Y coordinates are normalized based on the viewport provided to OpenGL 
  // and the Z coordinate is normalized based on the depth range, which is
  // ultimately what gives you your (0, 1) range for depth (unless you use glDepthRange to set a different range).

  // The window coordinates are computed with the given parameters of the above 2 functions; 
//glViewport(x, y, w, h); 
//glDepthRange(n, f); n= depth value that vertices at the near clipping  plane are mapped to; f= depth value that vertices at the far 
 // clipping plane are mapped to.


  gl_Position = uProjectionMatrix * tPosition;
   
  vTexCoord = aTexCoord;
 

}