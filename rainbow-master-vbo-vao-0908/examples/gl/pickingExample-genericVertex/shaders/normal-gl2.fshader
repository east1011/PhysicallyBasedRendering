#version 120

uniform sampler2D uTexColor;
uniform sampler2D uTexNormal;

// lights in eye space
uniform vec3 uLight1Pos;
//uniform vec3 uLight2Pos;

varying  vec2 vTexCoord;
varying  mat3 vNTMat;
varying vec3 vPosition; // varying geometry position in eye space
varying  vec3  vNormal;

//out vec4 fragColor;

void main() {
 // fragColor = vec4(0,1,0,1);
 // return;

  // TODO: replace the following line with loading of normal from uTexNormal
  //       transforming to eye space, and normalizing

  // for debugging
  //fragColor = vec4(1,0,0,1);
  //return;

  vec3 normal = vec3(0, 0, 1);

  vec3 viewDir = normalize(-vPosition);
  vec3 lightDir = normalize(uLight1Pos - vPosition);
  //vec3 lightDir2 = normalize(uLight2Pos - vPosition);

  float nDotL = dot(normal, lightDir);
  float diffuse = max(nDotL, 0.0);

  vec3 reflectionDir = normalize( 2.0 * normal * nDotL - lightDir);
  float rDotV = max(0.0, dot(reflectionDir, viewDir));
  float specular = pow(rDotV, 32.0);
  

  //nDotL = dot(normal, lightDir2);
  //diffuse += max(nDotL, 0.0);

 
  //reflection = normalize( 2.0 * normal * nDotL - lightDir2);
  //rDotV = max(0.0, dot(reflection, viewDir));
  
  //specular += pow(rDotV, 32.0);


  //vec3 color = texture(uTexColor, vTexCoord).xyz * diffuse + specular * vec3(0.6, 0.6, 0.6);
  
  vec3 color = texture2D(uTexColor, vTexCoord).xyz * diffuse; 

  gl_FragColor = vec4(color, 1);
}
