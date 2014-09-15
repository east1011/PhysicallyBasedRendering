#version 120

uniform vec3 uLight1Pos, uLight2Pos;

uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;

uniform vec3 uSunRayDir;
uniform vec3 uLightRGBColor;

//uniform sampler2D uTexUnit0;
varying vec2 vTexCoord;

varying  vec3 vNormal;
varying  vec3 vPosition;

//out vec4 fragColor;

void main() {


// 	 iproj=gl_ModelViewProjectionMatrix * intersection;
// iproj.z /= iproj.w;
// gl_FragDepth=(iproj.z+1.0)/2.0;

// float ndcDepth = ndcPos.z =
//    (2.0 * gl_FragCoord.z - gl_DepthRange.near - gl_DepthRange.far) /
//    (gl_DepthRange.far - gl_DepthRange.near);
//float clipDepth = ndcDepth / gl_FragCoord.w;
//gl_FragColor = vec4((clipDepth * 0.5) + 0.5); 

  //vec4 texColor0 = texture(uTexUnit0, vTexCoord);
  
  //vec3 toLight1 = normalize(uLight1Pos - vPosition);
  //vec3 toLight2 = normalize(uLight2Pos - vPosition);

  vec3 lightDir = normalize(uSunRayDir);

  vec3 normal = normalize(vNormal);

  float diffuseCoefficient = max( 0.0, dot(normal, lightDir) ) ;

 // vec3 diffuseColor = uLightRGBColor * diffuseCoefficient * vec3(uMaterialColor);
  vec3 diffuseColor =  diffuseCoefficient * vec3(uMaterialColor);
 
  vec3 toV = -normalize( vec3(vPosition) );
  vec3 h = normalize( toV + lightDir );

   
  float specularCoefficient  = pow( max(0.0, dot(h, normal) ), 30.0 );

  vec3 whiteLight = vec3( 0.6, 0.6, 0.6);

  vec3  specularColor = whiteLight * specularCoefficient;
     
  vec3  ambientColor = vec3(0.1, 0.1, 0.1);


  vec3 combinedColor = ambientColor + diffuseColor + specularColor;

  
  vec4 lightingColor  = vec4(combinedColor, 1.0);	

  //fragColor = lightingColor;
  gl_FragColor = vec4(diffuseColor,1);

}
