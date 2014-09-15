#version 130

uniform vec3 uLight1Pos, uLight2Pos;

uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;

uniform float uVertexScale;
uniform sampler2D uTexUnit0, uTexUnit1;
in vec2 vTexCoord0, vTexCoord1;


// uniform lowp float textureFlag;


//gl_FragColor = textureFlag * texture2D(texture, texcoordVarying) * colorVarying + 
//               (1.0 - textureFlag) * colorVarying;Or even:

//gl_FragColor = mix(
//                     colorVarying,
//                     texture2D(texture, texcoordVarying) * colorVarying,
 //                    textureFlag
// OR use different shaders for different materials


in vec3 vNormal;
in vec3 vPosition;

out vec4 fragColor;

void main() {

  vec4 texColor0 = texture(uTexUnit0, vTexCoord0);
  vec4 texColor1 = texture(uTexUnit1, vTexCoord1);
  float lerper = clamp(.5 *uVertexScale, 0., 1.);
  float lerper2 = clamp(.5 * uVertexScale + 1.0, 0.0, 1.0);



  vec3 toLight1 = normalize(uLight1Pos - vPosition);
  vec3 toLight2 = normalize(uLight2Pos - vPosition);

  vec3 normal = normalize(vNormal);

  float diffuseCoefficient = max( 0.0, dot(normal, toLight1) ) ;

  float diffuseIntensityR = uLight1Color.r * diffuseCoefficient * uMaterialColor.r;
  float diffuseIntensityG = uLight1Color.g * diffuseCoefficient * uMaterialColor.g;
  float diffuseIntensityB = uLight1Color.b * diffuseCoefficient * uMaterialColor.b;
  
  vec3  diffuseColor = vec3( diffuseIntensityR, diffuseIntensityG, diffuseIntensityB);

  vec3 toV = -normalize( vec3(vPosition) );
  vec3 h = normalize( toV + toLight1 );

   
  float specularCoefficient  = pow( max(0.0, dot(h, normal) ), 30.0 );

  vec3 whiteLight = vec3( 0.6, 0.6, 0.6);

  vec3  specularColor = whiteLight * specularCoefficient;
     
  vec3  ambientColor = vec3(0.1, 0.1, 0.1);


  vec3 combinedColor = ambientColor + diffuseColor + specularColor;

  
  vec4 lightingColor  = vec4(combinedColor.x, combinedColor.y, combinedColor.z, 1.0);	
  fragColor = ((lerper)*texColor1 + (1.0-lerper)*texColor0 ) * lerper2 + lightingColor * (1.0-lerper2);

}
