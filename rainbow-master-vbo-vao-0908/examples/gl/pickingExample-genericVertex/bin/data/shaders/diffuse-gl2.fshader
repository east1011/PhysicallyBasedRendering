uniform vec3 uLight1Pos, uLight2Pos;

uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;

varying vec3 vNormal;
varying vec3 vPosition;


void main() {
  vec3 toLight1 = normalize(uLight1Pos - vPosition);
  vec3 toLight2 = normalize(uLight2Pos - vPosition);
  vec3 normal = normalize(vNormal);

  float diffuseCoefficient = max( 0.0, dot(normal, toLight1)) ;

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

  
  gl_FragColor = vec4(combinedColor.x, combinedColor.y, combinedColor.z, 1.0);	
  
}
