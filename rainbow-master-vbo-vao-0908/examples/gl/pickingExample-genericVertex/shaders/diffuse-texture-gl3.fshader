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

//Remember how, in vertex shaders, you have inputs? And these inputs represent vertex attribute indices,
// the numbers you pass to glVertexAttribPointer and glEnableVertexAttribArray and so forth?
// You set up which input pulls from which attribute. In GLSL 3.30, you use this syntax:
//layout(location = 2) in color;


//This sets the color vertex shader input to come from attribute location 2.
// Before 3.30 (or without ARB_explicit_attrib_location), you would have to either set this up explicitly 
//with glBindAttrbLocation before linking or query the program for the attribute index with glGetAttribLocation. If you don't explicitly provide an attribute location, GLSL will assign a location arbitrarily (ie: in an implementation-defined manner).

//Setting it in the shader is almost always the better option.

//In any case, fragment shader outputs work almost exactly the same way. 
//Fragment shaders can write to multiple buffers. Therefore, you need to indicate
// which output goes to which buffer.

//This process begins with the fragment output location value. It's set very similarly to vertex shader 
//input locations:
//layout(location = 1) out secColor;


//There are also the API functions glBindFragDataLocation and glGetFragDataLocation, 
//which are analogous to glBindAttribLocation and glGetAttribLocation.

//But there's one important difference. When you don't explicitly specify an attribute location, 
//GL assigns them arbitrarily. When you don't explicitly specify a fragment shader output,
// GL 3.3 and above (that's a change from prior versions) will assign them all to location 0.

//And this makes sense. The most common case is rendering a single output to a single buffer. 
//Using the default draw buffer and the default framebuffer,
// location 0 is piped directly to the back buffer. 
//If you use glDrawBuffer to change what you render to, then you're changing what location 0 renders to.
// So everything just works.


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
