#version 130

uniform mat4 uEyeMatrix;
uniform sampler2D uTexColor;
uniform sampler2D uTexNormal;

uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;

// lights in eye space
uniform vec3 uLight1Pos;
uniform vec3 uLight2Pos;

// eyePos and aabb box for rainbow
//uniform vec3 uEye; // in eye space, so it should be (0,0,0): see below

uniform float uTolerance;
uniform float uRadius;
uniform vec3  uSunRay;

uniform vec3[2] uAABB;

in vec2 vTexCoord;
in mat3 vNTMat;

in vec3 vPosInEye; // varying geometry position in eye space

out vec4 fragColor;

vec3 uEye = vec3(0,0,0);


#define MAX_ITERATIONS 100

struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
}

Ray makeRay( vec3 origin, vec3 direction) {
 vec3 inv_direction = vec3(1.0) / direction;

 return Ray( 
      origin, direction, 
	  inv_direction,
	  { (inv_direction.x < 0): 1 ? 0,
	    (inv_direction.y < 0): 1 ? 0,
		(inv_direction.z < 0): 1 ? 0
	  }
 };

}

void intersection_distance_no_if ( 
   in Ray ray,   in vec3 aabb[2],
   out float tmin, out float tmax )
   
 {
   float tymin, tymax, tzmin, tzmax;
   tmin = ( aabb[ ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tmax = ( aabb[1- ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tymin = ( aabb[ ray.sign[0]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tymax = ( aabb[ 1 - ray.sign[0]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tzmin = ( aabb[ ray.sign[0] ].z - ray.origin.z ) * ray.inv_direction.z;
   tzmax = ( aabb[ 1- ray.sign[0] ].z - ray.origin.z ) * ray.inv_direction.z;
   tmin = max( max( tmin, tymin), tzmin );
   tmax = min( min( tmax, tymax), tzmax );
 
 }

 void shadeSurface( in vec3 vEyePos, in vec3 normal, out vec3 diffuse, out vec3 specular ) {
 
  vec3 viewDir = normalize(-vEyePos);
  vec3 lightDir = normalize(uLight1Pos - vEyePos);
  vec3 lightDir2 = normalize(uLight2Pos - vEyePos);

  float nDotL = dot(normal, lightDir);

  vec3 reflection = normalize( 2.0 * normal * nDotL - lightDir);
  float rDotV = max(0.0, dot(reflection, viewDir));
  float specular = pow(rDotV, 32.0);
  float diffuse = max(nDotL, 0.0);

  nDotL = dot(normal, lightDir2);
  reflection = normalize( 2.0 * normal * nDotL - lightDir2);
  rDotV = max(0.0, dot(reflection, viewDir));
  specular += pow(rDotV, 32.0);
  diffuse += max(nDotL, 0.0);

  }

vec3 computeRainbowColor( in float uRadius, in vec3 uEye, in vec3 dir,
                            in vec3 uSunRay, out float tmin, out float tmax) {
							
							
}

							
							
vec3 computeRainbowColorForCurrentFragment(  ) {

    vec3 pos = vPosInEye;
    // current position along the ray
	vec3 dir = normalize( pos - uEye );
    // ray direction
	
	Ray ray = makeRay( uEye, dir);

	// compute the distances at which the ray enters and leaves the ADF’s AABB
	
	float t, tmin, tmax;
	    
	intersection_distance_no_if( ray, AABB, tmin, tmax );		
	
	t = tmin;

	// traverse the distance ﬁeld until we hit the surface or leave the AABB
	
	for( int i = 0; i < MAX_ITERATIONS; ++i ) {
        //vec3 vb = query( level, pos );						// locate the cell’s voxel block
		//vec3 cp = toCellSpace( pos, level );				// pos into cell space
		//float d = reconstructDist( vb, cp ) - isovalue;

		//float d = getDist(pos) - isovalue;

		
		// Convert eyespace pos into the global position in AABB coordinate system
		//  because the AABB axes are parallel to the world system;
		//  It is needed because AABB texture is relative the world system.
		//  Also the global position should be normalized relative the AABB dimension 
		//  in order to be used as an index to the 3D texture

		vec3 globalPos = uEyeMatrix * pos;
		vec3 indexPos = vec3( globalPos[0] / (uAABB[1].x - uAABB[0].x),
		                      globalPos[1] / (uAABB[1].y - uAABB[0].x),
							  globalPos[2] / (uAABB[1].z - uAABB[0].z );

		float d = texture3D (uDistTex, indexPos).r;

        ///////////////////////////////////////////////

		t += d;
        // step along the ray
		pos = uEye + dir * rayLength;
        // update the current position
		
		if( d < uTolerance ) {
            // did we hit the isosurface?
			
			// compute the rainbow color along the direction from the hit position to the
			// front face of the AABB box.
			 
			vec3 volumeColor = computeRainbowColor( uRadius, uEye, dir, uSunRay, tmin, tmax);
            return volumeColor;
        }

		
		if( t >= tmax  )						// did we leave the AABB?
			discard;
        // ray does not intersect the isosurface
	} // endfor
	discard;
    // we have exceeded MAX ITERATIONS
}


}
  // it is assumed that the light color is (1,1,1)
void main() {
  // TODO: replace the following line with loading of normal from uTexNormal
  //       transforming to eye space, and normalizing
  vec3 normal = vec3(0, 0, 1);
  vec3 diffuse, specular;

  shadeSurface( vEyePos, normal, diffuse, specular);

  vec3 surfaceColor = texture(uTexColor, vTexCoord).xyz * diffuse + specular * vec3(0.6, 0.6, 0.6);

  vec3 volumeColor = computeRainbowColorForCurrentFragment();

  fragColor = vec4(surfaceColor, 1) + vec4(volumeColor,1);
}

