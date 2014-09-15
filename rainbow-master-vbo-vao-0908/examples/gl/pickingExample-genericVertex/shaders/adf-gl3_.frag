#version 130
# 
// Fragment Program:
//uniform vec3 uEyePos; // This is simply (0,0,0) in eye space => see below

// ray origin in world space
uniform float uIsovalue;
// value of the isosurface; usually 0

//uniform float uTolerance;

// intersection tolerance for sphere-tracing
uniform vec3 uLight, uLight2, uColor;
uniform sampler3D uDistTex;

uniform sampler2D uScatTex;

// AABB bounding box
uniform vec3 uAABBmax;
uniform vec3 uAABBmin;
// AABB[0] =(x0, y0, z0), AABB[1] = (x1, y1, z1)
uniform vec3 uSampleLength;
uniform mat4 uEyeMatrix;
uniform mat4 uInvEyeMatrix;
uniform float uZDistToAABB;


// from the vertex program, in world space
out vec3 vPosition;
out vec4 fragColor;

vec3 uEyePos = vec3(0,0,0);
float uTolerance = 0.1;

vec3 AABB[2];

#define MAX_ITERATIONS 1000

struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
};

Ray makeRay( vec3 origin, vec3 direction) {
 vec3 inv_direction = vec3(1.0) / direction;

 Ray retRay;

 retRay.origin = origin;
 retRay.direction = direction;
 retRay.inv_direction= inv_direction;
 retRay.sign[0] = (inv_direction.x < 0.0) ? int(1) : int(0);
 retRay.sign[1] = (inv_direction.y < 0.0) ? int(1) : int(0);
 retRay.sign[2] = (inv_direction.z < 0.0) ? int(1) : int(0);

 return retRay;

}

void intersection_distance_no_if ( 
   in Ray ray,   in vec3 aabb[2],
   out float tmin, out float tmax ) 
 {
   float tymin, tymax, tzmin, tzmax;
   float tmp;

   //swap min.z and max.z
   //tmp = aabb[1].z;
   //aabb[1].z = aabb[0].z;
   //aabb[0].z = tmp;

   tmin = ( aabb[ ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tmax = ( aabb[1- ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tymin = ( aabb[ ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tymax = ( aabb[ 1 - ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tzmin = ( aabb[ ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
   tzmax = ( aabb[ 1- ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
   tmin = max( max( tmin, tymin), tzmin );
   tmax = min( min( tmax, tymax), tzmax );
 
 }
    
// diffuse shading. toV is not being used for now.
vec4 shade (const vec3 pos, const vec3 normal, const vec3 toV)
{
    vec3 tolight = normalize(uLight - pos);
    vec3 tolight2 = normalize(uLight2 - pos);
    float diffuse = max(0.0, dot(normal, tolight));
    diffuse += max(0.0, dot(normal, tolight2));
    vec3 intensity = uColor * diffuse;
    return vec4(intensity, 1.0);
}

vec3 posInAABB (vec3 pos) {

	//vec3 AABB[2];
	//AABB[0] = uAABBmin;
	//AABB[1] = uAABBmax;
		
	float AABBwidth = abs(AABB[1].x - AABB[0].x);
	float AABBheight = abs(AABB[1].y - AABB[0].y);
	float AABBdepth = abs(AABB[1].z - AABB[0].z);

	//pos = vec3(pos[0], pos[1], pos[2] + uZDistToAABB);

	//vec3 scaledPos = vec3( (pos.x + AABBwidth / 2) / AABBwidth,
	 //           (pos.y + AABBheight / 2) / AABBheight,
	//	    	(pos.z + AABBdepth / 2) / AABBdepth);	// scale to [0,1]

	// pos = AABB[0] + posWithinAABB, where all coordinates are eye coordinates
	vec3 posWithinAABB = pos - AABB[0];
	
	vec3 scaledPos = posWithinAABB / vec3(AABBwidth, AABBheight, AABBdepth);

	return scaledPos;
}

vec3 calculateNormal(const vec3 pos)
{
  vec3 normal;
  const float bias = 0.5;	
  
  vec3 posLeft, posRight, posUp, posDown, posFront, posRear;
  float dLeft, dRight, dUp, dDown, dFront, dRear;
  /*
  posLeft = posInAABB(vec3(pos.x-bias, pos.y, pos.z));
  posRight = posInAABB(vec3(pos.x+bias, pos.y, pos.z));
  posUp = posInAABB(vec3(pos.x, pos.y+bias, pos.z));
  posDown = posInAABB(vec3(pos.x, pos.y-bias, pos.z));
  posFront = posInAABB(vec3(pos.x, pos.y, pos.z+bias));
  posRear = posInAABB(vec3(pos.x, pos.y, pos.z-bias));
  
  dLeft = texture3D (uDistTex, posLeft).r;
  dRight = texture3D (uDistTex, posRight).r;
  dUp = texture3D (uDistTex, posUp).r;
  dDown = texture3D (uDistTex, posDown).r;
  dFront = texture3D (uDistTex, posFront).r;
  dRear = texture3D (uDistTex, posRear).r;
  */
  
  //float deltaX = uSampleLength.x;
  //float deltaY = uSampleLength.y;
  //float deltaZ = uSampleLength.z;

  float deltaX = 1.0/256.0;
  float deltaY = 1.0/256.0;
  float deltaZ = 1.0/256.0;
  //float deltaX = 0.1;
  //float deltaY = 0.1;
  //float deltaZ = 0.1;


  vec3 indexPos = posInAABB(pos);
  dLeft = texture3D (uDistTex, vec3(indexPos.x-deltaX, indexPos.y, indexPos.z)).r;
  dRight = texture3D (uDistTex, vec3(indexPos.x+deltaX, indexPos.y, indexPos.z)).r;
  dUp = texture3D (uDistTex, vec3(indexPos.x, indexPos.y+deltaY, indexPos.z)).r;
  dDown = texture3D (uDistTex, vec3(indexPos.x, indexPos.y-deltaY, indexPos.z)).r;
  dFront = texture3D (uDistTex, vec3(indexPos.x, indexPos.y, indexPos.z+deltaZ)).r;
  dRear = texture3D (uDistTex, vec3(indexPos.x, indexPos.y, indexPos.z-deltaZ)).r;
  
  normal = vec3(dRight - dLeft, dUp - dDown, dFront - dRear);
  
  // temp normal
  //normal = vec3 (1,0,0);
    
  return normal;
}


void main()
{  
    // ray starts at the fragment’s position, on the surface of the ADF’s AABB

	// current position along the ray
	vec3 pos = vec3(vPosition.x, vPosition.y, vPosition.z);
    // ray direction
	vec3 dir = normalize( pos - uEyePos );

	//vec3 eyePosGlobal = vec3(uEyeMatrix * vec4(uEyePos, 1));
	//Ray ray = makeRay( eyePosGlobal, dir);

	Ray ray = makeRay( uEyePos, dir);
	float rayLength, maxRayLength;

	/*
	vec3 AABB[2];
	AABB[0] = uAABBmin;
	AABB[1] = uAABBmax;
	*/

	// For Debugging: paint the scat data 

	vec3 indexPos = posInAABB(pos);

	fragColor =  vec4( texture( uScatTex, vec2(indexPos)  ).r * 1.0e+6,  0,0,1);
	return;

	//intersectAABB( eye, dir, rayLength, maxRayLength );  // original code
 
	// compute the distances at which the ray enters leaves the ADF's AABB

	// uAABBmin and uAABBmax are in the world coordinates, change them to eye coordinates
	// because all the other coordinates are in eye coordinate system

	AABB[0] = vec3( uInvEyeMatrix * vec4( uAABBmin,1) );
	AABB[1] = vec3( uInvEyeMatrix * vec4(uAABBmax,1) );
	
	intersection_distance_no_if( ray, AABB, rayLength, maxRayLength );		

	//gl_FragColor = vec4 ( ( (rayLength - uZDistToAABB + abs(AABB[1].z - AABB[0].z) / 2)), 0, 0, 1);
	//return;
		
	//intersection not working properly now.
	//rayLength = uZDistToAABB - abs(AABB[1].z - AABB[0].z)/2 + 60.0;
	//rayLength = uZDistToAABB;
	//maxRayLength = 500.0;

	// test the result of tmin and tmax
  
	if ( rayLength > maxRayLength ) { // not hit AABB => paint black for debug
	  fragColor  = vec4(0,1,0,1);
	  return;

	}



	// traverse the distance ﬁeld until we hit the isosurface or leave the AABB
	
	int level = 1;			// initial level for queries

	for( int i = 0; i < MAX_ITERATIONS; ++i ) {

		// original code
        //vec3 vb = query( level, pos );			// locate the cell’s voxel block
		//vec3 cp = toCellSpace( pos, level );		// pos into cell space
		//float d = reconstructDist( vb, cp ) - isovalue;

		//vec4 globalPos = vec4(pos, 1);

		vec3 indexPos = posInAABB(pos);

		// get distances to D 
		float d = texture (uDistTex, indexPos).r - uIsovalue;

		// step along the ray
		rayLength += d;
        // update the current position
		pos = uEyePos + dir * rayLength;

		if( d < uTolerance ) {	// ray hit the surface
			
			fragColor = vec4( rayLength/ maxRayLength, 0,0,1); // paint red with different intensity

			// calculate the normal of the position
			vec3 normal = calculateNormal(pos);
			
			// paint that point
           // gl_FragColor = shade( pos, normal, -dir );
			//gl_FragColor = vec4(d/maxRayLength * 10,0,0,1);
			return;
        }

		
	   if ( rayLength > maxRayLength ) { // reached the back face of  AABB without
	                                     // "hitting" the surface => paint black for debug
	     fragColor = vec4(0,0,0,1);
	     return;

	   }

	} // for 

	fragColor  = vec4 (0,1,0,1);	// paint green  for debug
	return;
	//discard; // we have exceeded MAX ITERATIONS
} // endmain


/*
// main for debugging

void main() {

	float a;

	vec3 AABB[2];
	AABB[0] = uAABBmin;
	AABB[1] = uAABBmax;

	float AABBwidth = AABB[1].x - AABB[0].x;
	float AABBheight = AABB[1].y - AABB[0].y;
	float AABBdepth = AABB[1].z - AABB[0].z;

	//vec4 globalPos = vec4(vPosition[0], vPosition[1], vPosition[2] + uZDistToAABB - AABBdepth, 1);

	//vec3 indexPos = vec3( (globalPos[0] + AABBwidth / 2) / AABBwidth,
	 //                     (globalPos[1] + AABBheight / 2) / AABBheight,
	//					  (globalPos[2] + AABBdepth / 2) / AABBdepth);	// scale to [0,1]

	//if (indexPos.z == AABBdepth/2 ) // this test should be (indexPos.z == 1.0), because indexPos is scaled
	                                // the vPosition is on the front face of AABB box.
    
	//gl_FragColor = texture3D (uDistTex, vec3(indexPos.x, indexPos.y, indexPos.z)) / 55.0;

	vec3 indexPos = posInAABB(vPosition);
	
	if (indexPos.z == 1.0 )
		gl_FragColor = vec4(1,0,0,1);

	else
		gl_FragColor = vec4(0,0,1,1);
		

			
	//if (indexPos.x < .2)
	//	gl_FragColor = vec4(1,0,0,1);
	//else if (indexPos.x < .4)
	//	gl_FragColor = vec4(0,1,0,1);
	//else if (indexPos.x < .6)
	//	gl_FragColor = vec4(0,0,1,1);
	//else if (indexPos.x < .8)
	//	gl_FragColor = vec4(1,0,1,1);
	//else if (indexPos.x < 1.0)
	//	gl_FragColor = vec4(1,1,0,1);
	//else
	//	gl_FragColor = vec4(0,1,1,1);
		

}
*/

/* not being used */
//vec3 query( inout int level, const vec3 p ) {		// check if a cell exists at the current ‘level’ and position ‘p’
	
//	vec3 k = S(level, p);							// generate a key for the cell
//	ivec2 cell = texture3D( cellTex, H(k) ).ra;		// fetch cell
//	int voxels = cell.y;							// save the address of the cell’s voxel block
//	if( computeId(k) == cell.x ) {					// is this the right cell?
//		while( level <= maxLevel ) {				// cell exists, try to descend the tree
//			k = S( ++level, p );					// generate key for the child cell containin } p
//			cell = texture3D( cellTex, H(k) ).ra	// fetch cell
//			if( computeId(k) != cell.x ) {			// is this the right cell?
//				--level;							// no cell at this level, go back to the last valid level
//				break;								// the last visited cell is the smallest one
//			}
//			voxels = cell.y;						// keep track of the smallest cell
//		}
//	} else {										// no cell at the current ‘level’ and position ‘p’
//		while( level > 0 ) {						// ascend the octree until we ﬁnd a valid cell
//			k = S( --level, p );					// generate key for the parent cell
//			cell = texture3D( cellTex, H(k) ).ra;	// fetch cell
//			if( computeId(k) == cell.x ) {			// is this the right cell?
//				voxels = cell.y;					// the ﬁrst cell we ﬁnd upwards is the smallest
//				break;								// so we are done
//			}
//		}
//	} // map 1D address of the cell’s voxel block to an unnormalized 3D texture coord
//	return map1Dto3D( voxels, voxelTexDim ) + 0.5;
//}

//float reconstructDist( const vec3 vb, const vec3 cp ) {
//    return texture3D(voxelTex, (vb + cp) / voxelTexDim).a;

//}

/*
float RGBAToFloat( vec4 rgba ) {
  return dot( rgba, vec4(1.0, 1/255.0, 1/65025.0, 1/16581375.0) );
}

*/