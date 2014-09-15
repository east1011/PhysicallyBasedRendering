
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

inline float RGBAToFloat( vec4 rgba ) {
  return dot( rgba, vec4(1.0, 1/255.0, 1/65025.0, 1/16581375.0) );
}

float getDist (vec3 pos) {

	// TODO: calculate texIndex from pos

	vec3 distInRGBA = texture3D(uDistTex, texIndex);
	return RGBAToFloat(distInRGBA);		// return distance
  
}

// diffuse shading. toV is not being used for now.
vec4 shade ( const vec3 pos, const vec3 normal, const vec3 toV)
{
  vec3 tolight = normalize(uLight - pos);
  vec3 tolight2 = normalize(uLight2 - pos);

  float diffuse = max(0.0, dot(normal, tolight));
  diffuse += max(0.0, dot(normal, tolight2));
  vec3 intensity = uColor * diffuse;

  return vec4(intensity, 1.0);
}



// Fragment Program:

//uniform vec3 eye;						// ray origin in world space
uniform float isovalue;					// value of the isosurface; usually 0
uniform float tolerance;				// intersection tolerance for sphere-tracing
uniform vec3 uLight, uLight2, uColor;
uniform sampler3D uDistTex;

varying vec3 vPosition;					// from the vertex program, in world space

/*
void main()
{										
// ray starts at the fragment’s position, on the surface of the ADF’s AABB

	vec3 pos = vPosition;				// current position along the ray
	vec3 dir = normalize( pos - eye );	// ray direction
	
	// compute the distances at which the ray enters and leaves the ADF’s AABB
	//float rayLength, maxRayLength;
	float rayLength = 0;
	//intersectAABB( eye, dir, rayLength, maxRayLength );		
	
	// traverse the distance ﬁeld until we hit the isosurface or leave the AABB
	int level = 1;											// initial level for queries
	for( int i = 0; i < MAX_ITERATIONS; ++i ) {
		//vec3 vb = query( level, pos );						// locate the cell’s voxel block
		//vec3 cp = toCellSpace( pos, level );				// pos into cell space
		//float d = reconstructDist( vb, cp ) - isovalue;

		//float d = getDist(pos) - isovalue;

		///////////////////////////////////////////////
		// TODO: get distances to D
		// float d = getDistance (cp);
		///////////////////////////////////////////////

		rayLength += d;										// step along the ray
		pos = eye + dir * rayLength;						// update the current position
		
		if( d < tolerance ) {								// did we hit the isosurface?
			
			// TODO: calculate normals! how do I calculate normals??
			// vec3 normal = reconstructNormal( vb, cp );

			gl_FragColor = shade( pos, normal, -dir );
			return;
		}
		
		if( rayLength >= maxRayLength )						// did we leave the AABB?
			discard;										// ray does not intersect the isosurface
	} // endfor
	discard;											// we have exceeded MAX ITERATIONS
}
*/

void main() {

	gl_FragColor = texture3D(uDistTex, vPosition);

}