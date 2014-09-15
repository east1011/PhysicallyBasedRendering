vec3 query (inout int level, const vec3 p) {

	// check if a cell exists at the current 'level' and position P
	vec3 k = S(level, p);	// generate a key for the cell
	ivec2 cell = texture3D (cellTex, H(k)).ra	//fetch cell
	int voxels = cell.y;	// save the address of the cell's voxel block

	if (computeId(k) == cell.x) {	// is this the right cell?
		
		while (level <= maxLevel ) {	// cell exists, try to descend the tree
			
			k = S( ++level, p );	// generate key for the child cell contatining P
			cell = texture3D(cellTex, H(k)).ra;	//fetch cell
			
			if (computeId(k) != cell.x) {	// is this the right cell?
				--level;	// no cell at this level, go back to the last valid level
				break;
			}
			voxels = cell.y;	//keep track of the smallest cell
		}
	}
	else {	// no cell at the current 'level' and position 'p'
		while (level > -) {	//ascend for octree untill we find a valid cell
			k = S (--level, p);		//generate key for the parent cell
			cell = texture3D( cellTex, H(k)).ra;	// fetch cell
			if (computeId(k) == cell.x) {	// is this the right cell?
				voxels = cell.y;	// the first cell we find upwards is the smallest
				break;				so we are done
			}
		}
	}	// map ID address of the cell's voxel block to an unnormalized 3D texture coord
	return mapIdto3D (voxels, voxelTexDim) + 0.5;
}	


float reconstructDist (const vec3 vb, const vec3 op) {
	return texture3D(voxelTex, (vb + cp) / voxelTexDim).a;
}

// Fragment Program

uniform vec3 eye;			// ray origin in world space
uniform float isovalue;		// value of the isosurface; usually()
uniform float tolerance;	// intersection tolerance for sphere-tracing
varying vec3 fragment-pos;	// from the vertex program, in world space

void main()
{
	// ray starts at the fragment's position, on the surface of the ADF's AABB
	vec3 pos = fragment_pos;	// current position along the ray
	vec3 dir = normalize (pos - eye);	// ray direction
	// compute the distances at which the ray enters leaves the ADF's AABB
	float rayLength, maxRayLength;
	
	intersectAABB( eye, dir, rayLength, maxRayLength );
	// traverse the distance field until we hit the isosurface or leave the AABB
	int level = 1;	//initial level for queries

	for ( int i = 0; i < MAX_ITERATIONS; ++i ) {
		vec3 vb = query (level, pos);	//locate the cell's voxel block
		vec3 cp = toCellSpace (pos, level);	//pos into cell space
		float d = reconstructDist (vb, cp) - isovalue;
		rayLength += d;	//step along the ray
		pos = eye + dir * rayLength;	// update the current position
		if (d < tolerance ) {	//did we hit the isosurface?
			vec3 normal = reconstructNormal (vb, cp) ;
			gl_FragColor = shade (pos, normal, -dir);
			return;
		}
		if (rayLength >= maxRayLength)	// did we leave the AABB?
			discard;
	}
	discard;	//we have exceeded MAX_ITERATIONS
}
