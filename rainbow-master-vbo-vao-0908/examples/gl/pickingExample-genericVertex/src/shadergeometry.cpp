#include <stdexcept>
#include <string>
#include <cstddef>

#include "shadergeometry.h" // which includes "glsupport.h"

using namespace std;
using namespace tr1;

// extern QOpenGLFunctions_3_0 *m_funcs;  declared within geometry.h


const VertexFormat VertexPN::FORMAT = VertexFormat( sizeof(VertexPN) )
                                      .put("aPosition", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, p))
                                      .put("aNormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, n));

// FORMAT is a static member of VertexPNX class

const VertexFormat VertexPNX::FORMAT = VertexFormat(sizeof(VertexPNX))
                                       .put("aPosition", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, p))
                                       .put("aNormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, n))
                                       .put("aTexCoord", 2, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, x));

const VertexFormat VertexPNTBX::FORMAT = VertexFormat( sizeof(VertexPNTBX) )
                                         .put("aPosition", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, p))
                                         .put("aNormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, n))
                                         .put("aTangent", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNTBX, t))
                                         .put("aBinormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPNTBX, b))
                                         .put("aTexCoord", 2, GL_FLOAT, GL_FALSE, offsetof(VertexPNX, x));


BufferObjectGeometry::BufferObjectGeometry(const string geoName)
  : Geometry(geoName), wiringChanged_(true),
  primitiveType_(GL_TRIANGLES)
{}

BufferObjectGeometry& BufferObjectGeometry::wireAttribute (
  const string& targetAttribName,
  shared_ptr<FormattedVbo> source,
  const string& sourceAttribName) {

  wiringChanged_ = true;
  wiring_[targetAttribName] = make_pair(source, sourceAttribName); // sourceAttributeName ="aPosition", "aNormal". These are defined 
                                                                   // in Vertex::FORMAT
                                                                   // targetAttributeName refers to attribute names used in shaders.
                                                                   // In our case, both are the same.

  return *this;
}

BufferObjectGeometry& BufferObjectGeometry::wireAttribute (shared_ptr<FormattedVbo> source, const string& sourceAttribName) {

  return wireAttribute (sourceAttribName, source, sourceAttribName);
}

BufferObjectGeometry& BufferObjectGeometry::wireAttributes (shared_ptr<FormattedVbo> source) {

  const VertexFormat& vfd = source->getVertexFormat();

  for (int i = 0, n = vfd.getNumAttribs(); i < n; ++i) {
    wireAttribute ( source, vfd.getAttrib(i).name); // source: one formatted vertex buffer (e.g PN) => two attributes whose
	                                      // names are "aPosition" and "aNormal"
  }
  return *this;
}

BufferObjectGeometry& BufferObjectGeometry::indexedBy(shared_ptr<FormattedIbo> ib) {
  ib_ = ib;
  return *this;
}

BufferObjectGeometry& BufferObjectGeometry::primitiveType(GLenum primitiveType) {
  switch (primitiveType) {
  case GL_POINTS:
  case GL_LINE_STRIP:
  case GL_LINE_LOOP:
  case GL_LINES:
  case GL_TRIANGLE_STRIP:
  case GL_TRIANGLE_FAN:
  case GL_TRIANGLES:
  case GL_QUAD_STRIP:
  case GL_QUADS:
  case GL_POLYGON:
    break;
  default:
    assert(0);
  }
  primitiveType_ = primitiveType;
  return *this;
}

const vector<string>& BufferObjectGeometry::getVertexAttribNames() {
  if (wiringChanged_)

    processWiring(); // defined below

  return vertexAttribNames_;
}

void BufferObjectGeometry::draw(int attribIndices[]) {

  if (wiringChanged_)

    processWiring(); // defined below

  const unsigned int UNDEFINED_VB_LEN = 0xFFFFFFFF;
  unsigned int vboLen = UNDEFINED_VB_LEN;

  // bind the vertex buffers and set the vertex attribute pointers
  // for each vbo:
  for (int i = 0, n = perVbWirings_.size(); i < n; ++i) {

    const PerVbWiring &vertexWire = perVbWirings_[i];

    const VertexFormat& vertexFormat = vertexWire.vb->getVertexFormat();

    glBindBuffer(GL_ARRAY_BUFFER, *( vertexWire.vb)); // vertexWire.vb = a vertex buffer object
	checkGlErrors();
	

    vboLen = min(vboLen, (unsigned int) vertexWire.vb->length());

    for (size_t j = 0; j < vertexWire.vb2GeoIdx.size(); ++j) {

      int loc = attribIndices[ vertexWire.vb2GeoIdx[j].second];
      if (loc >= 0)

        vertexFormat.setGlVertexAttribPointer( vertexWire.vb2GeoIdx[j].first, loc); // loc = attribLocation
	    // vertexWire.vb2GeoIdx[j].first = index to an attribute
	   // => glVertexAttribPointer(glAttribLocation, ad.size,..) where loc = glAttribLocation
    } // for
  } // for

  if (isIndexed()) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *ib_);
	checkGlErrors();
	
   glDrawElements(primitiveType_, ib_->length(), ib_->getIndexFormat(), 0);
	checkGlErrors();
	 

  }
  else if (vboLen != UNDEFINED_VB_LEN) {
    glDrawArrays(primitiveType_, 0, vboLen);
	checkGlErrors();
  }
} // draw()

void BufferObjectGeometry::processWiring() {

  perVbWirings_.clear();
  vertexAttribNames_.clear();

  // maps from  vbo to index within perVbWiring_

  map< shared_ptr<FormattedVbo>, int> vbIdx;

  // go through all wiring definitions
  //  wiring_[targetAttribName] = make_pair(source, sourceAttribName); 
  
  // wiring_[] is created by  wireAttributes (vbo) in SimpleUnindexedGeometry() constructor.

  for (Wiring::const_iterator i = wiring_.begin(), e = wiring_.end(); i != e; ++i) {

    // wiring_[targetAttribute] =( source, sourceAttribute); but in our case, targetAttribute = sourceAttribute

    shared_ptr<FormattedVbo> vb = i->second.first; //  vbo 

    const VertexFormat& vertexFormat = vb->getVertexFormat(); //  vbo's vertex format = Vertex::FORMAT

    // see if this vbo is already in our vbIdx map

    map<shared_ptr<FormattedVbo>, int>::iterator j = vbIdx.find( vb );

    int idx = 0; // idx of vbo in perVbWiring_, to be set

    if (j == vbIdx.end()) {
      idx = perVbWirings_.size();
      vbIdx[vb] = idx;
      perVbWirings_.push_back( PerVbWiring( vb.get() ) ); // vb.get() returns the pointer to vb
    }
    else {
      idx = j->second;
    }
    const int globalIdx = vertexAttribNames_.size(); // idx within vertexAttribNames_, which refers to ALL the attributeNames of
	                                                 // a BufferObjectGeometry

    perVbWirings_[idx].vb2GeoIdx.push_back( make_pair( vertexFormat.getAttribIndexForName( i->second.second ), globalIdx ) );

    vertexAttribNames_.push_back( i->first ); // i->first = targetAttribute
  }

  wiringChanged_ = false; // a new wiring has been processed.
}
