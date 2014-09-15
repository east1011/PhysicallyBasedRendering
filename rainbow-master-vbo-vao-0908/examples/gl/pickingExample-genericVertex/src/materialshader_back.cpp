#include <cassert>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>

#include "glsupport.h"

#include "uniforms.h"
#include "materialshader.h"
#include "shadergeometry.h"

//#include "ofMain.h"


//extern QOpenGLFunctions_3_0 *m_funcs;
extern bool g_Gl2Compatible;

using namespace std;
using namespace tr1;

struct GlProgramDesc {

  struct UniformDesc {

    string name;
    GLenum type;
    GLint size; // number of elements (e.g. matrices)
    GLint location;
  };

  struct AttribDesc {
    string name;
    GLenum type;
    GLint size;
    GLint location;
  };


  GlProgram program; // => default constructor calls _handle = glCreateProgram();

  vector<UniformDesc> uniforms;
  vector<AttribDesc> attribs;

  // constructor  of GlProgramDesc: query the shader program and retrieve the information about 
  // the uniform variables and attribute variables, and store it in uniforms and attribs members.
  // This information will be used when drawing, to bind the data provided by the user to the appropriate
  // uniform and attribute variables of the shader program.


  GlProgramDesc(GLuint vsHandle, GLuint fsHandle, string vsFilename, string fsFilename)  {

    try {
		linkShader(program, vsHandle, fsHandle, vsFilename.c_str(), fsFilename.c_str()); // program object is type-converted to handle_
	}
	catch (const runtime_error & error ) {
     messageFile << error.what() << endl;
	 std::cout << error.what() << endl;
	 throw;

	}

    int numActiveUniforms, numActiveAttribs, uniformMaxLen, attribMaxLen;

    glGetProgramiv(program, GL_ACTIVE_UNIFORMS, &numActiveUniforms);

	//An attribute variable (either built-in or user-defined) is considered active if it is determined during the link operation
	// that it may be accessed during program execution. Therefore, program should have previously been the target 
	// of a call to glLinkProgram, but it is not necessary for it to have been linked successfully.

    glGetProgramiv(program, GL_ACTIVE_ATTRIBUTES, &numActiveAttribs);

    glGetProgramiv(program, GL_ACTIVE_UNIFORM_MAX_LENGTH, &uniformMaxLen);
    glGetProgramiv(program, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &attribMaxLen);

    const int bufSize = max(uniformMaxLen, attribMaxLen) + 1; //  the size of the string name of the uniform or attribute name;
	                                                          //  1 is added because the string terminates with null char.

    vector<GLchar> buffer(bufSize); // buffers containing the name of the uniform or attribute variable.

    uniforms.resize(numActiveUniforms);

	// void glGetActiveUniform( GLuint program,  GLuint index,  GLsizei bufSize,  	GLsizei *length,  GLint *size, GLenum *type, GLchar *name); 

    for (int i = 0; i < numActiveUniforms; ++i) { // i = index of the uniform variable to be queried

      GLsizei charsWritten;  // length of the buffer containing the name of uniform

      glGetActiveUniform(program, i, bufSize, &charsWritten, & (uniforms[i].size), & (uniforms[i].type), &buffer[0]);
      
	  // type: The type argument will return a pointer to the uniform variable's data type.
	  // The symbolic constants GL_FLOAT, GL_FLOAT_VEC2, GL_FLOAT_VEC3, GL_FLOAT_VEC4, GL_INT, GL_INT_VEC2, 
	  // GL_INT_VEC3, GL_INT_VEC4, GL_BOOL, GL_BOOL_VEC2, GL_BOOL_VEC3, GL_BOOL_VEC4, GL_FLOAT_MAT2, 
	  // GL_FLOAT_MAT3, GL_FLOAT_MAT4, GL_SAMPLER_2D, GL_SAMPLER_3D, or GL_SAMPLER_CUBE, etc.
	  
	  assert(charsWritten + 1 <= bufSize);

      uniforms[i].name = string(buffer.begin(), buffer.begin() + charsWritten);
      
	  uniforms[i].location = glGetUniformLocation(program, &buffer[0]); // buffer: shader variable name
    }

    attribs.resize(numActiveAttribs);

    for (int i = 0; i < numActiveAttribs; ++i) {
      GLsizei charsWritten;
      glGetActiveAttrib(program, i, bufSize, &charsWritten, &(attribs[i].size), &(attribs[i].type), &buffer[0]); // buffer = name of an attribute

      assert(charsWritten + 1 <= bufSize);
      attribs[i].name = string(buffer.begin(), buffer.begin() + charsWritten);
      attribs[i].location = glGetAttribLocation(program, &buffer[0]);
    }

	// As with uniform variables and vertex attributes, the output variables of the fragment shader 
	// also have a location. In the shader above, the variable colorOut is bound to location 0 (zero) 
	// by default. In the application we should bind the variable colorOut to output location zero, 
	// the default location when not using framebuffer objects, using the following function:
	// void glBindFragDataLocation( GLuint program,    GLuint colorNumber,   const char * name);
	// colorNumber <= GL_MAX_DRAW_BUFFERS; The bindings specified by glBindFragDataLocation have no effect 
	// until program is next linked. Bindings may be specified at any time after program has been created. 

	// Varying out varyings may have indexed locations assigned explicitly in the shader text using 
	// a location layout qualifier. If a shader statically assigns a location to a varying out variable 
	// in the shader text, that location is used and any location assigned with glBindFragDataLocation is ignored. 



 


    if (!g_Gl2Compatible) {
		// index 0, 1, 2 refers to the index of the array attachments in:
		// GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2 }; 
		// glDrawBuffers(3,  attachments);

      glBindFragDataLocation(program, 0, "fragColor");
	  glBindFragDataLocation(program, 1, "phaseSpectrum");
	  glBindFragDataLocation(program, 2, "rainbowColorSpectrum");
	}

    checkGlErrors();
  } //GlProgramDesc()


}; // GlprogamDesc;


class GlProgramLibrary {
  typedef map< pair<string, GLenum>, shared_ptr<GlShader> > GlShaderMap; // GlShader => shader handle
  typedef map< pair<string, string>, shared_ptr<GlProgramDesc> > GlProgramDescMap;

  GlShaderMap shaderMap;
  GlProgramDescMap programMap;

  GlProgramLibrary() {}

public:
  // static method = class method

  static GlProgramLibrary& getSingleton() {
    static GlProgramLibrary pl;
    return pl;
  }

  shared_ptr<GlProgramDesc> getProgramDesc(const string& vsFilename, const string& fsFilename) {

    GlProgramDescMap::key_type key(vsFilename, fsFilename);


	GlProgramDescMap::iterator i = programMap.find(key);

    if (i == programMap.end()) { // no program exists with the same vshader and fshader

      shared_ptr<GlProgramDesc> program( new GlProgramDesc( *getShader(vsFilename, GL_VERTEX_SHADER), 
		      *getShader(fsFilename, GL_FRAGMENT_SHADER) , vsFilename, fsFilename ) );

	  // getShader() compiles the shaders and GlProgramDesc() links the shaders

      programMap[key] = program;
      return program;
    }
    else {
      return i->second;
    }
  }

protected:
  shared_ptr<GlShader> getShader(const string& filename, GLenum shaderType) {
    string f = filename;
    if (g_Gl2Compatible) { // if compatible, change -gl3 to -gl2 in the end of the filename
      size_t pos = f.rfind("-gl3");
      if (pos != string::npos) {
        f[pos+3] = '2';
      }
    }

    GlShaderMap::key_type key(f, shaderType);

    GlShaderMap::iterator i = shaderMap.find(key);
    if (i == shaderMap.end()) {
      shared_ptr<GlShader> shader(new GlShader(shaderType));

	  try {
         readAndCompileSingleShader(*shader, f.c_str());
	  }
	  catch (const runtime_error & error ) {
	     std::cout << error.what() << endl;
	     messageFile << error.what() << endl;
	  }
    
      shaderMap[key] = shader;
      return shader;
    }
    else {
      return i->second;
    }
  }
}; //GlProgramLibrary 

MaterialShader::MaterialShader(const string shaderName, const string& vsFilename, const string& fsFilename)

  : shaderName ( shaderName),  programDesc_ ( GlProgramLibrary::getSingleton().getProgramDesc(vsFilename, fsFilename) ) 
  // set member programDesc_ of Material class.
  // shared_ptr<GlProgramDesc> programDesc_. which contains members: GlProgram program; vector<UniformDesc> uniforms; 
  // vector<AttribDesc> attribs;
  // programDesc_ may refer to the same program (vshader + fshader). But each MaterialShader has different programDesc_'s which
  // refer to the same program. In this case, different uniform and attribute values will be assgined to the same program. 
{}

static const char * getGlConstantName(GLenum c) {
  struct ValueNamePair {
    GLenum value;
    const char * name;
  };

  ValueNamePair valueNamePairs[] = {
    { GL_FLOAT, "GL_FLOAT" },
    { GL_FLOAT_VEC2, "GL_FLOAT_VEC2" },
    { GL_FLOAT_VEC3, "GL_FLOAT_VEC3" },
    { GL_FLOAT_VEC4, "GL_FLOAT_VEC4" },
    { GL_FLOAT_MAT2, "GL_FLOAT_MAT2" },
    { GL_FLOAT_MAT3, "GL_FLOAT_MAT3" },
    { GL_FLOAT_MAT4, "GL_FLOAT_MAT4" },
    { GL_FLOAT_MAT2x3, "GL_FLOAT_MAT2x3" },
    { GL_FLOAT_MAT2x4, "GL_FLOAT_MAT2x4" },
    { GL_FLOAT_MAT3x2, "GL_FLOAT_MAT3x2" },
    { GL_FLOAT_MAT3x4, "GL_FLOAT_MAT3x4" },
    { GL_FLOAT_MAT4x2, "GL_FLOAT_MAT4x2" },
    { GL_FLOAT_MAT4x3, "GL_FLOAT_MAT4x3" },
    { GL_INT, "GL_INT" },
    { GL_INT_VEC2, "GL_INT_VEC2" },
    { GL_INT_VEC3, "GL_INT_VEC3" },
    { GL_INT_VEC4, "GL_INT_VEC4" },
    { GL_UNSIGNED_INT_VEC2, "GL_UNSIGNED_INT_VEC2" },
    { GL_UNSIGNED_INT_VEC3, "GL_UNSIGNED_INT_VEC3" },
    { GL_UNSIGNED_INT_VEC4, "GL_UNSIGNED_INT_VEC4" },
    { GL_SAMPLER_1D, "GL_SAMPLER_1D" },
    { GL_SAMPLER_2D, "GL_SAMPLER_2D" },
	 { GL_SAMPLER_3D, "GL_SAMPLER_3D" },
    { GL_SAMPLER_CUBE, "GL_SAMPLER_CUBE" },
    { GL_SAMPLER_1D_SHADOW, "GL_SAMPLER_1D_SHADOW"},
    { GL_SAMPLER_2D_SHADOW, "GL_SAMPLER_2D_SHADOW"},
  };

  for (int i = 0, n = sizeof(valueNamePairs)/sizeof(valueNamePairs[0]); i < n; ++i) {
    if (valueNamePairs[i].value == c)
      return valueNamePairs[i].name;
  }
  return "Unkonwn";
}

void MaterialShader::draw(Geometry& geometry, const Uniforms& extraUniforms) {

  static GLint maxTextureImageUnits = 0;
  	// create and bind a vertex array object (vao) which is required in the core project of GL
    // We now use the core profile set up by openframeworks
	// The following added after a week's debugging by Moon Jung, 2014/8/9
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray( vao );

  // Initialize maxTextureImageUnits if this is called for the first time
  if (maxTextureImageUnits == 0) {

    glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &maxTextureImageUnits);
    assert(maxTextureImageUnits > 0); // GL spec says this has to be at least 2
  }

 // install a program object as part of current rendering state
// glUseProgram means that the given program object is the current program that will be
//	used for things that use programs (glUniform, rendering commands, etc). 
//	0 is a lot like NULL for OpenGL objects; it represents "not an object". 
//	Therefore, glUseProgram means that no program is current, 
//	and therefore no program will be used for things that use programs.
//  
 //  Note that in the core profile, there is no fixed-function pipeline, 
//  so you'll just get an OpenGL error (GL_INVALID_OPERATION) 
//  –  Robert Rouhani Nov 24 '12 at 23:41 
  glUseProgram( programDesc_->program ); // this->programDesc_, this->uniforms_
   try {
		 checkGlErrors();
	 }
	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

  renderStates_.apply();  // transit to current states

  // Step 1:
  // set the uniforms and bind the textures

  int currTextureUnit = 0;  // textureUnit accumulates while visiting the texture uniform variables, and a single
                        // texture uniform variable may refer to an array of texture units

  // programDesc_->uniforms = the uniform variables OF the shader programs
  for (int i = 0, n = programDesc_->uniforms.size(); i < n; ++i) { 
	    // programDesc_->uniforms refers to the uniforms  active in the SHADER  program.
	    

    // ud = uniforms[i], where uniforms[i].name and uniforms[i].location

    const GlProgramDesc::UniformDesc& shaderUniform = programDesc_->uniforms[i]; // shader program's uniforms

    const Uniforms* userUniformsList[] = { &uniforms_,  &extraUniforms}; // ALL the uniform variables SET BY THE USERser program
	
	 // uniformsList = uniforms provided by the USER program (main program)

	// check if the two lists match

    int j = 0;
    for (; j < 2; ++j) {  // uniformsList[0] = uniforms_, and uniformsList[1] = extraUniforms
    
      const Uniforms::Value* userUniformValue = userUniformsList[j]->get( shaderUniform.name); 
	      // get() returns the values of ud.name  from the uniforms list "put" by the USERER
	      // It is NULL if ud.name is not found in the list.
	                               
      // if the name looks like x[0], and the uniform is not found, we also try stripping the '[0]'
      if (userUniformValue == NULL && shaderUniform.name.length() >= 3 && shaderUniform.name.compare( shaderUniform.name.length() - 3,  3,  "[0]") == 0)
		                                              // compare(pos, size, str)
        userUniformValue = userUniformsList[j]->get( shaderUniform.name.substr( 0, shaderUniform.name.length() - 3) );

      if (userUniformValue) { // userUniformValue is a user uniform variable which matches the current shader uniform variable i

        if (userUniformValue->type == shaderUniform.type && userUniformValue->size <=  shaderUniform.size) { 
			// user variable's value u and shader variable's value  are of the same type
			// and user value's  size is  shorter  than the size allocated to the shader uniform variable. 
			// ==> This is a happy situation

          switch (userUniformValue->type) { // which is also the type of the shader value uShader

			  // in case of textures:
          case GL_SAMPLER_1D:
          case GL_SAMPLER_2D: // e.g. uniform sampler2D uTexColor;
          case GL_SAMPLER_3D: // e.g. uniform sampler3D uShaderistTex, uPhaseTex;
          case GL_SAMPLER_CUBE:
          case GL_SAMPLER_1D_SHADOW:
          case GL_SAMPLER_2D_SHADOW:
            {
              const shared_ptr<ShaderTexture> *textures = userUniformValue->getTextures(); 
			   // textures is a polymorphic pointer which points to the subclasses of ShaderTexture

			  
			    // since userUniformValue->type is a texture, userUniformValue is TexturesValue
			   // tex is  named texture(s) which has been created by glGenTexture()

              // If this assert hits, the Uniform::Value is incorrectly implemented
              assert(textures != NULL);

			  if ( textures == NULL ) {
				  cout  << "texture is NULL " << endl;
				  messageFile << "texture is NULL " << endl;
				  throw;
			  }

              static const int MAX_TEX_UNITS = 1024;

              GLint boundTexUnits[MAX_TEX_UNITS];
              int texIndex = 0; // refers to the texutre units that belong to a SINGLE uniform variable
       
			  // bind all the textures to the texture units incrementing from GL_TEXTURE0

			  for (; texIndex < userUniformValue->size; ++texIndex) { // userUniformValue is an  array of texture units

                if ( currTextureUnit == maxTextureImageUnits) { // local variable textureUnit was initialized to 0
                  stringstream s;

				  // this->shaderName
                  s << shaderName << ": " << "System allows a maximum of " << maxTextureImageUnits << ". The current shader is trying to use more than that.";
                  throw runtime_error(s.str());
                }

				// The texture  bind locations are named  GL_TEXTURE0​, GL_TEXTURE1​, etc. Alternatively, you can use bind locations of the form GL_TEXTURE0​ + i, 
				//where i is a texture unit number between 0 and the implementation-defined constant 
				//GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS​.

				// you may wish to use more than one texture at the same time to achieve certain effects such as bump mapping 
				// or gloss mapping. To do this generate, bind, set parameters, and load data into multiple texture names, 
				//  then you bind them to multiple texture image units and send their unit numbers to the shader like this:

                //
                // glActiveTexture(GL_TEXTURE0); //switch to texture image unit 0
                // glBindTexture(GL_TEXTURE_2D, textures[0]);   //Bind a texture to this unit
                // glActiveTexture(GL_TEXTURE1); //switch to texture image unit 1
                // glBindTexture(GL_TEXTURE_2D, textures[1]);   //Bind a different texture to this unit

                

                glActiveTexture(GL_TEXTURE0 + currTextureUnit);  
				  checkGlErrors();

                textures[ texIndex]->bind(); // make the texture indexed by texIndex bound to the currently active texture unit
				//  bind method calls  glBindTexture(GL_TEXTURE_2D, tex[texIndex]); 

				// You can have multiple textures loaded into your Rendering Context.
				                
				boundTexUnits[ texIndex] = currTextureUnit++; // assign currTextureUnit to the list of actviated texture units and increment currTextureUnit for another texture;
				                                  // This increment accumulates over all the texture units used in a single USER  program.
				
              } // for 

			  // send  the array of texture units  boundTexUnits provided by
			  // the user to the uniform variable "shaderUniform.location",  defined in the shader

              userUniformValue->apply(shaderUniform.location, userUniformValue->size, boundTexUnits); // upload the uniform variable which may refer to 
			                                            // an array of texture units, texUnits.
			  

			  // texUnits = an array of  texture units
			  // apply() => genericGlUniformv(location, userUniformValue->size, boundTexUnits);
            } // case of textures

            break;  // user uniform variable is a texture type

          default:   // uniform variables: other than texture types

            userUniformValue->apply(shaderUniform.location,  userUniformValue->size, NULL); // 

          } // switch
		  
		  //break; // the uniform variable has been found and processed; no need to loop over the other uniform variables.


        } //  if (the current shader uniform variable  matches a user  uniform variable in type and size

        else { // userUniformValue and  shader uniform variable do not  match in type and size, although their names are equal
          stringstream s;
          s << shaderName << ": " << "Uniform variable " << shaderUniform.name << ": supplied value and declared variable do not match in type and/or size."
            << "\nSupplied value: type = " << getGlConstantName(userUniformValue->type) << ", size = " << userUniformValue->size
            << "\nDeclared in shader: type = " << getGlConstantName( shaderUniform.type) << ", size = " << shaderUniform.size;
          throw runtime_error(s.str());

		  // break;

        }

        break; // the the current shader uniform variable has been found in the user uniform variables and processed; no need to loop over the other uniform variables.
		       // // once the name has been found, break from the search loop (j loop)

      } // if (the uniform variable  is provided  ) 

	  // userUniformValue is NULL or not

    } // j loop (j=0,1)

	/* for debuggig: comment out
    if (j == 2) { // Unless the j loop is broken, it means that the j loop has failed to search the current shader uniform variable in the user uniform variable list
		
      stringstream s;
      s << shaderName << ": " << "Uniform variable " << shaderUniform.name << ": used in the shader codes, but not supplied. Its Type = " << getGlConstantName(shaderUniform.type) << ", Size = " << shaderUniform.size;
      throw runtime_error(s.str()); // This error will be catched after draw() method, but is processed simply by writing the error message

    }

	*/

  } // i: for all uniform variables of the shader program.

  // Step 2:
  // see what attribs are provided by the geometry

  const vector<string>& geoAttribNames = geometry.getVertexAttribNames();

      //      calls processWiring(); return  vertexAttribNames_ which contains the vertex attribute names.


  const static int MAX_ATTRIB = 64;
  int attribIndices[MAX_ATTRIB];
  const size_t numAttribs = geoAttribNames.size();
  


  if (numAttribs > MAX_ATTRIB) {
	  stringstream s;
	  s << shaderName <<": " << "Number of attributes contained in geometry is greater than maximally supported number of attributes. Consider increasing MAX_ATTRIB.";

      throw runtime_error( s.str() );
  }

  for (size_t i = 0; i < numAttribs; ++i) {
    attribIndices[i] = -1;
  }

  // simple and stupid O(n^2) wiring, should use a hashtable to reduce to O(n)
  for (int i = 0, n = programDesc_->attribs.size(); i < n; ++i) {
    const GlProgramDesc::AttribDesc& ad = programDesc_->attribs[i];

    size_t j = 0;
    for (; j < numAttribs; ++j) { // numAttribs from geometry
      if (geoAttribNames[j] == ad.name) { // ad from shader
        attribIndices[j] = ad.location;
        break;
      }
    }

	// comment out for debugging, 2014/5/25

    if (j == numAttribs) {
      stringstream s;
	  s << shaderName << ": " << string("Vertex attribute ") + ad.name
                          + ": used in the shader codes, but not supplied." ;
      throw runtime_error(s.str() );
	  //for debugging purpose 14.04.22
    }

	

  }

  for (size_t i = 0; i < numAttribs; ++i) {
    if (attribIndices[i] >= 0) {  // for the active attributes whose values are given
       glEnableVertexAttribArray(attribIndices[i]);
  

       try {
		 checkGlErrors();
		 
	   }

	   catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled
		// GL_INVALID_ENUM
        // An unacceptable value is specified for an enumerated argument. 
		// The offending command is ignored and has no other side effect than to set the error flag. 

		//GL_INVALID_OPERATION
        // The specified operation is not allowed in the CURRENT STATE. The offending command is ignored and has no other side effect than to set the error flag. 


	  }

	} // if
  } // for

  // Now let the geometry draw its self
  
  //glUseProgram( programDesc_->program );

  geometry.draw(attribIndices); // upload the attribute data and call opengl Draw() function

  for (size_t i = 0; i < numAttribs; ++i) {
    if (attribIndices[i] >= 0) { 
      glDisableVertexAttribArray(attribIndices[i]);
 

     try {
		 checkGlErrors();
	 }
	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }
  } // if
 } // for

  glBindVertexArray( vao );

} // MaterialShader::draw
