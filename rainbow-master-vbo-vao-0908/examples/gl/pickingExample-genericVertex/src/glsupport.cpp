#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "glsupport.h"  

//#define GL_INFO_LOG_LENGTH 0x8B84


using namespace std;

extern std::ofstream messageFile;



void _check_gl_error( const char *file, int line) {
        GLenum err (glGetError());
 
        if (err !=GL_NO_ERROR ) {
                stringstream  errMsg;
                string error;

                switch(err) {
                        case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
                        case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
                        case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
                        case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
                        case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
                }
                //cout <<  "GL_" << error <<" - "<<file<<":"<<line<<endl;
                errMsg << "GL_" << error <<" - "<<file<<":"<<line<<endl;
               
				cout << errMsg.str() << endl;

				messageFile << errMsg.str() << endl;

				//throw runtime_error( errMsg.str() );
        }
}

//extern QOpenGLFunctions_3_0 *m_funcs;


 


// Dump text file into a character vector, throws exception on error
static void readTextFile(const char *fn, vector<char>& data) {
  // Sets ios::binary bit to prevent end of line translation, so that the
  // number of bytes we read equals file size

  static int noOfFiles = 0;

  noOfFiles ++ ;


  ifstream ifs (fn, ios::binary);

 // std::string fileName = std::string("messageFile") + std::to_string( noOfFiles ) + std::string(".txt");
  
 // ofstream messageFile (fileName);

  if (!ifs) // input file stream is not created.

    throw runtime_error(string("Cannot open file ") + fn);

  // Sets bits to report IO error using exception
  //ifs.exceptions(ios::eofbit | ios::failbit | ios::badbit);

  ifs.exceptions(  ios::failbit | ios::badbit);

  ifs.seekg(0, ios::end);  // move to the end of the file

  size_t len = ifs.tellg(); // tell me where you are; move the get pointer; seekp: put pointer

  data.resize(len);

  ifs.seekg(0,  ios::beg);  // move to the beginning of the file

  ifs.read( &data[0], len);

 

  // check if the content of the file is in good shape
  ifs.seekg( 0, ios::beg);
  //char ch;
  
  /*
  messageFile << fn << endl;
  cout << fn << endl;

  //istream& get (char* s, streamsize n,  delim);

  for (int i = 0; i < len; i++ ) {
	   cout << data[i];
	   messageFile << data[i];

  }


  messageFile.close();
  */

}

void readAndCompileSingleShader(GLuint shaderHandle, const char *fn) {
  vector<char> source;

  readTextFile(fn, source);

  const char *ptrs[] = {&source[0]}; // an array of pointers to strings
  const GLint lens[] = {source.size()}; // an array of string lengths

  GLsizei count = 1;
  // there is only one element in ptrs and lens arrays

  glShaderSource(shaderHandle, count, ptrs, lens);   // load the shader sources

  glCompileShader(shaderHandle); // normal-gl3.fshader caused compile errors "unhandled viloation"


  GLint compiled = 0;
  glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &compiled);

  if (!compiled) {

     GLint maxLength = 0;
	 
	 //#define GL_INFO_LOG_LENGTH 0x8B84
     glGetShaderiv(shaderHandle, GL_INFO_LOG_LENGTH, &maxLength);

	 //The maxLength includes the NULL character
	 char * errorLog = new char[maxLength];

	 glGetShaderInfoLog(shaderHandle, maxLength, &maxLength, errorLog);
	  checkGlErrors();

	 /*
	 char * charArray = new char [ maxLength ];
	 for ( int i =0; i < maxLength; i++ ) {
		  charArray[i] = errorLog[i];
	 }
	 */

	 string  errorMsg;
	 
	 char errorType [80];

	 sprintf( errorType, "%s %s: ", "Compile failure in shader", fn);

	 errorMsg = string(errorType ) +  string(errorLog );

	 /*
	 for ( int i = 0; i < maxLength;i++ ) {
	      messageFile << errorLog[i] << endl;
	 }
	 */

	

	 messageFile << errorMsg << endl;
	 std::cout << errorMsg << endl;



	 glDeleteShader( shaderHandle );

     throw runtime_error( errorMsg );

  }
}

void linkShader(GLuint programHandle, GLuint vs, GLuint fs, const char * vsFilename, const char * fsFilename) {

  glAttachShader(programHandle, vs);
  glAttachShader(programHandle, fs);

  glLinkProgram(programHandle);

  /* 
  As a result of a successful link operation, all active user-defined uniform variables belonging to program will 
  be initialized to 0, and each of the program object's active uniform variables will be assigned a location 
  that can be queried by calling glGetUniformLocation. Also, any active user-defined attribute variables 
  that have not been bound to a generic vertex attribute index will be bound to one at this time.
  
  When a program object has been successfully linked, 
  the program object can be made part of current state by calling glUseProgram
  */

  glDetachShader(programHandle, vs);
  glDetachShader(programHandle, fs);

  GLint linked = 0;
  glGetProgramiv(programHandle, GL_LINK_STATUS, &linked);
 

  if (!linked) {

     GLint maxLength = 0;
	 
	 //#define GL_INFO_LOG_LENGTH 0x8B84
     glGetShaderiv(programHandle, GL_INFO_LOG_LENGTH, &maxLength);

	 //The maxLength includes the NULL character

	  //The maxLength includes the NULL character
	 char * errorLog = new char[maxLength];
	
	 glGetProgramInfoLog(programHandle, maxLength, &maxLength, errorLog);
	 	 checkGlErrors();
	 // printF => F: formatted: printf() = fprintf(stdout, ...)
	 
	 string  errorMsg;
	 
	 char errorType [80];

	 sprintf( errorType, "%s %s + %s: ", "link failure in program ", vsFilename, fsFilename);


	 errorMsg = string(errorType ) +  string(errorLog );

	 /*
	 for ( int i = 0; i < maxLength;i++ ) {
	      messageFile << errorLog[i] << endl;
	 }
	 */



	 messageFile << errorMsg << endl;
	 std::cout << errorMsg << endl;

	 throw runtime_error( errorMsg );
	 

  }
}


void readAndCompileShader(GLuint programHandle, const char * vertexShaderFileName, const char * fragmentShaderFileName) {
  GlShader vs(GL_VERTEX_SHADER);
  GlShader fs(GL_FRAGMENT_SHADER);

  try {
    readAndCompileSingleShader(vs, vertexShaderFileName);
  }
  catch (const runtime_error & error ) {
	  std::cout << error.what() << endl;
	  messageFile << "comiple shader:" << error.what() << endl;

  }
  
  try {
    readAndCompileSingleShader(fs, fragmentShaderFileName);
  } 
  catch (const runtime_error & error ) {
	  std::cout << error.what() << endl;
	   messageFile <<"comiple shader:" << error.what() << endl;
  }

  try {
    linkShader(programHandle, vs, fs, vertexShaderFileName, fragmentShaderFileName);
  }
  catch (const runtime_error & error ) {
	  std::cout << error.what() << endl;
	   messageFile << "comiple shader:" << error.what() << endl;
	   throw;
  }

}