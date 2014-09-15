#ifndef UNIFORMS_H
#define UNIFORMS_H

#include <map>
#include <vector>
#include <memory>
#include <stdexcept>
#include <string>
#if __GNUG__
#   include <tr1/memory>
#endif


#include "matrix4.h"
#include "glsupport.h"
#include "shadertexture.h"

//extern QOpenGLFunctions_3_0 *m_funcs;

// Private namespace for some helper functions. You should ignore this unless you
// are interested in the internal implementation.
namespace _helper {
inline void genericGlUniformi(GLint location, int i) {
:: glUniform1i(location, i);  // bind texture unit
  checkGlErrors();

}
inline void genericGlUniformf(GLint location, float f) {
  ::glUniform1f(location, f);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const GLint *v) {
  ::glUniform1iv(location, size, v);  // bind texture units; size is more than 1, uniform variable
    checkGlErrors();
  // at location should be an array
}
inline void genericGlUniformv(GLint location, int size, const GLfloat *v) {
  ::glUniform1fv(location, size, v);
}
inline void genericGlUniformv(GLint location, int size, const Cvec<int, 1> *v) {
  ::glUniform1iv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<int, 2> *v) {
  ::glUniform2iv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<int, 3> *v) {
  ::glUniform3iv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<int, 4> *v) {
  ::glUniform4iv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<float, 1> *v) {
  ::glUniform1fv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<float, 2> *v) {
  ::glUniform2fv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<float, 3> *v) {
  ::glUniform3fv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformv(GLint location, int size, const Cvec<float, 4> *v) {
  ::glUniform4fv(location, size, &v[0][0]);
    checkGlErrors();
}
inline void genericGlUniformMatrix4v(GLint location, int size, const Cvec<float, 16> *m) {
  ::glUniformMatrix4fv(location, size, GL_FALSE, &m[0][0]);
  checkGlErrors();
}

template<typename T, int n>
inline GLenum getTypeForCvec();   // should replace with STATIC_ASSERT

template<>
inline GLenum getTypeForCvec<int, 1>() { return GL_INT; }
template<>
inline GLenum getTypeForCvec<int, 2>() { return GL_INT_VEC2; }
template<>
inline GLenum getTypeForCvec<int, 3>() { return GL_INT_VEC3; }
template<>
inline GLenum getTypeForCvec<int, 4>() { return GL_INT_VEC4; }
template<>
inline GLenum getTypeForCvec<float, 1>() { return GL_FLOAT; }
template<>
inline GLenum getTypeForCvec<float, 2>() { return GL_FLOAT_VEC2; }
template<>
inline GLenum getTypeForCvec<float, 3>() { return GL_FLOAT_VEC3; }
template<>
inline GLenum getTypeForCvec<float, 4>() { return GL_FLOAT_VEC4; }
template<>
inline GLenum getTypeForCvec<bool, 1>() { return GL_BOOL; }
template<>
inline GLenum getTypeForCvec<bool, 2>() { return GL_BOOL_VEC2; }
template<>
inline GLenum getTypeForCvec<bool, 3>() { return GL_BOOL_VEC3; }
template<>
inline GLenum getTypeForCvec<bool, 4>() { return GL_BOOL_VEC4; }
}

// The Uniforms keeps a map from strings to values
//
// Currently the value can be of the following type:
// - Single int, float, or Matrix4
// - Cvec<T, n> with T=int or float, and n = 1, 2, 3, or 4
// - shared_ptr<Texture>
// - arrays of any of the above
//
// You either use uniform.put("varName", val) or
// uniform.put("varArrayName", vals, numVals);
//
// A Uniforms instance will start off empty, and you can use
// its "put" member function to populate it.


//  UniformDesc is different from Uniforms class

class Uniforms {

public:

  Uniforms& put(const std::string& name, int value) { // put() is the "reference" to some value, not the value
	                                                  // itself or the pointer to the value

    Cvec<int, 1> v ( value );  // v: local vector with one element

    valueMap[name].reset( new CvecsValue<int, 1> (&v, 1) );  // one one-vector 
	// AHA: CvecsValue () constructor is called in Uniforms::put()

    return *this;
  }

  Uniforms& put(const std::string& name, float value) {
    Cvec<float, 1> v( value );

    valueMap[name].reset(new CvecsValue<float, 1>(&v, 1));
    return *this;
  }

  Uniforms& put(const std::string& name, const Matrix4& value) {
    valueMap[name].reset(new Matrix4sValue(&value, 1));
    return *this;
  }

  Uniforms& put( const std::string& name, const std::tr1::shared_ptr<ShaderTexture>& value ) {
	 

    valueMap[name].reset(new TexturesValue(&value, 1) );
    return *this;
  }

  template<int n>

  Uniforms& put(const std::string& name, const Cvec<int, n>& v) {
    valueMap[name].reset(new CvecsValue<int, n> (&v, 1));  // one n-vector
    return *this;
  }

  template<int n>
  Uniforms& put(const std::string& name, const Cvec<float, n>& v) {
    valueMap[name].reset(new CvecsValue<float, n>(&v, 1));
    return *this;
  }

  template<int n>
  Uniforms& put(const std::string& name, const Cvec<double, n>& v) {
    Cvec<float, n> u;
    for (int i = 0; i < n; ++i) {
      u[i] = float(v[i]);
    }
    valueMap[name].reset(new CvecsValue<float, n>(&u, 1));  // one n- float vector, because glsl 1.30 allows only floats

    return *this;
  }

  Uniforms& put(const std::string& name, const int *values, int count) {  // an array input
    valueMap[name].reset(new CvecsValue<int, 1>( reinterpret_cast<const Cvec<int, 1>*>(values), count) );
	 // an integer is regarded as Cvec<int,1>, a vector of one element
    return *this;
  }

  Uniforms& put(const std::string& name, const float *values, int count) {
    valueMap[name].reset(new CvecsValue<int, 1>(reinterpret_cast<const Cvec<float, 1>*>(values), count));
    return *this;
  }

  Uniforms& put(const std::string& name, const Matrix4 *values, int count) {
    valueMap[name].reset(new Matrix4sValue(values, count));
	
	// Matrix4sValue() converts a row-major matrix to a column-major

    return *this;
  }

  Uniforms& put(const std::string& name, const std::tr1::shared_ptr<ShaderTexture> *values, int count) {
    valueMap[name].reset(new TexturesValue(values, count));
    return *this;
  }

  template<int n>
  Uniforms& put(const std::string& name, const Cvec<int, n> *v, int count) {
    valueMap[name].reset(new CvecsValue<int, n>(v, count));
    return *this;
  }

  template<int n>
  Uniforms& put(const std::string& name, const Cvec<float, n> *v, int count) {
    valueMap[name].reset(new CvecsValue<float, n>(v, count));
    return *this;
  }


  template<int n>
  Uniforms& put(const std::string& name, const Cvec<double, n> *v, int count) {
    valueMap[name].reset(new CvecsValue<float, n>(v, count));
    return *this;
  }

  // Future work: add put for different sized matrices, and array of basic types
//protected: just for temporary debugging 2013/12/06

  // Ghastly implementation details follow. Viewer be warned.

  // MEMBER variables of Uniforms class, used above

  friend class Material; // Material can access private members of Uniform classs as its friend
  class ValueHolder; // forward decl, has value_ as the pointer to some value

  class Value;

  typedef std::map<std::string, ValueHolder> ValueMap;

  ValueMap valueMap;

  const Value* get(const std::string& name) const {

    std::map< std::string, ValueHolder >::const_iterator i = valueMap.find(name);

    return  i == valueMap.end() ? NULL : i->second.get(); // second==ValueHolder
  }

  class ValueHolder {
    Value *value_;

  public:
    ValueHolder() : value_(NULL) {}
    ValueHolder(Value* value) : value_(value) {}
    ValueHolder(const ValueHolder& u) : value_( u.value_ ? u.value_->clone() : NULL) {}
    ~ValueHolder() {
      if (value_)
        delete value_;
    }
    void reset(Value* value) {
      if (value_)
        delete value_;
      value_ = value;
    }

    Value *get() const {
      return value_;
    }

    ValueHolder& operator= (const ValueHolder& u) {
      reset(u.value_ ? u.value_->clone() : NULL);
      return *this;
    }
  }; // class ValueHolder


  class Value {
  public:
    
    const GLenum type;

    // 1 for non-array type, otherwise the number of elements in the array
    const GLint size;

    virtual Value* clone() const = 0;
    virtual ~Value() {}

    // If type is one of GL_SAMPLER_*, the method apply should use the "boundTexUnits" argument
	// as the argument for glUniform*.
    //
    // Otherwise, boundTexUnit should be ignored and whatever values contained in
    // the Value instance should be set to given "location" by glUniform*
    //
    // `count' specifies how many actural uniforms are specified by the shader, and
    // should be used as input parameter to glUniform*

    virtual void apply(GLint location,  GLsizei count, const GLint *boundTexUnits) const = 0; // virtual

	// If type is one of GL_SAMPLER_*, method getTextures()  should provide a pointer to
    // the array of shared_ptr<Texture> stored by the Uniform. 

    virtual const std::tr1::shared_ptr<ShaderTexture> * getTextures() const { return NULL; }; 
	       // reimplemented in subclasses

  protected:  // constructor
    Value(GLenum aType, GLint aSize) : type(aType), size(aSize) {}

  }; // class Value

  template<typename T, int n>

  class CvecsValue : public Value {
  
	  std::vector<Cvec<T, n> >  vs_;   // private member

  public:

    CvecsValue(const Cvec<T, n> *vs, int size)

      : Value(_helper::getTypeForCvec<T,n>(), size), vs_(vs, vs + size) {

      assert(size > 0);
    }

    // construct from Cvecs of another type

    template<typename S>

    CvecsValue(const Cvec<S, n> *vs, int size)
      : Value( _helper::getTypeForCvec<T,n>(), size), vs_(size) {
      assert(size > 0);
      for (int i = 0; i < size; ++i) { // size is the member of CvecsValue
        for (int d = 0; d < n; ++d) {
          vs_[i][d] = T(vs[i][d]);   
        }
      }
    }

    virtual Value* clone() const {
      return new CvecsValue(*this);
    }

    virtual void apply(GLint location, GLsizei count, const GLint *boundTexUnit) const {  // boundTextUnit not used here

      assert( count <= size);
      _helper::genericGlUniformv(location, count, &vs_[0]);
    }
  };  // class CvecsValue

  class Matrix4sValue : public Value {
    // we use cvecs here instead of matrix4s here since matrix4 (in glsl)
    // doesn't allow double element (yet), and to pass the data
    // into glUniformMatrix4fv, we need to have the internal buffer
    // to be typed float.

    std::vector< Cvec<float, 16> > ms_;

  public:

    Matrix4sValue(const Matrix4 *m, int size)
      : Value( GL_FLOAT_MAT4, size ), ms_ ( size )
    {
      assert(size > 0);
      for (int i = 0; i < size; ++i)

        m[i].writeToColumnMajorMatrix(&ms_[i][0]);
    }

    virtual Value* clone() const {
      return new Matrix4sValue(*this); // a copy constructor (default) is used
    }

    virtual void apply(GLint location, GLsizei count, const GLint *boundTexUnit) const { // boundTexUnit not used here

      assert( count <= size ); // Value::size
      _helper::genericGlUniformMatrix4v(location, count, &ms_[0]);
    }

  };  // class  Matrix4sValue 

  // vector copy constructor:
//  vector <T> v(otherVector); or
//vectpr <T> v = otherVector; (alternate usage syntax).

// vector<T> v(inIterBegin, inIterEnd);
  // Construct v containing values from the range [inIterBegin,inIterEnd) in another (not necessarily vector) container, 
  // but whose component type is the same as the component type of v.
  /*
  The vector constructor (and all STL constructors that accept pointer ranges, for that matter) are designed to take
  in a range of STL-style iterators. When using iterators, you specify a range by providing a pointer to the first- 
  and the past-the-end elements, not the first and last elements. If you wanted to create a vector as a copy of the subrange (0, 99)
  out of another vector, you could write

      vector<some> dest(source.begin(), source.begin() + 100);

   Note that this uses vector iterators to specify the slice of the first 100 elements rather than operator[], 
   
   which has undefined behavior when the provided index is out of bounds. 
   
   If you want to use raw C++ arrays as input to the vector constructor, you could do it like this:
   
      vector<some> dest(source, source + 100);
*/

  class TexturesValue : public Value {

    std::vector< std::tr1::shared_ptr<ShaderTexture> > texs_;

  public:
    TexturesValue( const std::tr1::shared_ptr<ShaderTexture> *tex, int size)
  
		: Value( tex[0]->getSamplerType(), size), texs_(tex, tex + size) // a vector of [tex, tex+1,..., tex+count],
	     // tex[0]->getSamplerType(): GL_SAMPLER_2D, GL_SAMPLER_3D, etc.                                                             // where tex+n is a pointer
		 // texs_ ( tex, tex+count ) is a range constructor of vector class ? Yes: the range is defined by a pair of iterators ( sort of pointer)

    {
      assert(size > 0);

      for (int i = 0; i < size; ++i) {

        assert( tex[i]->getSamplerType() == type ); // The type of each texture in the array should be the same as that of the first one

		// The value of type is set by initialization Value( tex[0]->getSamplerType(), count) 

		
		// type is a member of Value class (GL_SMAPLER_2D, GL_SAMPLER_3D, etc)
		//  type: One of the uniform type as returned by glGetActiveUniform, which is executed
		//  just after linking of the shader program.
      }
    }


    virtual Value* clone() const {
      return new TexturesValue(*this);
    }

    virtual void apply(GLint location,  GLsizei count, const GLint *boundTexUnits) const {
       assert(count <= size);
      _helper::genericGlUniformv(location, count, boundTexUnits);
    }

    virtual const std::tr1::shared_ptr<ShaderTexture> *getTextures() const {

      return &texs_[0];

    }
  }; // class TexturesValue

}; // class Uniforms 


#endif