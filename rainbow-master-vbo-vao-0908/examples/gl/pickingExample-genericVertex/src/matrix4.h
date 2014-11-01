#ifndef MATRIX4_H
#define MATRIX4_H

#include <cassert>
#include <cmath>


#include "cvec.h"
#include <string>
#include  <iostream>

/*Thing myThing("asdf");

Instead of
Thing myThing = Thing("asdf");

The first line creates a new object on the stack by calling a constructor of the format Thing(const char*). 

The second one is a bit more complex. It essentially does the following
1.Create an object of type Thing using the constructor Thing(const char*)
2.Create an object of type Thing using the constructor Thing(const Thing&): This is a copy constructor
3.Call ~Thing() on the object created in step #1

==> like Thing myThing( Thing() ); similar to MyThing *myThing = new MyThing("uiae"); But new creates objects on th heap, rather than
on the stack

 The syntax ClassName(constructor params) means to create a temporary object of the type ClassName, 
 using the given parameters for its constructor. So you aren't simply calling the constructor; you are creating an object.
 Constructors are  not functions that you can find and call at your whim.

 When the name of a type is followed by a parenthesized expression list, 
 as is the case here, it is an "Explicity type conversion (functional notation)"

 So impress_errors::Error(err_code) (or impress_errors::Error::Error(err_code) means convert err_code into 
 an impress_errors::Error. Which, in this case, will result in calling the constructor,
 since that's the only way to convert an int into an impress_errors::Error.

 You also sometimes explicitly use a constructor to build a temporary. For example, if you have some class with a constructor:
class Foo
{
    Foo(char* c, int i);
};

and a function
void Bar(Foo foo);

but you don't have a Foo around, you could do
Bar(Foo("hello", 5));

This is like a cast. Indeed, if you have a constructor that takes only one parameter, the C++ compiler will use that constructor to perform implicit casts.

It is not legal to call a constructor on an already-existing object. That is, you cannot do
Foo foo;
foo.Foo();  // compile error!

no matter what you do. But you can invoke a constructor without allocating memory - that's what placement new is for.
char buffer[sizeof(Foo)];      // a bit of memory
Foo* foo = new(buffer) Foo();  // construct a Foo inside buffer

You give new some memory, and it constructs the object in that spot instead of allocating new memory. 
This usage is considered evil, and is rare in most types of code, but common in embedded and data structure code. 

For example, std::vector::push_back uses this technique to invoke the copy constructor. 
That way, it only needs to do one copy, instead of creating an empty object and using the assignment operator.

You're not calling the constructor. The syntax is <typename>(ctor-arg list), and 
to be pedantic that is not the same as <typename>::<typename>(ctor-arg list). The syntax just makes it look like you're calling the constructor. 
In fact, you never "call" a constructor like a function. 
*/


//where the type std::string is defined along with the operator<< for its output (the diagnostic said "there is no operator<< that takes std::string")

//Many compilers include <string> implicitly when including <iostream>, which is why they would compile your program, but the C++ standard only guarantees compilation if #include <string> is present.
using namespace std;

// When an object is created by given its parameters (maybe empty), the default constructor is called.
//For example, MyClass myObject; 

// When an object is created by given another object as its parameter, the copy constructor is called.
// For example, MyClass myObject2(myObject1); 

//When an object is created by any way, and then assign another object to it, the assignment operator is called. 
//	For example, MyClass myObject2; myObject2 = myObject1; 

/*
Matrix m0; // call the default constructor.
	Matrix m1 = m0; // call the copy constructor, not the assignment operator.
	Matrix m2(m0); // call the copy constructor.
	m1 = m0;  // call the assignment operator.


*/
//An assignment operator doesnot invoke the copy constructor. It simply assigns the values of an object to another, member by member. 


/* Functors:
This can be used to create "functors", objects that act like functions:
class Multiplier {
public:
    Multiplier(int m): multiplier(m) {}
    int operator()(int x) { return multiplier * x; }
private:
    int multiplier;
};

Multiplier m(5);
cout << m(4) << endl; => print 5*4
*/


// Forward declaration of Matrix4 and transpose since those are used below
class Matrix4;

//std::ofstream& operator<< (std::ofstream& os, const Matrix4& m);

Matrix4 transpose(const Matrix4& m);

// A 4x4 Matrix.
// To get the element at ith row and jth column, use a(i,j)

class Matrix4 {
  double d_[16]; // layout is row-major, whereas in Open

public:

 Matrix4() { // default   constructor identity matrix construction

 for (int i = 0; i < 16; ++i) {
      d_[i] = 0;
    }
 for (int i = 0; i < 4; ++i) {
      (*this)(i,i) = 1;
    }
  }

 
 Matrix4(const float *a, int num) {
   for (int i = 0; i < num; ++i) {
      d_[i] = a[i];
    }
 }
 // copy constructor and assignment constructor
 //Matrix4 (const Matrix4 &other); we use the default copy constructor
// Matrix4 &operator= (const Matrix4 &other);

 Matrix4(const double a) {
    for (int i = 0; i < 16; ++i) {
      d_[i] = a;
    }
  }


  Matrix4(float t00, float t01, float t02, float t03,
              float t10, float t11, float t12, float t13,
              float t20, float t21, float t22, float t23,
              float t30, float t31, float t32, float t33){
		d_[0]=t00;
		d_[1]=t01;
		d_[2]=t02;
		d_[3]=t03;
		d_[4]=t10;
		d_[5]=t11;
		d_[6]=t12;
		d_[7]=t13;
		d_[8]=t20;
		d_[9]=t21;
		d_[10]=t22;
		d_[11]=t23;
		d_[12]=t30;
		d_[13]=t31;
		d_[14]=t32;
		d_[15]=t33;
		

   }
  double &operator () (const int row, const int col) {
    return d_[(row << 2) + col];
  }

  const double &operator () (const int row, const int col) const {
    return d_[(row << 2) + col];
  }

  double& operator [] (const int i) {
    return d_[i];
  }

  const double& operator [] (const int i) const {
    return d_[i];
  }

 

  template <class T>
  Matrix4& readFromColumnMajorMatrix(const T m[]) {
    for (int i = 0; i < 16; ++i) {
      d_[i] = m[i];
    }
    return *this = transpose(*this);
  }

  template <class T>
  void writeToColumnMajorMatrix( T m[]) const {
    Matrix4 t = transpose(*this);
    for (int i = 0; i < 16; ++i) {
      m[i] = T( t.d_[i] );
    }
  }

  Matrix4 & operator = (const Matrix4& m) { 
	for (int i = 0; i < 16; ++i) {
      d_[i] = m.d_[i];
    }
    return *this;   
   
  }
  Matrix4& operator += (const Matrix4& m) {
    for (int i = 0; i < 16; ++i) {
      d_[i] += m.d_[i];
    }
    return *this;
  }

  Matrix4& operator -= (const Matrix4& m) {
    for (int i = 0; i < 16; ++i) {
      d_[i] -= m.d_[i];
    }
    return *this;
  }

  Matrix4& operator *= (const double a) {
    for (int i = 0; i < 16; ++i) {
      d_[i] *= a;
    }
    return *this;
  }

  Matrix4& operator *= (const Matrix4& a) {
    return *this = *this * a;
  }

  Matrix4 operator + (const Matrix4& a) const {
    return Matrix4(*this) += a;
  }

  Matrix4 operator - (const Matrix4& a) const {
    return Matrix4(*this) -= a;
  }

  Matrix4 operator * (const double a) const {
    return Matrix4(*this) *= a;
  }

Cvec4 operator * (const Cvec4& v) const {
    Cvec4 r(0); // r = 0
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        r[i] += (*this)(i,j) * v(j);
      }
    }
    return r;
  }

Cvec4f operator * (const Cvec4f& v) const {
    Cvec4f r(0); // r = 0
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        r[i] += (*this)(i,j) * v(j);
      }
    }
    return r;
  }

  Matrix4 operator * (const Matrix4& m) const {
    Matrix4 r(0); // r = 0
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
          r(i,k) += (*this)(i,j) * m(j,k);
        }
      }
    }
    return r;
  }

	

  static Matrix4 makeXRotation(const double ang) {
    return makeXRotation(std::cos(ang * CS175_PI/180), std::sin(ang * CS175_PI/180));
  }

  static Matrix4 makeYRotation(const double ang) {
    return makeYRotation(std::cos(ang * CS175_PI/180), std::sin(ang * CS175_PI/180));
  }

  static Matrix4 makeZRotation(const double ang) {
    return makeZRotation(std::cos(ang * CS175_PI/180), std::sin(ang * CS175_PI/180));
  }

  static Matrix4 makeXRotation(const double c, const double s) {
    Matrix4 r; // id
    r(1,1) = r(2,2) = c;
    r(1,2) = -s;
    r(2,1) = s;
    return r;
  }

  static Matrix4 makeYRotation(const double c, const double s) {
    Matrix4 r;
    r(0,0) = r(2,2) = c;
    r(0,2) = s;
    r(2,0) = -s;
    return r;
  }

  static Matrix4 makeZRotation(const double c, const double s) {
    Matrix4 r;
    r(0,0) = r(1,1) = c;
    r(0,1) = -s;
    r(1,0) = s;
    return r;
  }

   
 static Matrix4 makeAxisRotation( double angle, const Cvec3& axis) { // typedef Cvec <double, 3> Cvec3;
    Cvec3  a = normalize( axis );
    //float s = sinf( angle);
	double s = sin( angle );
	//float c = cosf(angle);
	double c = cos(angle);

    Matrix4 m;

    m(0,0) = a[0] * a[0] + (1.f - a[0] * a[0]) * c;
    m(0,1) = a[0] * a[1] * (1.f - c) - a[2] * s;
    m(0,2) = a[0] * a[2] * (1.f - c) + a[1] * s;
    m(0,3) = 0;

    m(1,0) = a[0] * a[1] * (1.f - c) + a[2] * s;
    m(1,1) = a[1] * a[1] + (1.f - a[1] * a[1]) * c;
    m(1,2) = a[1] * a[2] * (1.f - c) - a[0] * s;
    m(1,3) = 0;

    m(2,0) = a[0] * a[2] * (1.f - c) - a[1] * s;
    m(2,1) = a[1] * a[2] * (1.f - c) + a[0] * s;
    m(2,2) = a[2] * a[2] + (1.f - a[2] * a[2]) * c;
    m(2,3) = 0;

    m(3,0) = 0;
    m(3,1) = 0;
    m(3,2) = 0;
    m(3,3) = 1;

    return m;

 }

  static Matrix4 makeTranslation(const Cvec3& t) {
    Matrix4 r;
    for (int i = 0; i < 3; ++i) {
      r(i,3) = t[i];
    }
    return r;
  }

  static Matrix4 makeScale(const Cvec3& s) {
    Matrix4 r;
    for (int i = 0; i < 3; ++i) {
      r(i,i) = s[i];
    }
    return r;
  }

  static Matrix4 makeProjection(
    const double top, const double bottom,
    const double left, const double right,
    const double nearClip, const double farClip) {
    Matrix4 r(0);
    // 1st row
    if (std::abs(right - left) > CS175_EPS) {
      r(0,0) = -2.0 * nearClip / (right - left);
      r(0,2) = (right+left) / (right - left);
    }
    // 2nd row
    if (std::abs(top - bottom) > CS175_EPS) {
      r(1,1) = -2.0 * nearClip / (top - bottom);
      r(1,2) = (top + bottom) / (top - bottom);
    }
    // 3rd row
    if (std::abs(farClip - nearClip) > CS175_EPS) {
      r(2,2) = (farClip+nearClip) / (farClip - nearClip);
      r(2,3) = -2.0 * farClip * nearClip / (farClip - nearClip);
    }
    r(3,2) = -1.0;
    return r;
  }

  // http://unspecified.wordpress.com/2012/06/21/calculating-the-gluperspective-matrix-and-other-opengl-matrix-maths/

  static Matrix4 makeProjection(const double fovy, const double aspectRatio, const double zNear, const double zFar) {
    Matrix4 r(0);
    const double ang = fovy * 0.5 * CS175_PI/180;
    const double f = std::abs(std::sin(ang)) < CS175_EPS ? 0 : 1/std::tan(ang);
    if (std::abs(aspectRatio) > CS175_EPS)
      r(0,0) = f/aspectRatio;  // 1st row

    r(1,1) = f;    // 2nd row

    if (std::abs(zFar - zNear) > CS175_EPS) { // 3rd row
      r(2,2) = (zFar+zNear) / (zFar - zNear);
      r(2,3) = -2.0 * zFar * zNear / (zFar - zNear);
    }

    r(3,2) = -1.0; // 4th row
    return r;
  }

  
// ordinary function ( not a member) which can access the private members of the class
     
friend std::ostream& operator<< (std::ostream& os, const Matrix4& m) {
	
	
	os <<  "[" <<   m(0,0) <<  setw(15) << m(0,1) << setw(15) << m(0,2)  <<  setw(15) << m(0,3) <<  "]" << endl;
	os <<  "[" <<  m(1,0) <<  setw(15) <<m(1,1) << setw(15) << m(1,2)  << setw(15) << m(1,3) <<  "]" << endl;
	os <<  "[" <<  m(2,0) <<  setw(15) << m(2,1) << setw(15) <<m(2,2)  <<  setw(15) << m(2,3) <<  "]" << endl;
	os <<  "[" <<  m(3,0) <<  setw(15) <<m(3,1) <<  setw(15) <<m(3,2)  <<   setw(15) <<m(3,3) <<  "]" << endl;
	

	return os;
	
} 

}; // class Matrix4


inline bool isAffine(const Matrix4& m) {
  return std::abs(m[15]-1) + std::abs(m[14]) + std::abs(m[13]) + std::abs(m[12]) < CS175_EPS;
}

inline double norm2(const Matrix4& m) {
  double r = 0;
  for (int i = 0; i < 16; ++i) {
    r += m[i]*m[i];
  }
  return r;
}

// computes inverse of affine matrix. assumes last row is [0,0,0,1]
inline Matrix4 inv(const Matrix4& m) {
  Matrix4 r;                                              // default constructor initializes it to identity
  assert(isAffine(m));
  double det = m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) +
               m(0,1)*(m(1,2)*m(2,0) - m(1,0)*m(2,2)) +
               m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));

  // check non-singular matrix
  assert(std::abs(det) > CS175_EPS3);

  // "rotation part"
  r(0,0) =  (m(1,1) * m(2,2) - m(1,2) * m(2,1)) / det;
  r(1,0) = -(m(1,0) * m(2,2) - m(1,2) * m(2,0)) / det;
  r(2,0) =  (m(1,0) * m(2,1) - m(1,1) * m(2,0)) / det;
  r(0,1) = -(m(0,1) * m(2,2) - m(0,2) * m(2,1)) / det;
  r(1,1) =  (m(0,0) * m(2,2) - m(0,2) * m(2,0)) / det;
  r(2,1) = -(m(0,0) * m(2,1) - m(0,1) * m(2,0)) / det;
  r(0,2) =  (m(0,1) * m(1,2) - m(0,2) * m(1,1)) / det;
  r(1,2) = -(m(0,0) * m(1,2) - m(0,2) * m(1,0)) / det;
  r(2,2) =  (m(0,0) * m(1,1) - m(0,1) * m(1,0)) / det;

  // "translation part" - multiply the translation (on the left) by the inverse linear part
  r(0,3) = -(m(0,3) * r(0,0) + m(1,3) * r(0,1) + m(2,3) * r(0,2));
  r(1,3) = -(m(0,3) * r(1,0) + m(1,3) * r(1,1) + m(2,3) * r(1,2));
  r(2,3) = -(m(0,3) * r(2,0) + m(1,3) * r(2,1) + m(2,3) * r(2,2));
  assert(isAffine(r) && norm2(Matrix4() - m*r) < CS175_EPS2);
  return r;
}

inline Matrix4 transpose(const Matrix4& m) {
  Matrix4 r(0);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      r(i,j) = m(j,i);
    }
  }
  return r;
}

inline Matrix4 normalMatrix(const Matrix4& m) {
  Matrix4 invm = inv(m);
  invm(0, 3) = invm(1, 3) = invm(2, 3) = 0;
  return transpose(invm);
}

inline Matrix4 transFact(const Matrix4& m) {
  // TODO
}

inline Matrix4 linFact(const Matrix4& m) {
  // TODO
}

// Matrix3



// Forward declaration of Matrix4 and transpose since those are used below
class Matrix3;

//std::ofstream& operator<< (std::ofstream& os, const Matrix4& m);

Matrix3 transpose(const Matrix3& m);

// A 3x3 Matrix.
// To get the element at ith row and jth column, use a(i,j)

class Matrix3 {
  double d_[9]; // layout is row-major, whereas in Open

public:

 Matrix3() { // default   constructor identity matrix construction

 for (int i = 0; i < 9; ++i) {
      d_[i] = 0;
    }
 for (int i = 0; i < 3; ++i) {
      (*this)(i,i) = 1;
    }
  }

 
 Matrix3(const float *a, int num) {
   for (int i = 0; i < num; ++i) {
      d_[i] = a[i];
    }
 }
 // copy constructor and assignment constructor
 //Matrix4 (const Matrix4 &other); we use the default copy constructor
// Matrix4 &operator= (const Matrix4 &other);

 Matrix3(const double a) {
    for (int i = 0; i < 9; ++i) {
      d_[i] = a;
    }
  }


  Matrix3(float t00, float t01, float t02, 
              float t10, float t11, float t12,
              float t20, float t21, float t22
              ){
		d_[0]=t00;
		d_[1]=t01;
		d_[2]=t02;
		
		d_[3]=t10;
		d_[4]=t11;
		d_[5]=t12;
		
		d_[6]=t20;
		d_[7]=t21;
		d_[8]=t22;
		
		

   }

  Matrix3( Cvec3 c0, Cvec3 c1, Cvec3 c2 ){
	  // set columun vector c_i to column j_i

		d_[0]=c0[0];
		d_[1]=c1[0];
		d_[2]=c2[0];
		
		d_[3]=c0[1];
		d_[4]=c1[1];
		d_[5]=c2[1];
		
		d_[6]=c0[2];
		d_[7]=c1[2];
		d_[8]=c2[2];
		
		

   }

  double &operator () (const int row, const int col) {
    return d_[ row * 3  + col];
  }

  const double &operator () (const int row, const int col) const {
    return d_[ row * 3 + col];
  }

  double& operator [] (const int i) {
    return d_[i];
  }

  const double& operator [] (const int i) const {
    return d_[i];
  }

 

  template <class T>
  Matrix3& readFromColumnMajorMatrix(const T m[]) {
    for (int i = 0; i < 9; ++i) {
      d_[i] = m[i];
    }
    return *this = transpose(*this);
  }

  template <class T>
  void writeToColumnMajorMatrix( T m[]) const {
    Matrix3 t = transpose(*this);
    for (int i = 0; i < 9; ++i) {
      m[i] = T( t.d_[i] );
    }
  }

  Matrix3 & operator = (const Matrix3& m) { 
	for (int i = 0; i < 9; ++i) {
      d_[i] = m.d_[i];
    }
    return *this;   
   
  }
  Matrix3& operator += (const Matrix3& m) {
    for (int i = 0; i < 9; ++i) {
      d_[i] += m.d_[i];
    }
    return *this;
  }

  Matrix3& operator -= (const Matrix3& m) {
    for (int i = 0; i < 9; ++i) {
      d_[i] -= m.d_[i];
    }
    return *this;
  }

  Matrix3& operator *= (const double a) {
    for (int i = 0; i <9; ++i) {
      d_[i] *= a;
    }
    return *this;
  }

  Matrix3& operator *= (const Matrix3& a) {
    return *this = *this * a;
  }

  Matrix3 operator + (const Matrix3& a) const {
    return Matrix3(*this) += a;
  }

  Matrix3 operator - (const Matrix3& a) const {
    return Matrix3(*this) -= a;
  }

  Matrix3 operator * (const double a) const {
    return Matrix3(*this) *= a;
  }

Cvec3 operator * (const Cvec3& v) const {
    Cvec3 r(0); // r = 0
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        r[i] += (*this)(i,j) * v(j);
      }
    }
    return r;
  }

Cvec3f operator * (const Cvec3f& v) const {
    Cvec3f r(0); // r = 0
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        r[i] += (*this)(i,j) * v(j);
      }
    }
    return r;
  }

  Matrix3 operator * (const Matrix3& m) const {
    Matrix3 r(0); // r = 0
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          r(i,k) += (*this)(i,j) * m(j,k);
        }
      }
    }
    return r;
  }

   static Matrix3 makeScale(const Cvec3& s) {
    Matrix3 r;
    for (int i = 0; i < 3; ++i) {
      r(i,i) = s[i];
    }
    return r;
  }

  
// ordinary function ( not a member) which can access the private members of the class
     
friend std::ostream& operator<< (std::ostream& os, const Matrix3& m) {
	
	
	os <<  "[" <<   m(0,0) <<  setw(15) << m(0,1) << setw(15) << m(0,2)  <<  setw(15)  <<  "]" << endl;
	os <<  "[" <<  m(1,0) <<  setw(15) <<m(1,1) << setw(15) << m(1,2)  << setw(15) <<  "]" << endl;
	os <<  "[" <<  m(2,0) <<  setw(15) << m(2,1) << setw(15) <<m(2,2)  <<  setw(15) << "]" << endl;
	

	return os;
	
} 

}; // class Matrix3


inline double norm2(const Matrix3& m) {
  double r = 0;
  for (int i = 0; i < 9; ++i) {
    r += m[i]*m[i];
  }
  return r;
}


//http://www.ryanjuckett.com/programming/rgb-color-space-conversion/

// computes inverse of  matrix3. 
inline  Matrix3 inv( const Matrix3& m) {
  Matrix3 r;                                              // default constructor initializes it to identity
  
   // calculate the minors for the first row
  double  minor00 = m(1,1)*m(2,2) - m(1,2)*m(2,1);
  double  minor01 = m(1,2)*m(2,0) - m(1,0)*m(2,2);
  double  minor02 = m(1,0)*m(2,1) - m(1,1)*m(2,0);
  
    // calculate the determinant
  double det =   m(0,0) * minor00
                        + m(0,1) * minor01
                        + m(0,2) * minor02;
  
  // check non-singular matrix
  assert(std::abs(det) > CS175_EPS3);

  // the inverse of inMat is (1 / determinant) * adjoint(inMat)
    double  invDet = 1.0f / det;
    r(0,0) = invDet * minor00;
    r(0,1) = invDet * (m(2,1)* m(0,2) - m(2,2)*m(0,1) );
    r(0,2) = invDet * ( m(0,1) * m(1,2) - m(0,2)* m(1,1) );
  
    r(1,0) = invDet * minor01;
    r(1,1)  = invDet * ( m(2,2)* m(0,0) - m(2,0)* m(0,2) );
    r(1,2) = invDet * ( m(0,2) * m(1,0) - m(0,0) * m(1,2) );
  
    r(2,0)  = invDet * minor02;
    r(2,1)  = invDet * ( m(2,0) * m(0,1)  - m(2,1) * m(0,0) );
    r(2,2) = invDet * ( m(0,0) * m(1,1)  - m(0,1) * m(1,0) );

  
  return r;
}

inline Matrix3 transpose(const Matrix3& m) {
  Matrix3 r(0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r(i,j) = m(j,i);
    }
  }
  return r;
}


// GAMMA Correction: 
// the electron guns of a CRT monitor have a non-linear output and must be supplied with gamma corrected 
// input such that the appropriate real-world luminance is created.

 

//******************************************************************************
// Convert an sRGB color channel to a linear sRGB color channel.
//******************************************************************************
inline float gammaExpandsRGB(float nonlinear) // expand = decoding
{    
    return   ( nonlinear <= 0.04045f )
           ? ( nonlinear / 12.92f )
           : ( pow( (nonlinear+0.055f)/1.055f, 2.4f ) );
}
  
//******************************************************************************
// Convert a linear sRGB color channel to a sRGB color channel.
//******************************************************************************
inline float gammaCompresssRGB(float linear) // gamma compress = gamma encoding = gamma correction = linear^{1/g), where g > 1.
{    
    return   ( linear <= 0.0031308f )
           ? ( 12.92f * linear )
           : ( 1.055f * pow( linear, 1.0f/2.4f ) - 0.055f );
}
  
//******************************************************************************
// Convert an sRGB color to a linear sRGB color.
//******************************************************************************
inline Cvec3  gammaExpandsRGB(Cvec3  pColor)
{   Cvec3 c;
   c [0]  = gammaExpandsRGB( pColor  [0] );
   c [1]  = gammaExpandsRGB( pColor[1] );
   c [2]  = gammaExpandsRGB( pColor [2] );
   return c;

}
  
//******************************************************************************
// Convert a linear sRGB color to an sRGB color.
//******************************************************************************
inline Cvec3 gammaCompresssRGB(Cvec3 pColor)
{
	Cvec3 c;
    c[0]  = gammaCompresssRGB( pColor[0]);
    c[1]  = gammaCompresssRGB( pColor[1] );
    c[2]  = gammaCompresssRGB( pColor[2] );
}


//



//******************************************************************************
// Convert a linear sRGB color to an sRGB color 
//******************************************************************************
inline Matrix3  makeConversionMatrixRGBtoXYZ
(
    
    const Cvec2  & red_xy,   // xy chromaticity coordinates of the red primary
    const Cvec2  & green_xy, // xy chromaticity coordinates of the green primary
    const Cvec2  & blue_xy,  // xy chromaticity coordinates of the blue primary
    const Cvec2 & white_xy  // xy chromaticity coordinates of the white point
)
{
    // generate xyz chromaticity coordinates (x + y + z = 1) from xy coordinates
    Cvec3 r ( red_xy[0],   red_xy[1],   1.0f - (red_xy[0] + red_xy[1]) );
    Cvec3 g ( green_xy[0], green_xy[1], 1.0f - (green_xy[0] + green_xy[1]) );
    Cvec3 b  ( blue_xy[0],  blue_xy[1],  1.0f - (blue_xy[0] + blue_xy[1]) );
    Cvec3 w ( white_xy[0], white_xy[1], 1.0f - (white_xy[0] + white_xy[1]) );
  
    // Convert white xyz coordinate to XYZ coordinate by letting that the white
    // point has  XYZ relative luminance of 1.0. Relative luminance is the Y
    // component of and XYZ color.
    //   XYZ = xyz * (Y / y) => XYZ = xyz / y = (x/y, 1, z/y)
    w[0]  /= white_xy[1];
    w[1]  /= white_xy[1];
    w[2]  /= white_xy[1];
  
    // Solve for the transformation matrix 'M' from RGB to XYZ
    // * We know that the columns of M are equal to the unknown XYZ values of r, g and b.
    // * We know that the XYZ values of r, g and b are each a scaled version of the known
    //   corresponding xyz chromaticity values.
    // * We know the XYZ value of white based on its xyz value and the assigned relative
    //   luminance of 1.0.
    // * We know the RGB value of white is (1,1,1).
    //                  
    //   white_XYZ = M * white_RGB
    //
    //       [r.x g.x b.x]
    //   N = [r.y g.y b.y]
    //       [r.z g.z b.z]
    //
    //       [sR 0  0 ]
    //   S = [0  sG 0 ]
    //       [0  0  sB]
    //
    //   M = N * S
    //   white_XYZ = N * S * white_RGB
    //   N^-1 * white_XYZ = S * white_RGB = (sR,sG,sB)
    //
    // We now have an equation for the components of the scale matrix 'S' and
    // can compute 'M' from 'N' and 'S'
  
    Matrix3 N (r, g, b );
  
    Matrix3  invN = inv( N );
	  
    Cvec3 scale = invN *  w;
  
	Matrix3 S = Matrix3::makeScale( scale );
    Matrix3 M = N * S;
	
	return M;
} // makeConversionMatrixRGBtoXYZ
 

//******************************************************************************
// Example of using the color space conversion functions
//******************************************************************************
inline Matrix3 makeConversionMatrixsRGBtoXYZ()
{
    // define chromaticity coordinates for sRGB space
    Cvec2  sRGB_red_xy (0.64f, 0.33f );
    Cvec2  sRGB_green_xy ( 0.30f, 0.60f );
    Cvec2  sRGB_blue_xy  ( 0.15f, 0.06f );
    Cvec2  sRGB_white_xy ( 0.3127f, 0.3290f );
  
    // generate conversion matrix from linear sRGB space to XYZ space
    Matrix3 sRGBtoXYZ = makeConversionMatrixRGBtoXYZ(
                                         sRGB_red_xy,
                                         sRGB_green_xy,
                                         sRGB_blue_xy,
                                         sRGB_white_xy );
    return sRGBtoXYZ;

	/*
    // generate conversion matrix from XYZ space to linear sRGB space
    Matrix3 g_XYZ_to_sRGB = inv( sRGB_to_XYZ );
  
    // define a color in sRGB space
    Cvec3 myColor ( 0.2f, 0.5f, 0.8f );
  
    // convert form sRGB to XYZ
    // convert from gamma-corrected sRGB to linear sRGB
    Cvec3 linsRGB =gammaExpandsRGB( myColor );
  
        // convert from linear sRGB to XYZ
    Cvec3 XYZ  = sRGB_to_XYZ * linsRGB;
      
    // convert form XYZ back to sRGB
            // convert from XYZ to linear sRGB
    linsRGB = g_XYZ_to_sRGB * XYZ;
  
        // convert from linear sRGB to gamma-corrected sRGB
     Cvec3 sRGB = gammaCompresssRGB( linsRGB  );
    */

} // makeConversionMatrixRGBtoXYZ

inline Matrix3 makeConversionMatrixXYZtosRGB()
{
    // define chromaticity coordinates for sRGB space
    Cvec2  sRGB_red_xy  ( 0.64f, 0.33f );
    Cvec2  sRGB_green_xy ( 0.30f, 0.60f );
    Cvec2  sRGB_blue_xy  ( 0.15f, 0.06f );
    Cvec2  sRGB_white_xy ( 0.3127f, 0.3290f );
  
    // generate conversion matrix from linear sRGB space to XYZ space
    Matrix3 sRGBtoXYZ = makeConversionMatrixRGBtoXYZ(
                                         sRGB_red_xy,
                                         sRGB_green_xy,
                                         sRGB_blue_xy,
                                         sRGB_white_xy );
    

	
    // generate conversion matrix from XYZ space to linear sRGB space
    Matrix3 XYZtosRGB = inv( sRGBtoXYZ );
    
	return XYZtosRGB;

	/*
    // define a color in sRGB space
    Cvec3 myColor ( 0.2f, 0.5f, 0.8f );
  
    // convert form sRGB to XYZ
    // convert from gamma-corrected sRGB to linear sRGB
    Cvec3 linsRGB =gammaExpandsRGB( myColor );
  
        // convert from linear sRGB to XYZ
    Cvec3 XYZ  = sRGB_to_XYZ * linsRGB;
      
    // convert form XYZ back to sRGB
            // convert from XYZ to linear sRGB
    linsRGB = g_XYZ_to_sRGB * XYZ;
  
        // convert from linear sRGB to gamma-corrected sRGB
     Cvec3 sRGB = gammaCompresssRGB( linsRGB  );
    */

} // makeConversionMatrixXYZtoRGB
 

 

#endif

