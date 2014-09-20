#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <cassert>
#include <string>  // to use: operator<<(std::ostream&, std::string)
#include <stdio.h>
#include <iostream>
#include <iomanip> // for std::setw

#include <fstream>
#include <istream>


using namespace std;

static const double CS175_PI = 3.14159265358979323846264338327950288;
static const double CS175_EPS = 1e-8;
static const double CS175_EPS2 = CS175_EPS * CS175_EPS;
static const double CS175_EPS3 = CS175_EPS * CS175_EPS * CS175_EPS;

//template<typename T, int n> class Cvec;  // pre-declare the template class itself
  
template <typename T, int n>
class Cvec {

  T d_[n];

public:

  Cvec() {
    for (int i = 0; i < n; ++i) {
      d_[i] = 0;
    }
  }

  Cvec(const T& t) {
    for (int i = 0; i < n; ++i) {
      d_[i] = t;
    }
  }

   Cvec(const T * t, int m) {
    assert ( m <= n);
    for (int i = 0; i < m; ++i) {  // m <= n should be the case
      d_[i] = t[i];
    }
  }



  Cvec(const T& t0, const T& t1) {
    assert(n == 2); // better to use static_assert from c++11
    d_[0] = t0, d_[1] = t1;
  }

  Cvec(const T& t0, const T& t1, const T& t2) {
    assert(n == 3); // better to use static_assert from c++11
    d_[0] = t0, d_[1] = t1, d_[2] = t2;
  }

  Cvec(const T& t0, const T& t1, const T& t2, const T& t3) {
    assert(n == 4); // better to use static_assert from c++11
    d_[0] = t0, d_[1] = t1, d_[2] = t2, d_[3] = t3;
  }

  // either truncate if m < n, or extend with extendValue
  template<int m>
  explicit Cvec(const Cvec<T, m>& v, const T& extendValue = T(0)) {
    for (int i = 0; i < std::min(m, n); ++i) {
      d_[i] = v[i];
    }
    for (int i = std::min(m, n); i < n; ++i) {  // n < n is false.
      d_[i] = extendValue;
    }
  }

  T& operator [] (const int i) {
    return d_[i];
  }

  const T& operator [] (const int i) const {
    return d_[i];
  }

  T& operator () (const int i) {
    return d_[i];
  }

  const T& operator () (const int i) const {
    return d_[i];
  }

  Cvec operator - () const {
    return Cvec(*this) *= -1;
  }

  Cvec& operator += (const Cvec& v) {
    for (int i = 0; i < n; ++i) {
      d_[i] += v[i];
    }
    return *this;
  }

 
  Cvec& operator -= (const Cvec& v) {
    for (int i = 0; i < n; ++i) {
      d_[i] -= v[i];
    }
    return *this;
  }

  Cvec& operator *= (const T a) {
    for (int i = 0; i < n; ++i) {
      d_[i] *= a;
    }
    return *this;
  }


  Cvec&  operator *= (const Cvec  &v) { 
	for (int i = 0; i < n; ++i) {
      d_[i] *= v[i];
    }
	return *this;
  }

  
  
  Cvec& operator /= (const T a) {
    const T inva(1/a);
    for (int i = 0; i < n; ++i) {
      d_[i] *= inva;
    }
    return *this;
  }

  Cvec operator + (const Cvec& v) const {
    return Cvec(*this) += v;  // Cvec(*this) sets d[] array
  }

  Cvec operator - (const Cvec& v) const {
    return Cvec(*this) -= v;
  }

  Cvec operator * (const T a) const {
    return Cvec(*this) *= a;
  }

  Cvec operator / (const T a) const {
    return Cvec(*this) /= a;
  }

  // Normalize self and returns self
  Cvec& normalize() {
    assert(dot(*this, *this) > CS175_EPS2);
    return *this /= std::sqrt(dot(*this, *this));
  }
  
  bool operator == (const Cvec& v) {
    for (int i = 0; i < n; ++i) {
      if ( (*this)[i] != v[i] ) return false;
    }
    return true;
  }

  // Overload an ordinary operator  << ( not a member) as a friend of Cvec so that it 
   // can access the private members of class Cvec
  //    friend std::ostream& operator<< <> (std::ostream& o, const Foo<T>& x);
   friend std::ostream& operator<< ( std::ostream& os, const Cvec &v) {
	if ( n == 3) {
	     os << "(" << v[0] << "," << v[1] << "," << v[2] << ") " << endl;
	}
	else if (n == 4) {
         os << "(" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << ") " << endl;
	}

	return os;	
   }
	 
};

/* 
Inline member functions (C++ only)
You may either define a member function inside its class definition, or you may define it outside if you have already declared (but not defined) the member function in the class definition.

A member function that is defined inside its class member list is called an inline member function. Member functions containing a few lines of code are usually declared inline. In the above example, add() is an inline member function. If you define a member function outside of its class definition, it must appear in a namespace scope enclosing the class definition. You must also qualify the member function name using the scope resolution (::) operator.

An equivalent way to declare an inline member function is to either declare it in the class with the inline keyword (and define the function outside of its class) or to define it outside of the class declaration using the inline keyword.

In the following example, member function Y::f() is an inline member function:

struct Y {
private:
  char* a;
public:
  char* f() { return a; }
};The following example is equivalent to the previous example; Y::f() is an inline member function: 

struct Y {
private:
  char* a;
public:
  char* f();
};

inline char* Y::f() { return a; }The inline specifier does not affect the linkage of a member or nonmember function: linkage is external by default.

Member functions of a local class must be defined within their class definition. As a result, member functions of a local class are implicitly inline functions. These inline member functions have no linkage.



*/
template<typename T>
inline Cvec<T,3> cross(const Cvec<T,3>& a, const Cvec<T,3>& b) {
  return Cvec<T,3>(a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0));
}

// using inline allows you to put that function in a header file, 
// where it can be included in multiple source files.
// Using inline makes the identifier in file scope, much like declaring it static. 
// Without using inline, you would get a multiple symbol definition error from the linker.



template<typename T, int n>
inline T dot(const Cvec<T,n>& a, const Cvec<T,n>& b) {
  T r(0);
  for (int i = 0; i < n; ++i) {
    r += a(i)*b(i);
  }
  return r;
}

template<typename T, int n>
inline T norm2(const Cvec<T, n>& v) {
  return dot(v, v);
}

template<typename T, int n>
inline T norm(const Cvec<T, n>& v) {
  return std::sqrt(dot(v, v));
}

// Return a normalized vector without modifying the input (unlike the member
// function version v.normalize() ).
template<typename T, int n>
inline Cvec<T, n> normalize(const Cvec<T,n>& v) {
  assert(dot(v, v) > CS175_EPS2);
  return v / norm(v);
}

// element of type double precision float
typedef Cvec <double, 2> Cvec2;
typedef Cvec <double, 3> Cvec3;
typedef Cvec <double, 4> Cvec4;

// element of type single precision float
typedef Cvec <float, 2> Cvec2f;
typedef Cvec <float, 3> Cvec3f;
typedef Cvec <float, 4> Cvec4f;

// elements of type unsigned byte
typedef Cvec <unsigned char, 2> Cvec2ub;
typedef Cvec <unsigned char, 3> Cvec3ub;
typedef Cvec <unsigned char, 4> Cvec4ub;


 
#endif
