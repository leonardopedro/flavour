#include <stdarg.h>
#include <vector>

namespace std{
	
/// A vector of vectors of vectors of... (N times) of class T objects.
template < class T, int N>
class multivector: public vector< multivector<T, N-1> >{
public:
    typedef vector< multivector< T,N-1> > v;
    
    /// Default constructor
    multivector(): v() {}
	/// Copy constructor
	multivector(const multivector& m): v(m){}
///Recommended constructor
/** Example: multivector<double, 2> m(1.5,4,6), m is a matrix of doubles with dimensions 4x6, with all doubles initialized to 1.5
 * \param value the value with which every objects are initialized
 * \param ... list with the number of dimensions of each vector
 * \see multivector(const T&, va_list &)
 * */
multivector(const T& value, ...){
		va_list listPointer;
		va_start(listPointer,value);
		int n=va_arg(listPointer,int);
		v::insert(v::begin(),n,multivector<T,N-1>(value,listPointer));
		va_end(listPointer);
		}

///Auxiliary constructor (recursive) 
/**\see multivector(const T&, ...)
 * \see multivector<T,1>
 * */
multivector(const T& value, va_list & listPointer)
		{
			int n=va_arg(listPointer,int);
			v::insert(v::begin(),n,multivector<T,N-1>(value,listPointer));
			}

};

/// Specialization template class of \ref multivector<T,N> for N=1
/**\see \ref multivector<T,N>
 * */
template< class T >
class multivector<T,1>: public vector<T>{
public:
 typedef vector< T > v;
 /// Default constructor
 multivector(): v() {}
 /// Copy constructor
	multivector(const multivector& m): v(m){}
 ///Recommended constructor
 /**\param value the value with which every objects are initialized
 * \param x number of dimensions of the vector
 * */
 multivector(const T& value, int x): v(x,value){}
 ///Auxiliary constructor
/**It is the last constructor to be called in the recursive constructor 
 * multivector<T,N>::multivector(const T&,va_list &).
 * \see multivector<T,N>::multivector(const T&,va_list &)
 * */
 multivector(const T& value, va_list & listPointer):
		v(va_arg(listPointer,int), value){}
};

class Matrix: public multivector< ex, 2>{
public:

Matrix(): multivector< ex,2>(0,3,3) {}

Matrix(const Matrix& m): multivector< ex,2>(m) {}
///constructs a symbolic matrix with the symbols names given by the argument
Matrix(const char * m[3][3]): multivector< ex,2>(0,3,3){
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++) at(i)[j]=symbol(m[i][j]);
	}
///constructs a symbolic matrix with the symbols names given by the arguments
Matrix(const char * name,const char ** index1, const char ** index2): multivector< ex,2>(0,3,3){
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++){
			 string res=string(name)+"_{"+string(index1[i])+" "+string(index2[j])+"}";
			 //cout<<res<<endl;
			 at(i)[j]=symbol(res.c_str());
			 }
	}
///constructs a diagonal matrix
Matrix(ex m1, ex m2, ex m3): multivector< ex,2>(0,3,3) {
	at(0)[0]=m1;
	at(1)[1]=m2;
	at(2)[2]=m3;	
}

///constructs a diagonal matrix with all diagonal elements equal
Matrix(ex m1): multivector< ex,2>(0,3,3){
	Matrix();
	at(0)[0]=m1;
	at(1)[1]=m1;
	at(2)[2]=m1;
	}
///constructs a unitary matrix in the standard form
Matrix(ex t12, ex t13, ex t23, ex d13): multivector< ex,2>(0,3,3) {
	Matrix();
	ex c12=cos(t12), c13=cos(t13), c23=cos(t23);
	ex s12=sin(t12), s13=sin(t13), s23=sin(t23);
	ex e13=exp(I*d13);
	ex e13t=ex(1)/e13;
	
	ex aux[3][3]={
			{c12*c13,s12*c13,s13*e13t},
			{-s12*c23-c12*s23*s13*e13,c12*c23-s12*s23*s13*e13,s23*c13},
			{s12*s23-c12*c23*s13*e13,-c12*s23-s12*c23*s13*e13,c13*c23}				
			};
	for(uint i=0;i<3; i++) at(i).assign(aux[i],aux[i]+3);
	}
///used in the unitary constructor
ex cs(ex t12){
	return (exp(I*t12)+1/exp(I*t12))/2;
	}
///used in the unitary constructor
ex sn(ex t12){
	return -I*(exp(I*t12)-1/exp(I*t12))/2;
	}
///computes the hermitian conjugate of the matrix 
Matrix conjugate() const {
	Matrix res;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			res[i][j]=at(j)[i].conjugate();
	
	return res;
}
};

///computes the matrix product
Matrix operator*(const Matrix & m1,const Matrix & m2){
	Matrix res;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
			for(uint k=0;k<3;k++)
				res[i][j]=res[i][j]+m1[i][k]*m2[k][j];
	return res;
}

///computes the matrix sum
Matrix operator+(const Matrix & m1,const Matrix & m2){
	Matrix res;
	for(uint i=0;i<3;i++)
		for(uint j=0;j<3;j++)
				res[i][j]=m1[i][j]+m2[i][j];
	return res;
}



}
