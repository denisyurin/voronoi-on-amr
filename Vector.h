#ifndef CLASS_VECTOR
#define CLASS_VECTOR

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <random>

template <typename T=double> 
class Vector {
	public:
		
		T x, y, z;
		
		Vector(): x(0), y(0), z(0) {}

		Vector(T s): x(s), y(s), z(s) {}

		Vector(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}
		
		
		T & operator [](int i){ return ((T*)this)[i]; }
		
		Vector(T *a){ *this = *(Vector*)a; }
		
		double norm() const {
			return sqrt(x*x + y*y + z*z);
		}
		
		Vector unit(){
			T m = sqrt(x*x + y*y + z*z);
			return Vector(x/m, y/m, z/m);
		}
		
		Vector operator + (Vector V){
			return Vector(x+V.x, y+V.y, z+V.z);
		}
		
		Vector operator / (T s){
			return Vector(x/s, y/s, z/s);
		}
		
		template <typename TT>
		Vector operator - (Vector<TT> V) const {
			return Vector(x-V.x, y-V.y, z-V.z);
		}
		
		Vector operator % (Vector V){ // Cross Product
			return Vector(y*V.z-V.y*z, z*V.x-V.z*x, x*V.y-V.x*y);
		} 
		
		T operator *(Vector V) { // Dot Product
			return x*V.x + y*V.y + z*V.z;
		};  
		
		template <typename TT>
		void operator +=(Vector<TT> V){
			x+=V.x;	y+=V.y;	z+=V.z;
		}	
		
		void operator -=(Vector V){
			x-=V.x;	y-=V.y;	z-=V.z;
		}	
		
		void operator *=(T s){
			x*=s;	y*=s;	z*=s;
		}	
		
		void operator /=(T s){
			x/=s;	y/=s;	z/=s;
		}	

		friend Vector operator *(T s, Vector V){
			return Vector(s*V.x, s*V.y, s*V.z);
		}
		
		friend Vector operator *(Vector V, T s){
			return Vector(s*V.x, s*V.y, s*V.z);
		};
		
		//template <typename T1, typename T2>
		double angle(Vector V){
			T tmp =  (this->unit()*V)/V.norm();
			return tmp<=1 ? acos(tmp) : 0;
		};
		
		operator T*(){
			return (T*)this;
		}
		
		template <typename T2>
		operator Vector<T2>(){
			
			Vector<T2> VT;
			VT.x = x; VT.y = y;VT.z = z;
			
			return VT;
		}
		
		template <typename T1, typename T2>
		Vector rotate(Vector<T1> V, T2 a){ // the direction vector V must be perpendicular to the rotatable vector 'this'!
			return this->norm()*(cos(a)*this->unit()+ sin(a)*V.unit());		
		}
		
		T a1(){
			T tmp = atan2(y,x);
			return tmp<0 ? (tmp + 2*M_PI) : tmp;
		}

		T a2(){
			return asin(z/norm());
		}

		T sqr(){
			return x*x+y*y+z*z;
		}

		Vector setPolar(T a1, T a2, T R){
			*this   = Vector(cos(a1), sin(a1), 0);
			*this*= cos(a2);
			this->z = sin(a2);
			*this*=R; 
			return *this;
		}
		
		Vector rnd(double mod=1){
			
			//std::mt19937 gen(rand()+taskid);
			//std::uniform_real_distribution<> uni_dist(0, 1);
			
			//double cosa1 = 1-2*uni_dist(gen);
			//double a1 = acos(cosa1); 
			//double a2 = 2*M_PI*uni_dist(gen); 
			//double sina1 =	sin(a1);
			
			double cosa1 = 1-2*drand48();
			double a1 = acos(cosa1); 
			double a2 = 2*M_PI*drand48(); 
			double sina1 =	sin(a1);
			
			
			x = mod*sina1*cos(a2);
			y = mod*sina1*sin(a2);
			z = mod*cosa1; 	

			return *this;
		}
		
		void print(const char *str = 0){
			printf("%s{% .15e, % .15e, % .15e}\n", str, x, y, z);
		}
		
		bool check() {
			
			return ( x != x || y != y || z != z /*|| std::isinf(x) || std::isinf(y) || std::isinf(z)*/ );
		}

};

typedef  Vector<double> DVector;
typedef  Vector<float>  FVector; 
typedef  Vector<double> dVector;
typedef  Vector<float>  fVector; 


#endif
