#ifndef VECTOR3_H        // Include guard (prevents multiple includes)
#define VECTOR3_H

#include <array>

///A general re-usable template structure for 3-vectors containing doubles
//Typical to declare AND define template structure in a single header file
//Will include all the standard vector operations

template <typename Derived>
struct Vector3Template{
    //Accessing elements by reference in constant or non-constant fashion
    double& x() {return static_cast<Derived*>(this)->x();};
    const double& x() const {return static_cast<const Derived*>(this)->x();};
    double& y() {return static_cast<Derived*>(this)->y();};
    const double& y() const {return static_cast<const Derived*>(this)->y();};
    double& z() {return static_cast<Derived*>(this)->z();};
    const double& z() const {return static_cast<const Derived*>(this)->z();};

    //In place addition (vector and scalar)
    template<typename OtherDerived>
    Derived& operator+=(const Vector3Template<OtherDerived>& other){
        x()+=other.x();
        y()+=other.y();
        z()+=other.z();
        return static_cast<Derived&>(*this);
    }
    Derived& operator+=(const double scalar){
        x()+=scalar;
        y()+=scalar;
        z()+=scalar;
        return static_cast<Derived&>(*this);
    }
    
    //In place subtraction (vector and scalar)
    template<typename OtherDerived>
    Derived& operator-=(const Vector3Template<OtherDerived>& other){
        x()-=other.x();
        y()-=other.y();
        z()-=other.z();
        return static_cast<Derived&>(*this);
    }
    Derived& operator-=(const double scalar){
        x()-=scalar;
        y()-=scalar;
        z()-=scalar;
        return static_cast<Derived&>(*this);
    }

    //In place scalar multiplication
    Derived& operator*=(const double scalar){
        x()*=scalar;
        y()*=scalar;
        z()*=scalar;
        return static_cast<Derived&>(*this);
    }

    //In place scalar division
    Derived& operator/=(const double scalar){
        x()/=scalar;
        y()/=scalar;
        z()/=scalar;
        return static_cast<Derived&>(*this);
    }

    //Get the magnitude squared
    double magsquared() const{
        return dot(static_cast<Derived&>(*this),static_cast<Derived&>(*this));
    }
};


struct Vector3 : public Vector3Template<Vector3>{
    //A copy of the source data, completely self-owned mutable etc.
    std::array<double,3> data;
  
    //Constructor from parameter values
    Vector3(double x=0, double y=0, double z=0);
    //Constructor by conversion
    template<typename Derived>
    explicit Vector3(const Vector3Template<Derived>& other);
    
    //Access elements by name
    double& x();
    double& y();
    double& z();

    const double& x() const;
    const double& y() const;
    const double& z() const; 

};

struct Vector3View : public Vector3Template<Vector3View>{
    //A simple read and write version of the full Vector3D that avoids copying
    //It returns a mutable view of the x,y,z components 
    std::array<double*,3> data;
    
    //Constructor
    Vector3View(double* x= nullptr, double* y= nullptr, double* z= nullptr);

    //Access elements by name
    double& x(); //non-constant member functions allow read and write
    double& y();
    double& z();

    const double& x() const; //constant member functions for read only
    const double& y() const;
    const double& z() const;

};

struct Vector3ConstView : public Vector3Template<Vector3ConstView>{
    //A read only immutable view of the data
    //Methods such as += will obviously not be valid but + and - still work ok
    std::array<const double*,3> data;
    
    //Constructor
    Vector3ConstView(const double* x=nullptr, const double* y=nullptr, const double* z= nullptr);


    const double& x() const; //only constant member functions
    const double& y() const;
    const double& z() const;

};


//NON MEMBER FUNCTIONS (always want to return a real vector object not a view)


//Addition and subtraction for vectors and scalars (creates a new Vector) uses NumPy broadcasting rules
template<typename DerivedL, typename DerivedR>
Vector3 operator+(const Vector3Template<DerivedL>& lhs, const Vector3Template<DerivedR>& rhs)  {
        Vector3 result(lhs);//Enforce its a real vector type
        result+=rhs;//Makes a copy of the lhs object 
        return result;
    }
template<typename Derived>
Vector3 operator+(const double scalar, const Vector3Template<Derived>& Vec)  {
        Vector3 result(Vec);
        result.x()+=scalar;
        result.y()+=scalar;
        result.z()+=scalar;
        return result;
    }
template<typename Derived>
Vector3 operator+(const Vector3Template<Derived>& Vec,const double scalar)  {
        return scalar+Vec;
    }

template<typename DerivedL, typename DerivedR>
Vector3 operator-(const Vector3Template<DerivedL>& lhs, const Vector3Template<DerivedR>& rhs)  {
        Vector3 result(lhs);
        result-=rhs;
        return result;
    }
template<typename Derived>
Vector3 operator-(const double scalar,const Vector3Template<Derived>& rhs)  {
        Vector3 result(rhs);
        result.x()-=scalar;
        result.y()-=scalar;
        result.z()-=scalar;
        return result;
    }
template<typename Derived>
Vector3 operator-(const Vector3Template<Derived>& lhs,const double scalar)  {
        return scalar-lhs;
    }
 
 
//Scalar multiplication and division
template<typename Derived>
Vector3 operator*(const Vector3Template<Derived>& lhs,const double scalar)  {
        Vector3 result(lhs);
        result*=scalar;
        return result;
    }
template<typename Derived>
Vector3 operator*(const double scalar,const Vector3Template<Derived>& lhs)  {
        return lhs*scalar;
    }
template<typename Derived>
Vector3 operator/(const Vector3Template<Derived>& lhs,const double scalar)  {
        Vector3 result(lhs);
        result/=scalar;
        return result;
    }
template<typename Derived>
Vector3 operator/(const double scalar,Vector3Template<Derived> lhs)  {
        return lhs/scalar;
    }


//Inner (dot) product
template<typename DerivedL, typename DerivedR>
double dot(const Vector3Template<DerivedL>& a, const Vector3Template<DerivedR>& b){
    return a.x()*b.x()+ a.y()*b.y() + a.z()*b.z();
 }


//Cross product
template<typename DerivedL, typename DerivedR>
Vector3 cross(const Vector3Template<DerivedL>& a, const Vector3Template<DerivedR>& b){
    return Vector3((a.y()*b.z()-a.z()*b.y()),(a.z()*b.x()-a.x()*b.z()),(a.x()*b.y()-a.y()*b.x()));
 }


#endif