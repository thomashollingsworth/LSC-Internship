#ifndef STATEVECTOR_H
#define STATEVECTOR_H

#include <array>
#include "Vector3.h"


template <typename Derived>
struct StateVectorTemplate{
    //This will be a template for a general vector of doubles with 9 components
    //The template will describe accessing by index but not specify any names
    //The template will describe general vector addition, scalar multiplication etc.
    //ConservedStateVector and PrimitiveStateVector will be a special case of this template
    //StateVector will be just a basic copy of this that can be used in main

    std::array<double,9> data; //This will be a shared member for all child structures
    
    //Construction by individual parameters as doubles and construction by conversion
    StateVectorTemplate(
        double val0=0,
        double val1=0,
        double val2=0,
        double val3=0,
        double val4=0,
        double val5=0,
        double val6=0,
        double val7=0,
        double val8=0
    ): data{{val0,val1,val2,val3,val4,val5,val6,val7,val8}} {}
    

    //Allows for explicit conversion between different types of StateVectos
    template <typename Other>
    explicit StateVectorTemplate(const StateVectorTemplate<Other>& other)
    : data(other.data) {}

    // Access by index is a shared member function 
    double& operator[](std::size_t i) { return data[i]; }
    const double& operator[](std::size_t i) const { return data[i]; }

    //In place addition (vector and scalar)
    template<typename OtherDerived>
    Derived& operator+=(const Vector3Template<OtherDerived>& other){
        for(size_t i=0;i<9;i++){
        data[i]+=other[i];
    }
        return static_cast<Derived&>(*this);//Return values enable chaining of operators!
    }
    Derived& operator+=(const double scalar){
        for(size_t i=0;i<9;i++){
        data[i]+=scalar;
    }
        return static_cast<Derived&>(*this);
    }
    
    //In place subtraction (vector and scalar)
    template<typename OtherDerived>
    Derived& operator-=(const Vector3Template<OtherDerived>& other){
        for(size_t i=0;i<9;i++){
        data[i]-=other[i];
    }
        return static_cast<Derived&>(*this);
    }
    Derived& operator-=(const double scalar){
        for(size_t i=0;i<9;i++){
        data[i]-=scalar;
    }
        return static_cast<Derived&>(*this);
    }

    //In place scalar multiplication
    Derived& operator*=(const double scalar){
        for(size_t i=0;i<9;i++){
        data[i]*=scalar;
    }
        return static_cast<Derived&>(*this);
    }

    //In place scalar division
    Derived& operator/=(const double scalar){
        for(size_t i=0;i<9;i++){
        data[i]/=scalar;
    }
        return static_cast<Derived&>(*this);
    }

};

//Forward declarations
struct ConservedStateVector;
struct StateVector;
struct PrimitiveStateVector;

struct StateVector : public StateVectorTemplate<StateVector>{
    //Will just inherit all feautures of the template but can be used
    using StateVectorTemplate::StateVectorTemplate; //(inherited constructors)
};


struct PrimitiveStateVector : public StateVectorTemplate<PrimitiveStateVector>{
    
    using StateVectorTemplate::StateVectorTemplate; //(inherited constructors)
    //Access scalar parameters by name
    double& density();   
    double& pressure();    
    double& psi();  
    const double& density() const;   
    const double& pressure() const;    
    const double& psi() const;      
    //Access vector3 parameters by name
    //B and momentum are accessed using Vector3View structure
    //Vector3View function allows direct read and write access
    //Vector3ConstView function is read only and can be included in other const member functions
    Vector3View B();
    Vector3View velocity();
    const Vector3ConstView B() const;
    const Vector3ConstView velocity() const;

    ConservedStateVector prim_to_con(const double gamma) const;
    //TOTAL PRESSURE p_T 
    double pressure_T();
    double pressure_T() const;
  

    
};

struct ConservedStateVector : public StateVectorTemplate<ConservedStateVector>{
    
    using StateVectorTemplate::StateVectorTemplate; //(inherited constructors)
    
    //Access scalar array elements by name
    double& density();   
    double& energy();    
    double& psi();     
    const double& density() const;   
    const double& energy() const;    
    const double& psi() const;
    
    //Access vector3 parameters by name
    //B and momentum are accessed using Vector3View structure
    //Vector3View function allows direct read and write access
    //Vector3ConstView function is read only and can be included in other const member functions
    Vector3View B();
    Vector3View momentum();
    const Vector3ConstView B() const;
    const Vector3ConstView momentum() const;

    PrimitiveStateVector con_to_prim(const double gamma) const;
    
};

//NON MEMBER FUNCTIONS
//These create a new object and act on all children of the StateVectorTemplate structure
//Any permutation of the child structures are allowed
//All child structures store data as an array of doubles
//Output can be any of the child structures

template<typename T>
using SVT=StateVectorTemplate<T>;


//Vector addition and subtraction
template<typename DerivedL, typename DerivedR>
SVT<DerivedL> operator+(SVT<DerivedL> lhs,const SVT<DerivedR>& rhs){
    lhs+=rhs;
    return lhs;
}

template<typename DerivedL, typename DerivedR>
SVT<DerivedL> operator-(SVT<DerivedL> lhs,const SVT<DerivedR>& rhs){
    lhs-=rhs;
    return lhs;
}
//Scalar multiplication and division

template<typename Derived>
SVT<Derived> operator*(SVT<Derived> vector,const double scalar){
    vector*=scalar;
    return vector;
}
template<typename Derived>
SVT<Derived> operator*(const double scalar,SVT<Derived> vector){
    vector*=scalar;
    return vector;
}

template<typename Derived>
SVT<Derived> operator/(SVT<Derived> vector,const double scalar){
    vector/=scalar;
    return vector;
}
template<typename Derived>
SVT<Derived> operator/(const double scalar,SVT<Derived> vector){
    vector/=scalar;
    return vector;
}

#endif