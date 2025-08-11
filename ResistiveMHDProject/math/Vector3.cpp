#include "Vector3.h"

//ALL GENERIC METHODS ARE ENCOMPASSED IN TEMPLATE STRUC IN HEADER FILE
//These methods are all unique to particular variants of the Vector3Template
//Vector3 stores 3 doubles
//Vector3View stores 3 refs to doubles



//Constructor Definitions
//construct by parameters
Vector3::Vector3(double x, double y, double z)
        : data{{x,y,z}} {} //Initialiser
//convert from different vector type
template<typename Derived>
Vector3::Vector3(const Vector3Template<Derived>& other)
        : data{{other.x(),other.y(),other.z()}} {} 


Vector3View::Vector3View(double* x, double* y, double* z)
        : data{{x,y,z}} {} 

Vector3ConstView::Vector3ConstView(const double* x, const double* y, const double* z)
        : data{{x,y,z}} {} 

//Defining index access functions

double& Vector3::x(){return data[0];};
double& Vector3::y(){return data[1];};
double& Vector3::z(){return data[2];};

double& Vector3View::x(){return *data[0];};
double& Vector3View::y(){return *data[1];};
double& Vector3View::z(){return *data[2];};

const double& Vector3::x() const{return data[0];};
const double& Vector3::y() const{return data[1];};
const double& Vector3::z() const{return data[2];};

const double& Vector3View::x() const{return *data[0];};
const double& Vector3View::y() const{return *data[1];};
const double& Vector3View::z() const{return *data[2];};

const double& Vector3ConstView::x() const{return *data[0];};
const double& Vector3ConstView::y() const{return *data[1];};
const double& Vector3ConstView::z() const{return *data[2];};