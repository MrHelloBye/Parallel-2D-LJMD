#include "vec2.h"
#include "random.h"
//#include <cmath>
#include <math.h>

using std::cout; using std::endl;

vec2::vec2()
{
    components[0] = 0;
    components[1] = 0;
}

vec2::vec2(const vec2 &copy)
{
    components[0] = copy.x();
    components[1] = copy.y();
}

vec2::vec2(double x, double y)
{
    components[0] = x;
    components[1] = y;
}

void vec2::print()
{
    // Will print matlab syntax vector. Output will be like: [2.09, 5.3, 9.1];
    cout << "[" << components[0] << ", " << components[1] << "]" << endl;
}

void vec2::print(string name)
{
    // Will print matlab syntax vector with a name. Output will be like: A = [2.09, 5.3, 9.1];
    cout << name << " = ";
    print();
}

/*
vec2 vec2::cross(vec2 aVector)
{
    return vec2(y()*aVector.z()-z()*aVector.y(), z()*aVector.x()-x()*aVector.z(), x()*aVector.y()-y()*aVector.x());
}
*/

double vec2::lengthSquared()
{
    // Returns the square of the length (or norm) of the vector
    return components[0]*components[0]+components[1]*components[1];
}

double vec2::length()
{
    // Returns the length (or norm) of the vector
    return sqrt(lengthSquared());
}

void vec2::set(double x, double y)
{
    setX(x);
    setY(y);
}

void vec2::randomGaussian(double mean, double standardDeviation)
{
    components[0] = Random::nextGaussian(mean, standardDeviation);  //nextGaussian is in random.h
    components[1] = Random::nextGaussian(mean, standardDeviation);
}

vec2 &vec2::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    return *this;
}

vec2 &vec2::operator+=(vec2 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    return *this;
}

vec2 &vec2::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    return *this;
}

vec2 &vec2::operator*=(int rhs)   //add definition  for multiply by int
{
    components[0] *= rhs;
    components[1] *= rhs;
    return *this;
}

vec2 &vec2::operator*=(vec2 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    return *this;
}

vec2 &vec2::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    return *this;
}

vec2 &vec2::operator-=(vec2 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    return *this;
}

vec2 &vec2::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    return *this;
}

vec2 &vec2::operator/=(vec2 rhs)
{
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    return *this;
}

bool vec2::operator==(vec2 rhs)
{
    if(components[0]==rhs[0] && components[1]==rhs[1]){
    return true;
    }else return false;
}

std::ostream &operator<<(std::ostream &os, const vec2 &myVector) // Allows cout << myVector << endl;
{
    os << "[" << myVector.x() << ", " << myVector.y() << "]";
    return os;
}
