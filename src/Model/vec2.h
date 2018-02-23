#ifndef VEC2_H
#define VEC2_H
#include <string>
#include <vector>
#include <iostream>
using std::string;

class vec2
{
public:
    double components[2];
    vec2();
    vec2(vec2 const&copy);
    vec2(double x, double y);
    void print();
    void print(string name);
    vec2 cross(vec2 aVector);
    double lengthSquared();
    double length();
    double x() const { return components[0]; }
    double y() const { return components[1]; }
    void setX(double x) { components[0] = x; }
    void setY(double y) { components[1] = y; }
    void set(double x, double y);
    void zeros() { components[0] = 0; components[1] = 0; }
    double &operator()(int index) { return components[index]; } // Allows access like myVector(0)
    double &operator[](int index) { return components[index]; } // Allows access like myVector[0]
    void randomGaussian(double mean, double standardDeviation);
    vec2 &operator+=(double rhs); // Componentwise addition with scalar
    vec2 &operator+=(vec2 rhs);   // Componentwise addition with vector
    vec2 &operator*=(double rhs); // Componentwise multiplication with scalar
    vec2 &operator*=(vec2 rhs);   // Componentwise multiplicationwith vector
    vec2 &operator*=(int rhs);//Componentwise multiplication with int scalar
    vec2 &operator-=(double rhs); // Componentwise subtraction with scalar
    vec2 &operator-=(vec2 rhs);   // Componentwise subtraction with vector
    vec2 &operator/=(double rhs); // Componentwise division with scalar
    vec2 &operator/=(vec2 rhs);   // Componentwise division with vector
    bool operator==(vec2 rhs);  //Check if two vectors are equal componentwise
    friend std::ostream& operator<<(std::ostream& os, const vec2& myVector); // Allows cout << myVector << endl;
};

inline vec2 operator+(vec2 lhs, double rhs) {
    lhs += rhs;
    return lhs;
}

inline vec2 operator+(double lhs, vec2 rhs) {
    rhs += lhs;
    return rhs;
}

inline vec2 operator+(vec2 lhs, vec2 rhs) {
    lhs += rhs;
    return lhs;
}


inline vec2 operator-(vec2 lhs, double rhs) {
    lhs -= rhs;
    return lhs;
}

inline vec2 operator-(double lhs, vec2 rhs) {
    rhs -= lhs;
    return rhs;
}

inline vec2 operator-(vec2 lhs, vec2 rhs) {
    lhs -= rhs;
    return lhs;
}


inline vec2 operator*(vec2 lhs, double rhs) {
    lhs *= rhs;
    return lhs;
}

inline vec2 operator*(double lhs, vec2 rhs) {
    rhs *= lhs;
    return rhs;
}

inline vec2 operator*(vec2 lhs, vec2 rhs) {
    lhs *= rhs;
    return lhs;
}


inline vec2 operator/(vec2 lhs, double rhs) {
    lhs /= rhs;
    return lhs;
}

inline vec2 operator/(double lhs, vec2 rhs) {
    rhs /= lhs;
    return rhs;
}

inline vec2 operator/(vec2 lhs, vec2 rhs) {
    lhs /= rhs;
    return lhs;
}

#endif // vec2_H
