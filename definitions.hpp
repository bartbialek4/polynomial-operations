#include <iostream>
#include <complex>
#include <vector>

/*Needed for FFT purposes*/
using cd = std::complex<double>;
const double PI = acos(-1);



class Polynomial
{
private:
    double* coefficients;
    int powerofpolynomial;
public:
    Polynomial(int power);
    void createPolynomial(double arr[]);
    void printPolynomial();
    Polynomial operator+(Polynomial p1);
    void operator+=(Polynomial p1);
    Polynomial operator-(Polynomial p1);
    void operator-=(Polynomial p1);
    double horner(double x);
    void fft(std::vector<cd> & a, bool invert);
    std::vector<double> multiply(std::vector<double> const& a, std::vector<double> const& b);
    Polynomial operator*(Polynomial p1);
    void operator*=(Polynomial p1);
};




void welcome();
int max(int a, int b);
void helper();
void add();
void get();
void sumPolynomial();
void increasePolynomial();
void subtractPolynomial();
void decreasePolynomial();
void getValue();
void information();
void multiplyOne();
void multiplyTwo();
void start();
