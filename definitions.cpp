
#include "definitions.hpp"

std::vector <Polynomial> database;

void welcome()
{
    std::cout << "Program implements polynomial operations\n";
    std::cout << "Possible operations to perform are: addition, subtraction, multiplication and calculation of polynomial values (Horner's algorithms).\n\n\n";
}
int max(int a, int b)
{
    if(a > b)
        return a;
    else
        return b;
}

Polynomial::Polynomial(int power)
{
    coefficients = new double[power + 1];
    powerofpolynomial = power;
    for(int i = 0; i <= powerofpolynomial; i++)
        coefficients[i] = 0;
}
void Polynomial::createPolynomial(double arr[])
{
    for(int i = 0; i <= powerofpolynomial; i++)
    {
        coefficients[i] = arr[i];
    }
}
void Polynomial::printPolynomial()
{
    for(int i = 0; i < powerofpolynomial + 1; i++)
    {
        if(coefficients[i] == 0)
            continue;
        if(i != 0 && coefficients[i] > 0)
            std::cout << '+' << " ";
        std::cout << coefficients[i];
        if(i > 0)
            std::cout << "x^" << i;
        std::cout << " ";
    }
    std::cout << '\n';
}
Polynomial Polynomial::operator+(Polynomial p1)
{
    int maxi = max(powerofpolynomial, p1.powerofpolynomial);
    Polynomial sum(maxi);
    for(int i = 0; i <= p1.powerofpolynomial; i++)
        sum.coefficients[i] += p1.coefficients[i];
    for(int i = 0; i <= powerofpolynomial; i++)
        sum.coefficients[i] += coefficients[i];
    return sum;
}
void Polynomial::operator+=(Polynomial p1)
{
    if(p1.powerofpolynomial > powerofpolynomial)
    {
        double* coeffs = new double[p1.powerofpolynomial + 1];
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coeffs[i] = 0;
        for(int i = 0; i <= powerofpolynomial; i++)
            coeffs[i] += coefficients[i];
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coeffs[i] += p1.coefficients[i];
        double* save = coefficients;
        coefficients = coeffs;
        powerofpolynomial = p1.powerofpolynomial;
        delete save;
    }
    else
    {
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coefficients[i] += p1.coefficients[i];
    }
}
void Polynomial::operator-=(Polynomial p1)
{
    if(p1.powerofpolynomial > powerofpolynomial)
    {
        double* coeffs = new double[p1.powerofpolynomial + 1];
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coeffs[i] = 0;
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coeffs[i] -= p1.coefficients[i];
        for(int i = 0; i <= powerofpolynomial; i++)
            coeffs[i] += coefficients[i];
        double* save = coefficients;
        coefficients = coeffs;
        powerofpolynomial = p1.powerofpolynomial;
        delete save;
    }
    else
    {
        for(int i = 0; i <= p1.powerofpolynomial; i++)
            coefficients[i] -= p1.coefficients[i];
    }
}
Polynomial Polynomial::operator-(Polynomial p1)
{
    int maxi = max(powerofpolynomial, p1.powerofpolynomial);
    Polynomial temp(maxi);
    for(int i = 0; i <= powerofpolynomial; i++)
        temp.coefficients[i] += coefficients[i];
    for(int i = 0; i <= p1.powerofpolynomial; i++)
        temp.coefficients[i] -= p1.coefficients[i];
    return temp;
}
double Polynomial::horner(double x)
{
    double ans = coefficients[powerofpolynomial];
    for(int i = powerofpolynomial - 1; i >= 0; i--)
        ans = ans * x + coefficients[i];
    return ans;
}
void Polynomial::fft(std::vector<cd> & a, bool invert) 
{
    int n = a.size();
    if (n == 1)
        return;
    std::vector<cd> a0(n / 2), a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double ang = 2 * PI / n * (invert ? -1 : 1);
    cd w(1), wn(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n/2] = a0[i] - w * a1[i];
        if (invert) {
            a[i] /= 2;
            a[i + n/2] /= 2;
        }
        w *= wn;
    }
}
std::vector<double> Polynomial::multiply(std::vector<double> const& a, std::vector<double> const& b) 
{
    std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    int m = a.size() + b.size();
    while (n < m) 
        n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    fft(fa, true);
    std::vector<double> result(n);
    for (int i = 0; i < n; i++)
        result[i] = round(fa[i].real());
    return result;
}
Polynomial Polynomial::operator*(Polynomial p1)
{
    std::vector <double> a(powerofpolynomial + 1);
    std::vector <double> b(p1.powerofpolynomial + 1);
    for(int i = 0; i <= powerofpolynomial; i++)
        a[i] = coefficients[i];
    for(int i = 0; i <= p1.powerofpolynomial; i++)
        b[i] = p1.coefficients[i];
    std::vector <double> ans = multiply(a, b);
    Polynomial p2(powerofpolynomial + p1.powerofpolynomial);
    for(int i = 0; i <= powerofpolynomial + p1.powerofpolynomial; i++)
        p2.coefficients[i] = ans[i];
    return p2;
}
void Polynomial::operator*=(Polynomial p1)
{
    std::vector <double> a(powerofpolynomial + 1);
    std::vector <double> b(p1.powerofpolynomial + 1);
    for(int i = 0; i <= powerofpolynomial; i++)
        a[i] = coefficients[i];
    for(int i = 0; i <= p1.powerofpolynomial; i++)
        b[i] = p1.coefficients[i];
    std::vector <double> ans = multiply(a, b);
    double* save = coefficients;
    powerofpolynomial += p1.powerofpolynomial;
    coefficients = new double[powerofpolynomial + 1];
    for(int i = 0; i <= powerofpolynomial; i++)
        coefficients[i] = ans[i];
    delete save;
}




void helper()
{
    std::cout << "To add a polynomial to the database, enter 'add' and follow the instructions" << '\n';
    std::cout << "To obtain the coefficients of a given polynomial, enter 'get' and follow the instructions" << '\n';
    std::cout << "To create a new polynomial by adding two polynomials, enter '+' and follow the instructions" << '\n';
    std::cout << "To increase polynomial a by the coefficients of polynomial b, enter '+=' and follow the instructions" << '\n';
    std::cout << "To create a new polynomial by subtracting two polynomials, enter '-' and follow the instructions" << '\n';
    std::cout << "To reduce polynomial a from the coefficients of polynomial b, enter '-=' and follow the instructions" << '\n';
    std::cout << "To create a polynomial by multiplying two polynomials, enter '*' and follow the instructions" << '\n';
    std::cout << "To modify the value of polynomial a by the product of polynomial b, enter *= and follow the instructions" << '\n';
    std::cout << "To obtain the value of a polynomial at a point (Horner's method), enter 'value' and follow the instructions" << '\n';
}

void add()
{
    std::cout << "Function for adding a new polynomial to the database\n";
    std::cout << "Please enter the degree of the polynomial: \n";
    int power;
    std::cin >> power;
    Polynomial p1(power);
    double arr[power + 1];
    std::cout << "Please enter the coefficients of the polynomial: \n";
    for(int i = 0; i <= power; i++)
    {
        double coeff;
        std::cin >> coeff;
        arr[i] = coeff;
    }
    p1.createPolynomial(arr);
    database.push_back(p1);
    std::cout << "Polynomial: ";
    p1.printPolynomial();
    std::cout << "has been added to the database!\n\n";
}

void get()
{
    std::cout << "The number of polynomials in the database is: " << database.size() << '\n';
    std::cout << "Please select a polynomial number from 1 to " << database.size() << " to get its coefficients\n";
    int number;
    std::cin >> number;
    database[number - 1].printPolynomial();
    std::cout << '\n';
}

void sumPolynomial()
{
    std::cout << "Please select two polynomial numbers from the database to sum them. Please enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    Polynomial p1 = database[num1 -1] + database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    p1.printPolynomial();
    std::cout << "Do you want to add the obtained polynomial to the database? If so, please enter y, otherwise n.\n";
    char c;
    std::cin >> c;
    if(c == 'y')
    {
        database.push_back(p1);
        std::cout << "Polynomial has been added to the database! \n\n";
    }
}
void increasePolynomial()
{
    std::cout << "Please select two polynomial numbers from the database, the values of the first one will be increased by the coefficients of the second one. \nPlease enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    database[num1 - 1] += database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    database[num1 - 1].printPolynomial();
}
void subtractPolynomial()
{
    std::cout << "Please select two polynomial numbers from the database in order to subtract the second one from the first one. Please enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    Polynomial p1 = database[num1 -1] - database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    p1.printPolynomial();
    std::cout << "Do you want to add the obtained polynomial to the database? If so, please enter y, otherwise n.\n";
    char c;
    std::cin >> c;
    if(c == 'y')
    {
        database.push_back(p1);
        std::cout << "Polynomial has been added to the database! \n\n";
    }
}
void decreasePolynomial()
{
    std::cout << "Please select two polynomial numbers from the database, the values of the first one will be reduced by the coefficients of the second one. \nPlease enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    database[num1 - 1] -= database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    database[num1 - 1].printPolynomial();
}

void getValue()
{
    std::cout << "Please select a polynomial from the database, please enter a number from 1 to " << database.size() << '\n';
    int num;
    std::cin >> num;
    std::cout << "Chosen polynomial: ";
    database[num - 1].printPolynomial();
    std::cout << "Please provide an argument to get the function value that the function takes for the given argument: \n";
    double val;
    std::cin >> val;
    std::cout << "For the argument " << val << " the function takes the value: " << database[num - 1].horner(val) << '\n' << '\n';
}


void information()
{
    std::cout << "The current number of polynomials in the database is: " << database.size() << '\n';
    if(database.size() == 0)
    {
        std::cout << "Please add the first polynomial to the database.\n";
        add();
    }
    std::cout << "For information about available functions, please enter 'helper'." << '\n';
    std::cout << "Please enter the operation to be performed: \n";
}

void multiplyOne()
{
    std::cout << "Please select two polynomial numbers from the database, the values ​​of the first one will be multiplied by the coefficients of the second one. \nPlease enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    database[num1 - 1] *= database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    database[num1 - 1].printPolynomial();
}

void multiplyTwo()
{
    std::cout << "Please select two polynomial numbers from the database in order to multiply the first one with the second one. Please enter two digits from 1 to " << database.size() << '\n';
    int num1, num2;
    std::cin >> num1 >> num2;
    Polynomial p1 = database[num1 -1] * database[num2 - 1];
    std::cout << "The resulting polynomial is: ";
    p1.printPolynomial();
    std::cout << "Do you want to add the obtained polynomial to the database? If so, please enter y, otherwise n.\n";
    char c;
    std::cin >> c;
    if(c == 'y')
    {
        database.push_back(p1);
        std::cout << "The polynomial has been added to the database! \n\n";
    }
}

void start()
{
    welcome();
    while(1)
    {
        information();
        std::string s;
        std::cin >> s;
        if(s == "add")
            add();
        else if(s == "get")
            get();
        else if(s == "+")
            sumPolynomial();
        else if(s == "+=")
            increasePolynomial();
        else if(s == "-")
            subtractPolynomial();
        else if(s == "-=")
            decreasePolynomial();
        else if(s == "value")
            getValue();
        else if(s == "*=")
            multiplyOne();
        else if(s == "*")
            multiplyTwo();
        else 
            helper();
    }
}