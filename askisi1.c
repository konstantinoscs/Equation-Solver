#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double bolzano_method(double a, double b, double (*f)(double x), int *count);
double Newton_Raphson_method(double xn, double (*f)(double x),
    double (*df)(double x), int * count);
double f1(double x);
double df1(double x);
double f2(double x);
double df2(double x);
double power(double base, int exp);

int main (void)
{
    double xn = 0, prevxn=0, a=1.2356563, b=5.239480;
    double en = 0;
    int count=0;

    printf("For function f1 in[%lf, %lf]:\n\n", a,b);
    xn = Newton_Raphson_method(bolzano_method(a, b, &f1, &count), &f1, &df1, &count);
    printf ("Combination of Bolzano and NR method for function f1 found ");
    printf("solution xn=%lf with %d repetitions in total.\n\n", xn, count);

    a = -2.349093;
    b = 0.0987876;
    count = 0;

    printf("For function f1 in[%lf, %lf]:\n\n", a,b);
    xn = Newton_Raphson_method(bolzano_method(a, b, &f1, &count), &f1, &df1, &count);
    printf ("Combination of Bolzano and NR method for function f1 found ");
    printf("solution xn=%lf with %d repetitions in total.\n\n", xn, count);

    a = 1;
    b = 2;
    count = 0;

    printf("For function f2 in[%lf, %lf]:\n\n", a,b);
    xn = Newton_Raphson_method(bolzano_method(a, b, &f2, &count), &f2, &df2, &count);
    printf ("Combination of Bolzano and NR method for function f1 found ");
    printf("solution xn=%lf with %d repetitions in total.\n", xn, count);
}

double bolzano_method(double a, double b, double(*f)(double x), int* count)
{
    double xn=0, prevxn=0, en=0, preven=0;
    int rep = 0;

    assert ((*f)(a)*(*f)(b)<0.0); //solution in [a,b], bolzano theorem

    xn = (a+b)/2.0;
    printf("x0=%lf\n", xn);

    if ((*f)(xn) == 0)
        return xn;      //xn is a solution

    if ((*f)(a)*(*f)(xn)< 0.0)
    {
        b = xn;
    }
    else
    {
        assert ((*f)(b)*(*f)(xn) < 0.0);
        a = xn;
    }

    do {
        prevxn = xn;
        xn = (a+b)/2.0;

        if ((*f)(xn) == 0)
            break;      //xn is a solution

        if ((*f)(a)*(*f)(xn)< 0.0)
        {
            b = xn;
        }
        else
        {
            assert ((*f)(b)*(*f)(xn) < 0.0);
            a = xn;
        }

        preven =en;
        en = (xn - prevxn);
        if (en < 0.0)
            en *= -1;   //absolute value
        rep++;
        printf("for n = %d, en is %.10lf, xn is %.10lf\n", rep,en, xn);
    } while(en > 0.005);

    printf ("Bolzano method ended with xn %lf and error en=%.10lf\n", xn, en);
    printf ("Made %d repetitions\n", rep);
    (*count) += rep;
    return xn;
}

double Newton_Raphson_method(double starting_xn, double (*f)(double x),
    double (*df)(double x), int *count)
{
    double prevxn =0, xn = starting_xn, en=0;
    int rep=0;

    prevxn = xn;
    xn = xn - (*f)(xn)/(*df)(xn);
    en = (xn - prevxn);
    if (en < 0.0)
        en *= -1;   //absolute value
    printf("for n = %d, en is %lf, xn is %.15lf\n", rep+(*count), en, xn);

    while (en > 0.0000005)
    {
        prevxn = xn;
        xn = xn - (*f)(xn)/(*df)(xn);
        en = (xn - prevxn);
        if (en < 0.0)
            en *= -1;   //absolute value
        rep++;
        printf("for n = %d, en is %lf, xn is %.15lf\n", rep+(*count), en, xn);
    }
    printf("Newton-Raphson method ended with xn=%lf and error en=%.10lf\n", xn, en);
    printf("Made %d repetitions\n", rep);
    (*count) += rep;
    return xn;
}

double f1(double x)
{
    return (x+1)*(x+1)*(x+1)*(x-2);
}

double df1(double x)
{
    return 3*(x+1)*(x+1)*(x-2) + (x+1)*(x+1)*(x+1);
}

double f2(double x)
{
    return exp(x) - x*x - 2;
}

double df2(double x)
{
    return exp(x) - 2*x;
}

double power(double base, int exp)
{
    double ret = base;
    for (int i=0; i<exp; i++)
        ret*=ret;
    return ret;
}
