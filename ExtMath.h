#pragma once

inline double sqr(double x)
{
 return x*x;
}

inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

inline int min(int a, int b)
{
 return (a<b) ? a : b;
}

inline int max(int a, int b)
{
 return (a>b) ? a : b;
}

inline double sign(double x)
{
 return (x<0.0) ? -1.0 : 1.0;
}

#ifndef LINUX
#define finite _finite
#else
#define finite isfinite
#endif

#define dNaN (double(HUGE_VAL))
#define JMAX 20
#define MAXIT 20
#define BrentMAXIT 100

class Spline
{
 int N;
 double *x_arr, *y_arr, *y2_arr;
 public:
 Spline(int _N, double *x, double *y);
 ~Spline();
 void Interpolate(double x, double *y, double *y1);
};

class Spline2D
{
 int Nx, Ny;
 double *x_arr, *y_arr;
 double **f_arr, **f2_xx_arr, **f2_yy_arr, **f4_xxyy_arr;
 public:
 Spline2D(int N_x, int N_y, double *x, double *y, double **f);
 ~Spline2D();
 void Interpolate(double x, double y, double *f, double *f_x, double *f_y, double *f2_yy);
};

class IntegrableFunction
{
 public:
 virtual double F(double x)=0;
};

void LQInterpolate(double x, int N, double *x_arr, double *y_arr, double *y, double *y1);
void LQInterpolate2D(double x, double y, 
	                 int Nx, int Ny, double *x_arr, double *y_arr, double **f_arr, 
	                 double *f, double *f_x, double *f_y, double *f2_yy);
double IntTabulated(double *x, double *y, int N);
double IntTabulatedLog(double *x, double *y, int N);
double ExpBesselK(int, double);
double qromb(IntegrableFunction *F, double a, double b, double EPS, int *err);
double trapzd(IntegrableFunction *F, double a, double b, int N);
double qrombLog(IntegrableFunction *F, double a, double b, double EPS, int *err);
double trapzdLog(IntegrableFunction *F, double a, double b, int N);
double SecantRoot(IntegrableFunction *F, double x1, double x2, double EPS);
double BrentRoot(IntegrableFunction *F, double x1, double x2, double EPS);
double InterpolateBilinear(double *arr, double i1, double i2, int N1, int N2, double missing);
double InterpolBilinear(double *arr, double *x1arr, double *x2arr, double x1, double x2, int N1, int N2);
double Erf(double);
void FindBesselJ(double, int, double*, double*);
void FindBesselJ_WH(double, double, double*, double*);
double Gamma(double z);