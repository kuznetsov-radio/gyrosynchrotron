#define M_PI 3.14159265358979323846

inline double sqr(double x)
{
 return x*x;
}

inline double sign(double x)
{
 return (x<0.0) ? -1.0 : 1.0;
}

#undef max
inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

#undef min
inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

#define finite _finite
#define dNaN (double(HUGE_VAL))

void FindBesselJ(double, int, double*, double*);
void FindBesselJ_WH(double, double, double*, double*);
double Erf(double);
double ExpBesselK(int, double);
double Gamma(double z);

class IntegrableFunction
{
 public:
 virtual double F(double x)=0;
};

class Arr2D
{
 public:
 virtual double A(int i, int j)=0;
};

#define JMAX 20
double qromb(IntegrableFunction *F, double a, double b, double EPS, int *err);
double trapzd(IntegrableFunction *F, double a, double b, int N);
double qrombLog(IntegrableFunction *F, double a, double b, double EPS, int *err);
double trapzdLog(IntegrableFunction *F, double a, double b, int N);
double IntTabulated(double *x, double *y, int N);
double IntTabulatedLog(double *x, double *y, int N);

#define MAXIT 20
#define BrentMAXIT 100
double SecantRoot(IntegrableFunction *F, double x1, double x2, double EPS);
double BrentRoot(IntegrableFunction *F, double x1, double x2, double EPS);

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

void LQInterpolate(double x, int N, double *x_arr, double *y_arr, double *y, double *y1);
void LQInterpolate2D(double x, double y, 
	                 int Nx, int Ny, double *x_arr, double *y_arr, double **f_arr, 
	                 double *f, double *f_x, double *f_y, double *f2_yy);