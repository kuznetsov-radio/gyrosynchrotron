#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "ExtMath.h"

void spline_init(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
 int i, k;
 double p, qn, sig, un, *u;
 u=(double*)malloc(sizeof(double)*n);
 if (!finite(yp1)) y2[0]=u[0]=0.0; 
 else 
 { 
  y2[0]=-0.5;
  u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
 }
 for (i=1; i<n-1; i++) 
 { 
  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p=sig*y2[i-1]+2.0;
  y2[i]=(sig-1.0)/p;
  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
 }
 if (!finite(ypn)) qn=un=0.0;
 else 
 { 
  qn=0.5;
  un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
 }
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for (k=n-2; k>=0; k--) y2[k]=y2[k]*y2[k+1]+u[k]; 
 free(u);
}

void spline_short(double *x, double *y, int n, double *y2)
{
 double dxl, dxr;
 double K1[3];

 dxl=x[1]-x[0];
 dxr=x[2]-x[1];
 K1[0]=-(2.0*dxl+dxr)/dxl/(dxl+dxr);
 K1[1]=(dxr+dxl)/(dxr*dxl);
 K1[2]=-dxl/dxr/(dxr+dxl);
 double y1l=K1[0]*y[0]+K1[1]*y[1]+K1[2]*y[2];

 dxl=x[n-2]-x[n-3];
 dxr=x[n-1]-x[n-2];
 K1[0]=dxr/dxl/(dxr+dxl);
 K1[1]=-(dxr+dxl)/(dxr*dxl);
 K1[2]=(2.0*dxr+dxl)/dxr/(dxr+dxl);
 double y1r=K1[0]*y[n-3]+K1[1]*y[n-2]+K1[2]*y[n-1];

 spline_init(x, y, n, y1l, y1r, y2);
}

Spline :: Spline(int _N, double *x, double *y)
{
 N=_N;

 x_arr=(double*)malloc(sizeof(double)*N);
 y_arr=(double*)malloc(sizeof(double)*N);
 y2_arr=(double*)malloc(sizeof(double)*N);

 for (int i=0; i<N; i++)
 {
  x_arr[i]=x[i];
  y_arr[i]=y[i];
 }

 spline_short(x_arr, y_arr, N, y2_arr);
}

void spline_interp(double *xa, double *ya, double *y2a, int n, double x, double *y, double *y1)
{
 int klo, khi, k;
 double h, b, a;

 if (x<=xa[0])
 {
  klo=0;
  khi=1;
 }
 else if (x>=xa[n-1])
 {
  klo=n-2;
  khi=n-1;
 }
 else
 {
  klo=0; 
  khi=n-1;
  while (khi-klo>1) 
  {
   k=(khi+klo)>>1;
   if (xa[k]>x) khi=k;
   else klo=k;
  } 
 }

 h=xa[khi]-xa[klo];
 a=(xa[khi]-x)/h; 
 b=(x-xa[klo])/h; 
 if (y)  *y =a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
 if (y1) *y1=(ya[khi]-ya[klo])/h+((1.0-3.0*a*a)*y2a[klo]+(3.0*b*b-1)*y2a[khi])*h/6.0;
}

Spline :: ~Spline()
{
 free(x_arr);
 free(y_arr);
 free(y2_arr);
}

void Spline :: Interpolate(double x, double *y, double *y1)
{
 spline_interp(x_arr, y_arr, y2_arr, N, x, y, y1);
}

Spline2D :: Spline2D(int N_x, int N_y, double *x, double *y, double **f)
{
 Nx=N_x;
 Ny=N_y;

 x_arr=(double*)malloc(sizeof(double)*Nx);
 y_arr=(double*)malloc(sizeof(double)*Ny);

 f_arr=(double**)malloc(sizeof(double*)*Nx);
 f2_xx_arr=(double**)malloc(sizeof(double*)*Nx);
 f2_yy_arr=(double**)malloc(sizeof(double*)*Nx);
 f4_xxyy_arr=(double**)malloc(sizeof(double*)*Nx);
 for (int i=0; i<Nx; i++)
 {
  f_arr[i]=(double*)malloc(sizeof(double)*Ny);
  f2_xx_arr[i]=(double*)malloc(sizeof(double)*Ny);
  f2_yy_arr[i]=(double*)malloc(sizeof(double)*Ny);
  f4_xxyy_arr[i]=(double*)malloc(sizeof(double)*Ny);
 }

 for (int i=0; i<Nx; i++) x_arr[i]=x[i];
 for (int j=0; j<Ny; j++) y_arr[j]=y[j];
 for (int i=0; i<Nx; i++) for (int j=0; j<Ny; j++) f_arr[i][j]=f[i][j];

 double *xprof=(double*)malloc(sizeof(double)*Nx);
 double *xprof2=(double*)malloc(sizeof(double)*Nx);
 for (int j=0; j<Ny; j++)
 {
  for (int i=0; i<Nx; i++) xprof[i]=f_arr[i][j];
  spline_short(x, xprof, Nx, xprof2);
  for (int i=0; i<Nx; i++) f2_xx_arr[i][j]=xprof2[i];
 }
 free(xprof);
 free(xprof2);

 for (int i=0; i<Nx; i++) 
 {
  spline_short(y, f_arr[i], Ny, f2_yy_arr[i]); 
  spline_short(y, f2_xx_arr[i], Ny, f4_xxyy_arr[i]);
 } 
}
 
Spline2D :: ~Spline2D()
{
 free(x_arr);
 free(y_arr);

 for (int i=0; i<Nx; i++)
 {
  free(f_arr[i]);
  free(f2_xx_arr[i]);
  free(f2_yy_arr[i]);
  free(f4_xxyy_arr[i]);
 }

 free(f_arr);
 free(f2_xx_arr);
 free(f2_yy_arr);
 free(f4_xxyy_arr);
}

void Spline2D :: Interpolate(double x, double y, 
                             double *f, double *f_x, double *f_y, double *f2_yy)
{
 int i1, i2;
 if (x<=x_arr[0])
 {
  i1=0;
  i2=1;
 }
 else if (x>=x_arr[Nx-1])
 {
  i1=Nx-2;
  i2=Nx-1;
 }
 else
 {
  i1=0; 
  i2=Nx-1;
  while (i2-i1>1) 
  {
   int k=(i1+i2)>>1;
   if (x_arr[k]>x) i2=k;
   else i1=k;
  } 
 }
 double hx=x_arr[i2]-x_arr[i1];
 double ax=(x_arr[i2]-x)/hx; 
 double bx=(x-x_arr[i1])/hx; 

 int j1, j2;
 if (y<=y_arr[0])
 {
  j1=0;
  j2=1;
 }
 else if (y>=y_arr[Ny-1])
 {
  j1=Ny-2;
  j2=Ny-1;
 }
 else
 {
  j1=0; 
  j2=Ny-1;
  while (j2-j1>1) 
  {
   int k=(j1+j2)>>1;
   if (y_arr[k]>y) j2=k;
   else j1=k;
  } 
 }
 double hy=y_arr[j2]-y_arr[j1];
 double ay=(y_arr[j2]-y)/hy; 
 double by=(y-y_arr[j1])/hy; 

 double f_lo=ay*f_arr[i1][j1]+by*f_arr[i1][j2]+((ay*ay*ay-ay)*f2_yy_arr[i1][j1]+(by*by*by-by)*f2_yy_arr[i1][j2])*(hy*hy)/6.0;
 double f_y_lo=(f_arr[i1][j2]-f_arr[i1][j1])/hy+((1.0-3.0*ay*ay)*f2_yy_arr[i1][j1]+(3.0*by*by-1)*f2_yy_arr[i1][j2])*hy/6.0;
 double f2_yy_lo=ay*f2_yy_arr[i1][j1]+by*f2_yy_arr[i1][j2];
 double f2_xx_lo=ay*f2_xx_arr[i1][j1]+by*f2_xx_arr[i1][j2]+((ay*ay*ay-ay)*f4_xxyy_arr[i1][j1]+(by*by*by-by)*f4_xxyy_arr[i1][j2])*(hy*hy)/6.0;
 double f3_xxy_lo=(f2_xx_arr[i1][j2]-f2_xx_arr[i1][j1])/hy+((1.0-3.0*ay*ay)*f4_xxyy_arr[i1][j1]+(3.0*by*by-1)*f4_xxyy_arr[i1][j2])*hy/6.0;
 double f4_xxyy_lo=ay*f4_xxyy_arr[i1][j1]+by*f4_xxyy_arr[i1][j2];

 double f_hi=ay*f_arr[i2][j1]+by*f_arr[i2][j2]+((ay*ay*ay-ay)*f2_yy_arr[i2][j1]+(by*by*by-by)*f2_yy_arr[i2][j2])*(hy*hy)/6.0;
 double f_y_hi=(f_arr[i2][j2]-f_arr[i2][j1])/hy+((1.0-3.0*ay*ay)*f2_yy_arr[i2][j1]+(3.0*by*by-1)*f2_yy_arr[i2][j2])*hy/6.0;
 double f2_yy_hi=ay*f2_yy_arr[i2][j1]+by*f2_yy_arr[i2][j2];
 double f2_xx_hi=ay*f2_xx_arr[i2][j1]+by*f2_xx_arr[i2][j2]+((ay*ay*ay-ay)*f4_xxyy_arr[i2][j1]+(by*by*by-by)*f4_xxyy_arr[i2][j2])*(hy*hy)/6.0;
 double f3_xxy_hi=(f2_xx_arr[i2][j2]-f2_xx_arr[i2][j1])/hy+((1.0-3.0*ay*ay)*f4_xxyy_arr[i2][j1]+(3.0*by*by-1)*f4_xxyy_arr[i2][j2])*hy/6.0;
 double f4_xxyy_hi=ay*f4_xxyy_arr[i2][j1]+by*f4_xxyy_arr[i2][j2];

 if (f) *f=ax*f_lo+bx*f_hi+((ax*ax*ax-ax)*f2_xx_lo+(bx*bx*bx-bx)*f2_xx_hi)*(hx*hx)/6.0;
 if (f_x) *f_x=(f_hi-f_lo)/hx+((1.0-3.0*ax*ax)*f2_xx_lo+(3.0*bx*bx-1)*f2_xx_hi)*hx/6.0;
 if (f_y) *f_y=ax*f_y_lo+bx*f_y_hi+((ax*ax*ax-ax)*f3_xxy_lo+(bx*bx*bx-bx)*f3_xxy_hi)*(hx*hx)/6.0;
 if (f2_yy) *f2_yy=ax*f2_yy_lo+bx*f2_yy_hi+((ax*ax*ax-ax)*f4_xxyy_lo+(bx*bx*bx-bx)*f4_xxyy_hi)*(hx*hx)/6.0;
}

void CreateLQInterpolationKernel(double x, int N, double *x_arr, 
	                             int *i1, int *i2, double *K, double *K1, double *K2)
{
 int i0;

 if (x<=x_arr[0]) i0=0;
 else if (x>=x_arr[N-1]) i0=N-1;
 else
 {
  int k1=0;
  int k2=N-1;
  while ((k2-k1)>1)
  {
   int l=(k1+k2) >> 1;
   if (x_arr[l]>x) k2=l;
   else k1=l;
  }
  i0=((x-x_arr[k1])<(x_arr[k2]-x)) ? k1 : k2;
 }

 double dx=x-x_arr[i0];
 double dxl, dxr;

 if (i0==0)
 {
  *i1=0;
  *i2=2;
  dxl=x_arr[1]-x_arr[0];
  dxr=x_arr[2]-x_arr[1];
  K[0]=1.0-dx/dxl;
  K[1]=dx/dxl;
  K[2]=0.0;
  if (K1)
  {
   K1[0]=(2.0*dx-2.0*dxl-dxr)/dxl/(dxl+dxr);
   K1[1]=-(2.0*dx-dxr-dxl)/(dxr*dxl);
   K1[2]=(2.0*dx-dxl)/dxr/(dxr+dxl);
  }
 }
 else if (i0==(N-1))
 {
  *i1=N-3;
  *i2=N-1;
  dxl=x_arr[N-2]-x_arr[N-3];
  dxr=x_arr[N-1]-x_arr[N-2];
  K[0]=0.0;
  K[1]=-dx/dxr;
  K[2]=1.0+dx/dxr;
  if (K1)
  {
   K1[0]=(2.0*dx+dxr)/dxl/(dxr+dxl);
   K1[1]=-(2.0*dx+dxr+dxl)/(dxr*dxl);
   K1[2]=(2.0*dx+2.0*dxr+dxl)/dxr/(dxr+dxl);
  }
 }
 else
 {
  *i1=i0-1;
  *i2=i0+1;
  dxr=x_arr[i0+1]-x_arr[i0];
  dxl=x_arr[i0]-x_arr[i0-1];
  if (dx>0)
  {
   K[0]=0.0;
   K[1]=1.0-dx/dxr;
   K[2]=dx/dxr;
  }
  else
  {
   K[0]=-dx/dxl;
   K[1]=1.0+dx/dxl;
   K[2]=0.0;
  }
  if (K1)
  {
   K1[0]=(2.0*dx-dxr)/dxl/(dxr+dxl);
   K1[1]=-(2.0*dx-dxr+dxl)/(dxr*dxl);
   K1[2]=(2.0*dx+dxl)/dxr/(dxr+dxl);
  }
 }

 if (K2)
 {
  K2[0]=2.0/dxl/(dxr+dxl);
  K2[1]=-2.0/(dxr*dxl);
  K2[2]=2.0/dxr/(dxr+dxl);
 }
}

void LQInterpolate(double x, int N, double *x_arr, double *y_arr, double *y, double *y1)
{
 double K[3], K1[3];
 int i1, i2;

 CreateLQInterpolationKernel(x, N, x_arr, &i1, &i2, K, K1, 0);

 *y=0.0;
 *y1=0.0;
 for (int i=i1; i<=i2; i++)
 {
  *y+=(y_arr[i]*K[i-i1]);
  *y1+=(y_arr[i]*K1[i-i1]);
 }
}

void LQInterpolate2D(double x, double y, 
	                 int Nx, int Ny, double *x_arr, double *y_arr, double **f_arr, 
	                 double *f, double *f_x, double *f_y, double *f2_yy)
{
 double Kx[3], Ky[3], Kx1[3], Ky1[3], Ky2[3];
 int i1, i2, j1, j2;

 CreateLQInterpolationKernel(x, Nx, x_arr, &i1, &i2, Kx, f_x ? Kx1 : 0, 0);
 CreateLQInterpolationKernel(y, Ny, y_arr, &j1, &j2, Ky, f_y ? Ky1 : 0, f2_yy ? Ky2 : 0);

 if (f) *f=0.0;
 if (f_x) *f_x=0.0;
 if (f_y) *f_y=0.0;
 if (f2_yy) *f2_yy=0.0;
 for (int i=i1; i<=i2; i++) for (int j=j1; j<=j2; j++)
 {
  if (f) *f+=(f_arr[i][j]*Kx[i-i1]*Ky[j-j1]);
  if (f_x) *f_x+=(f_arr[i][j]*Kx1[i-i1]*Ky[j-j1]);
  if (f_y) *f_y+=(f_arr[i][j]*Kx[i-i1]*Ky1[j-j1]);
  if (f2_yy) *f2_yy+=(f_arr[i][j]*Kx[i-i1]*Ky2[j-j1]); 
 }
}

double IntTabulated(double *x, double *y, int N)
{
 double s=0.0;
 for (int i=0; i<(N-1); i++) s+=(y[i]+y[i+1])*fabs(x[i+1]-x[i])/2;
 return s;
}

double IntTabulatedLog(double *t, double *u, int N)
{
 double s=0.0;
 for (int i=0; i<(N-1); i++) 
  s+=(exp(t[i+1]+u[i+1])-exp(t[i]+u[i]))/((t[i+1]+u[i+1])-(t[i]+u[i]))*(t[i+1]-t[i]);
 return s;
}

double bessi0(double x)
{
 double ax, ans;
 double y; 
 if ((ax=fabs(x))<3.75) 
 { 
  y=x/3.75;
  y*=y;
  ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
      +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
 } 
 else 
 {
  y=3.75/ax;
  ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
      +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
      +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
      +y*0.392377e-2))))))));
 }
 return ans;
}

double bessi1(double x)
{
 double ax, ans;
 double y; 
 if ((ax=fabs(x))<3.75) 
 { 
  y=x/3.75;
  y*=y;
  ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
      +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
 } 
 else 
 {
  y=3.75/ax;
  ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
      -y*0.420059e-2));
  ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
      +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
  ans*=(exp(ax)/sqrt(ax));
 }
 return (x<0.0) ? -ans : ans;
}

double expbessk0(double x)
{
 double y, ans; 
 if (x<=2.0) 
 { 
  y=x*x/4.0;
  ans=((-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
      +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
      +y*(0.10750e-3+y*0.74e-5)))))))*exp(x);
 } 
 else 
 {
  y=2.0/x;
  ans=(1.0/sqrt(x))*(1.25331414+y*(-0.7832358e-1
      +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
      +y*(-0.251540e-2+y*0.53208e-3))))));
 }
 return ans;
}

double expbessk1(double x)
{
 double y, ans;
 if (x<=2.0) 
 { 
  y=x*x/4.0;
  ans=((log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
      +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
      +y*(-0.110404e-2+y*(-0.4686e-4))))))))*exp(x);
 } 
 else 
 {
  y=2.0/x;
  ans=(1.0/sqrt(x))*(1.25331414+y*(0.23498619
      +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
      +y*(0.325614e-2+y*(-0.68245e-3)))))));
 }
 return ans;
}

double ExpBesselK(int n, double x)
//returns exp(x)*K_n(x)
{
 int j;
 double bk, bkm, bkp, tox;
 tox=2.0/x;
 bkm=expbessk0(x); 
 bk=expbessk1(x);
 for (j=1; j<n; j++) 
 { 
  bkp=bkm+j*tox*bk;
  bkm=bk;
  bk=bkp;
 }                  
 return bk;
}

void polint(double *xa, double *ya, int n, double x, double *y, double *dy)
{
 int i, m, ns=1;
 double den, dif, dift, ho, hp, w;
 double c[10], d[10]; //fixed-size arrays; n must be <=10
 dif=fabs(x-xa[1]);
 for (i=1; i<=n; i++)
 {
  if ((dift=fabs(x-xa[i]))<dif)
  {
   ns=i;
   dif=dift;
  }
  c[i]=ya[i];
  d[i]=ya[i];
 }
 *y=ya[ns--];
 for (m=1; m<n; m++)
 {
  for (i=1; i<=n-m; i++)
  {
   ho=xa[i]-x;
   hp=xa[i+m]-x;
   w=c[i+1]-d[i];
   den=ho-hp; 
   den=w/den;
   d[i]=hp*den;
   c[i]=ho*den;
  }
  *y+=(*dy=(2*ns<(n-m) ? c[ns+1] : d[ns--]));
 }
}

double trapzdQ(IntegrableFunction *F, double a, double b, int n, double *s)
{
 double x, tnm, sum, del;
 int it, j;
 if (n==1)
 {
  return ((*s)=0.5*(b-a)*(F->F(a)+F->F(b)));
 }
 else
 {
  for (it=1, j=1; j<n-1; j++) it<<=1;
  tnm=it;
  del=(b-a)/tnm;
  x=a+0.5*del;
  for (sum=0.0, j=1; j<=it; j++, x+=del) sum+=F->F(x);
  *s=0.5*(*s+(b-a)*sum/tnm);
  return *s;
 }
}

#define JMAXP (JMAX+1)
#define romK 6

double qromb(IntegrableFunction *F, double a, double b, double EPS, int *err)
{
 *err=0;               
 double ts, ss, dss;
 double s[JMAXP], h[JMAXP+1];
 int j;
 h[1]=1.0;
 for (j=1; j<=JMAX; j++)
 {
  s[j]=trapzdQ(F, a, b, j, &ts);
  if (j>=romK)
  {
   polint(&h[j-romK], &s[j-romK], romK, 0.0, &ss, &dss);
   if (fabs(dss)<=EPS*fabs(ss) || !finite(ss)) return ss;
  }
  h[j+1]=0.25*h[j];
 }
 *err=1; 
 return ss;
}

double trapzd(IntegrableFunction *F, double a, double b, int N)
{
 double s=0.0;
 double x=a;
 double dx=(b-a)/N;
 for (int i=0; i<=N; i++)
 {
  double u=F->F(x);
  if ((i==0) || (i==N)) u/=2;
  s+=u;
  x+=dx;
 }
 return s*dx;
}

class IntegrableFunctionLog : public IntegrableFunction
{
 public:
 IntegrableFunction *oldF{};
 double F(double t);
};

double IntegrableFunctionLog :: F(double t)
{
 double x=exp(t);
 return oldF->F(x)*x;
}

double qrombLog(IntegrableFunction *F, double a, double b, double EPS, int *err)
{
 IntegrableFunctionLog ifl;
 ifl.oldF=F;
 return qromb(&ifl, log(a), log(b), EPS, err);
}

double trapzdLog(IntegrableFunction *F, double a, double b, int N)
{
 IntegrableFunctionLog ifl;
 ifl.oldF=F;
 return trapzd(&ifl, log(a), log(b), N);
}

double InterpolateBilinear(double *arr, double i1, double i2, int N1, int N2, double missing)
/* Interpolation on an equidistant grid (like the IDL interpolate function),
   i1 and i2 - fractional indices of the required point. */
{
 if (i1<0 || i1>(N1-1) || i2<0 || i2>(N2-1)) return missing;

 int j=int(i1);
 int k=int(i2);
 double t=i1-j;
 double u=i2-k;

 double y1=arr[N2*j+k];
 double y2=arr[N2*(j+1)+k];
 double y3=arr[N2*(j+1)+k+1];
 double y4=arr[N2*j+k+1];

 return (1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
}

double InterpolBilinear(double *arr, double *x1arr, double *x2arr, double x1, double x2, int N1, int N2)
/* Interpolation on an arbitrary grid (like the IDL interpol function). 
   Performs extrapolation, if the point is outside the range. */   
{
 int j, j1, k, k1, l;

 if (x1<x1arr[0])
 {
  j=0;
  j1=1;
 }
 else if (x1>x1arr[N1-1])
 {
  j=N1-2;
  j1=N1-1;
 }
 else
 {
  j=0;
  j1=N1-1;
  while ((j1-j)>1)
  {
   l=(j1+j) >> 1;
   if (x1arr[l]>x1) j1=l;
   else j=l;
  }
 }
 double dx1=x1arr[j1]-x1arr[j];
 double t=(x1-x1arr[j])/dx1;

 if (x2<x2arr[0])
 {
  k=0;
  k1=1;
 }
 else if (x2>x2arr[N2-1])
 {
  k=N2-2;
  k1=N2-1;
 }
 else
 {
  k=0;
  k1=N2-1;
  while ((k1-k)>1)
  {
   l=(k1+k) >> 1;
   if (x2arr[l]>x2) k1=l;
   else k=l;
  }
 }
 double dx2=x2arr[k1]-x2arr[k];
 double u=(x2-x2arr[k])/dx2;
                                                           
 double y1=arr[N2*j+k];
 double y2=arr[N2*j1+k];
 double y3=arr[N2*j1+k1];
 double y4=arr[N2*j+k1];

 return (1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
}

double ErfC(double x)
//Returns the complementary error function erfc(x)
{
 double t, z, ans;
 z=fabs(x);
 t=1.0/(1.0+0.5*z);
 ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
     t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
     t*(-0.82215223+t*0.17087277)))))))));
 return (x>=0.0) ? ans : 2.0-ans;
}

double SecantRoot(IntegrableFunction *F, double x1, double x2, double EPS)
{
 double fl, f, dx, swap, xl, rts;

 fl=F->F(x1);
 f=F->F(x2);
 if (fabs(fl)<fabs(f)) 
 { 
  rts=x1; 
  xl=x2;
  swap=fl;
  fl=f;
  f=swap;
 } 
 else 
 {
  xl=x1;
  rts=x2;
 }

 int j=0;
 do
 { 
  dx=(xl-rts)*f/(f-fl); 
  xl=rts;
  fl=f;
  rts+=dx;
  f=F->F(rts);
  j++;
 }
 while (fabs(dx)>EPS && f!=0.0 && j<MAXIT);

 return (j<MAXIT) ? rts : dNaN;
}

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double BrentRoot(IntegrableFunction *F, double x1, double x2, double tol)
{
 int iter;
 double a=x1, b=x2, c=x2, d, e, min1, min2;
 double fa=F->F(a), fb=F->F(b), fc, p, q, r, s, tol1, xm;

 if (!finite(fa) || !finite(fb)) return dNaN;

 if (fa*fb>0.0) return dNaN;
 else
 {
  fc=fb;
  for (iter=1; iter<=BrentMAXIT; iter++) 
  {
   if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) 
   {
    c=a; 
    fc=fa;
    e=d=b-a;
   }
   if (fabs(fc)<fabs(fb)) 
   {
    a=b;
    b=c;
    c=a;
    fa=fb;
    fb=fc;
    fc=fa;
   }
   tol1=0.5*tol; 
   xm=0.5*(c-b);
   if (fabs(xm)<=tol1 || fb==0.0) return b;
   if (fabs(e)>=tol1 && fabs(fa)>fabs(fb)) 
   {
    s=fb/fa; 
    if (a==c) 
	{
     p=2.0*xm*s;
     q=1.0-s;
    } 
	else 
	{
     q=fa/fc;
     r=fb/fc;
     p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
     q=(q-1.0)*(r-1.0)*(s-1.0);
    }
    if (p>0.0) q=-q; 
    p=fabs(p);
    min1=3.0*xm*q-fabs(tol1*q);
    min2=fabs(e*q);
    if (2.0*p<(min1<min2 ? min1 : min2)) 
	{
     e=d; 
     d=p/q; 
    } 
	else 
	{
     d=xm; 
     e=d;
    }
   } 
   else 
   { 
    d=xm;
    e=d;
   }
   a=b; 
   fa=fb;
   if (fabs(d)>tol1) b+=d;
   else b+=SIGN(tol1, xm);
   fb=F->F(b); 
   if (!finite(fb)) return dNaN;
  }
  return dNaN; 
 }
}

double Erf(double x)
//Returns the error function erf(x)
{
 return 1.0-ErfC(x);
}

double bessj0(double x)
{
 double ax, z;
 double xx, y, ans, ans1, ans2;
 if ((ax=fabs(x))<8.0)
 {
  y=x*x;
  ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
  ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
  ans=ans1/ans2;
 }
 else
 {
  z=8.0/ax;
  y=z*z;
  xx=ax-0.785398164;
  ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
  ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)));
  ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
 }
 return ans;
}

double bessj1(double x)
{
 double ax, z;
 double xx, y, ans, ans1, ans2;
 if ((ax=fabs(x))<8.0)
 {
  y=x*x;
  ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
  ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
  ans=ans1/ans2;
 }
 else
 {
  z=8.0/ax;
  y=z*z;
  xx=ax-2.356194491;
  ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
  ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
  ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  if (x<0.0) ans=-ans;
 }
 return ans;
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
void FindBesselJ(double x, int n, double *Js, double *Js1)
{
 int j, jsum, m;
 double ax, bj, bjm, bjp, sum, tox, Bsn, Bsn1;
 if (n==1)
 {
  Bsn=bessj1(x);
  Bsn1=bessj0(x);
 }
 else
 {
  ax=fabs(x);
  if (ax>(double)n)
  {
   tox=2.0/ax;
   bjm=bessj0(ax);
   bj=bessj1(ax);
   for (j=1; j<n; j++)
   {
    bjp=tox*bj*j-bjm;
    bjm=bj;
    bj=bjp;
   }
   Bsn=bj;
   Bsn1=bjm;
  }
  else
  {
   tox=2.0/ax;
   m=2*((n+(int)sqrt(ACC*n))/2);
   jsum=0;
   bjp=Bsn=Bsn1=sum=0.0;
   bj=1.0;
   for (j=m; j>0; j--)
   {
    bjm=tox*bj*j-bjp;
    bjp=bj;
    bj=bjm;
    while (fabs(bj)>BIGNO)
    {
     bj*=BIGNI;
     bjp*=BIGNI;
     Bsn*=BIGNI;
     Bsn1*=BIGNI;
     sum*=BIGNI;
    }
    if (jsum) sum+=bj;
    jsum=jsum==0;
    if (j==n) Bsn=bjp;
    if (j==(n-1)) Bsn1=bjp;
   }
   sum=2.0*sum-bj;
   Bsn/=sum;
   Bsn1/=sum;
   if (n==2) Bsn1=bessj1(x);
  }
 }
 *Js=Bsn;
 *Js1=Bsn1-Bsn*n/x;
 if (x<0.0)
 {
  if (n & 1) *Js=-(*Js);
  else *Js1=-(*Js);
 }
 return;
}

void FindBesselJ_WH(double Sx, double S, double *JS, double *JS1)
{
 #define A 0.503297
 #define B 1.193000
 const double p16=1.0/6;
 const double pm23=-2.0/3;

 double x=Sx/S;
 double t1=sqrt(1.0-sqr(x));
 double t2=t1*t1*t1;
 double a=pow(t2+A/S, p16);
 double b=pow(t2+B/S, p16)*(1.0-0.2*pow(S, pm23));
 double F=a*b;
 double Z=x*exp(t1)/(1.0+t1);
 *JS=pow(Z, S)/sqrt(2.0*M_PI*S)/a;
 *JS1=F*(*JS)/x;
}

double LnGamma(double xx)
{
 double x, y, tmp, ser;
 static double cof[6]={ 76.18009172947146,
	                   -86.50532032941677,
                        24.01409824083091,
					   -1.231739572450155,
                        0.1208650973866179e-2,
					-0.5395239384953e-5};
 int j;
 y=x=xx;
 tmp=x+5.5;
 tmp-=(x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0; j<=5; j++) ser+=cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}

double Gamma(double z)
{
 return exp(LnGamma(z));
}