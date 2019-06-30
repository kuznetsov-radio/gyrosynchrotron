extern "C" __declspec(dllexport) double GET_MW(int argc, void **argv);

#define D1(arr, s, i) (((i)<(s)) ? arr[i] : 0)
#define D2(arr, s1, i1, i2) arr[(i1)+(i2)*(s1)]
#define D3(arr, s1, s2, i1, i2, i3) arr[(i1)+((i2)+(i3)*(s2))*(s1)]

class IDL_DFarr : public Arr2D
{
 int Nsteps, NE, Nmu;
 int step;
 double *f;
 public:
 IDL_DFarr(int _Nsteps, int _NE, int _Nmu, int _step, double *_f);
 double A(int i, int j);
};