## Originally written by Alexey Kuznetsov, April 2021
Fast codes for computing gyrosynchrotron radio emission; analytical and numerically defined electron distribution functions are supported. See the files CallingConventions.pdf and Diagram.pdf for more details and calling conventions, and the folder Examples for the usage examples.

## Modified some files to compile under MacOS (tested under 10.15.7 Catalina). Bin Chen, June 1, 2021
* Add a definition flag in makefile "-D OSX". If defined, change all calls for "malloc.h" to "stdlib.h", as the former is not available in MacOS. The flag is automatically added in makefile if make discovers that the system is running MacOS
* Modified makefile to include the "-std=c++11" flag to the g++ compiler. Otherwise the compiler make multiple complains such as
```
    ./DF.h:6:6: error: function definition does not declare parameters
	int N_intervals{};
```

* Note another error appears when compiling IDLinterface.cpp.  Need to install the "libomp" package using, e.g., "brew install libomp"
```
    g++ -c -std=c++11 -O3 -fPIC -D LINUX IDLinterface.cpp
    IDLinterface.cpp:11:10: fatal error: 'omp.h' file not found
    #include <omp.h>
```
