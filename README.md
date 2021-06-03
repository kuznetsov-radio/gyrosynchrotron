## Originally written by Alexey Kuznetsov, April 2021
Fast codes for computing gyrosynchrotron radio emission; analytical and numerically defined electron distribution functions are supported. See the files CallingConventions.pdf and Diagram.pdf for more details and calling conventions, and the folder Examples for the usage examples.

## Modified some files to compile under MacOS (tested under 10.15.7 Catalina and 11.2.2 Big Sur). Bin Chen, June 1, 2021
* Add a definition flag in makefile "-D OSX". If defined, change all calls for "malloc.h" to "stdlib.h", as the former is not available in MacOS. The flag is automatically added in makefile if make discovers that the system is running MacOS
* Modified makefile to include the "-std=c++11" flag to the g++ compiler. Otherwise the compiler make multiple complains such as
```
    ./DF.h:6:6: error: function definition does not declare parameters
	int N_intervals{};
```

See the files CallingConventions.pdf and Diagram.pdf for more details and calling conventions, and the folder Examples for the usage examples.

Compilation instructions:

for Windows users: use MS Visual Studio and open the gscodes.vcxproj project file;

for Linux or MacOS users: use make.

Note for MacOS users: you may need to install the parallel OMP library firstly, e.g., "brew install libomp" (or "arch -x86-64 brew install libomp" for Macs with ARM processor).
