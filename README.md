Fast codes for computing gyrosynchrotron radio emission; analytical and numerically defined electron distribution functions are supported.

Alexey Kuznetsov, April 2021.

See the files CallingConventions.pdf and Diagram.pdf for more details and calling conventions, and the folder Examples for the usage examples.

For Python users, one need to import pygs located in the Examples folder and then call the respective functions, e.g., pygs.get_mw(). Similar calling conventions apply. [Updated by B. Chen, June 2021]

Compilation instructions:

for Windows users: use MS Visual Studio and open the gscodes.vcxproj project file;

for Linux or MacOS users: use make.

Note for MacOS users: you may need to install the parallel OMP library firstly, e.g., "brew install libomp" (or "arch -x86-64 brew install libomp" for Macs with ARM processor).


