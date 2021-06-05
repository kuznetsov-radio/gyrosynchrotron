Fast codes for computing gyrosynchrotron radio emission; analytical and numerically defined electron distribution functions are supported.

Alexey Kuznetsov, April 2021.

See the files CallingConventions.pdf and Diagram.pdf for more details and calling conventions, and the folder Examples for the usage examples.

Compilation instructions:

For Windows users: use MS Visual Studio and open the gscodes.vcxproj project file; if necessary, adjust the paths to the IDL library (idl.lib) and header (idl_export.h) files in the project settings, according to your IDL installation.

For Linux or MacOS users: use make. Note for MacOS users: you may need to install the parallel OMP library firstly, e.g., "brew install libomp" (or "arch -x86-64 brew install libomp" for Macs with ARM processor).
