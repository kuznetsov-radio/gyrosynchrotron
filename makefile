MWTransferArr	:	Arr_DF.o Coulomb.o ExtMath.o FF.o getparms.o GS.o IDLinterface.o PythonInterface.o Messages.o Plasma.o MWmain.o Std_DF.o Zeta.o
			g++ $^ -shared -o MWTransferArr.so
Arr_DF.o	:	Arr_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Arr_DF.cpp
Coulomb.o	:	Coulomb.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Coulomb.cpp
ExtMath.o	:	ExtMath.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX ExtMath.cpp
FF.o		:	FF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX FF.cpp
getparms.o	:	getparms.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX getparms.cpp
GS.o		:	GS.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX GS.cpp
IDLinterface.o	:	IDLinterface.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX IDLinterface.cpp
PythonInterface.o :	PythonInterface.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX PythonInterface.cpp			
Messages.o	:	Messages.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Messages.cpp
MWmain.o	:	MWmain.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX MWmain.cpp
Plasma.o	:	Plasma.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Plasma.cpp
Std_DF.o	:	Std_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Std_DF.cpp
Zeta.o		:	Zeta.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Zeta.cpp