MWTransferArr	:	Arr_DF.o Coulomb.o ExtMath.o FF.o getparms.o GS.o IDLinterface.o Messages.o Plasma.o MWmain.o Std_DF.o Zeta.o
			g++ -shared -o MWTransferArr.so Arr_DF.o Coulomb.o ExtMath.o FF.o getparms.o GS.o IDLinterface.o Messages.o Plasma.o MWmain.o Std_DF.o Zeta.o
Arr_DF.o	:	Arr_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Arr_DF.cpp
Coulomb.o	:	Coulomb.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Coulomb.cpp
ExtMath.o	:	ExtMath.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS ExtMath.cpp
FF.o		:	FF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS FF.cpp
getparms.o	:	getparms.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS getparms.cpp
GS.o		:	GS.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS GS.cpp
IDLinterface.o	:	IDLinterface.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS IDLinterface.cpp
Messages.o	:	Messages.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Messages.cpp
MWmain.o	:	MWmain.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS MWmain.cpp
Plasma.o	:	Plasma.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Plasma.cpp
Std_DF.o	:	Std_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Std_DF.cpp
Zeta.o		:	Zeta.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX -D MACOS Zeta.cpp
