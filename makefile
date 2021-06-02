ifeq ($(OS), Windows_NT)
    target = MWTransferArr.dll
else
    UNAME_S := $(shell uname -s)
    target = MWTransferArr.so
    # add additional flag if the system is MacOS
    ifeq ($(UNAME_S), Darwin)
        osflag = -D OSX
    endif
endif
print-%  : ; @echo $* = $($*)

MWTransferArr	:	Arr_DF.o Coulomb.o ExtMath.o FF.o getparms.o GS.o IDLinterface.o Messages.o Plasma.o MWmain.o Std_DF.o Zeta.o
			g++ $^ -shared -o $(target)
Arr_DF.o	:	Arr_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Arr_DF.cpp
Coulomb.o	:	Coulomb.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Coulomb.cpp
ExtMath.o	:	ExtMath.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) ExtMath.cpp
FF.o		:	FF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) FF.cpp
getparms.o	:	getparms.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) getparms.cpp
GS.o		:	GS.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) GS.cpp
IDLinterface.o	:	IDLinterface.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) IDLinterface.cpp
Messages.o	:	Messages.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Messages.cpp
MWmain.o	:	MWmain.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) MWmain.cpp
Plasma.o	:	Plasma.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Plasma.cpp
Std_DF.o	:	Std_DF.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Std_DF.cpp
Zeta.o		:	Zeta.cpp Arr_DF.h Coulomb.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX $(osflag) Zeta.cpp
