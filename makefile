MWTransferArr	:	Arr_DF.o Astrophys.o ExtMath.o getparms.o GS.o IDLInterface.o Messages.o MWmain.o Std_DF.o
			g++ -shared -o MWTransferArr.so Arr_DF.o Astrophys.o ExtMath.o getparms.o GS.o IDLInterface.o Messages.o MWmain.o Std_DF.o
Arr_DF.o	:	Arr_DF.cpp ExtMath.h IDLInterface.h DF.h Astrophys.h Messages.h
			g++ -c -O3 -fPIC -D LINUX Arr_DF.cpp
Astrophys.o	:	Astrophys.cpp ExtMath.h Astrophys.h
			g++ -c -O3 -fPIC -D LINUX Astrophys.cpp
ExtMath.o	:	ExtMath.cpp ExtMath.h
			g++ -c -O3 -fPIC -D LINUX ExtMath.cpp
getparms.o	:	getparms.cpp
			g++ -c -O3 -fPIC -D LINUX getparms.cpp
GS.o		:	GS.cpp ExtMath.h IDLInterface.h DF.h GS.h Astrophys.h Messages.h
			g++ -c -O3 -fPIC -D LINUX GS.cpp
Messages.o	:	Messages.cpp
			g++ -c -O3 -fPIC -D LINUX Messages.cpp
MWmain.o	:	MWmain.cpp ExtMath.h IDLInterface.h Astrophys.h DF.h GS.h Messages.h
			g++ -c -O3 -fPIC -D LINUX MWmain.cpp
IDLInterface.o	:	IDLInterface.cpp Messages.h MWmain.h ExtMath.h IDLInterface.h
			g++ -c -O3 -fPIC -D LINUX IDLInterface.cpp
Std_DF.o	:	Std_DF.cpp ExtMath.h IDLInterface.h DF.h Astrophys.h Messages.h
			g++ -c -O3 -fPIC -D LINUX Std_DF.cpp
