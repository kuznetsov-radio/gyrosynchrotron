#include <stdarg.h>
#include <stdio.h>

#ifndef LINUX
#include <idl_export.h>

void IDLmsg(const char *fmt, ...)
{
 char arr[256];

 va_list argptr;
 va_start(argptr, fmt);
 vsprintf_s(arr, 256, fmt, argptr);
 va_end(argptr);

 IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, arr);
}

#else

void IDLmsg(const char *fmt, ...)
{}

#endif

int LOGinit=0;
char LOGname[32];

void LOGout(const char *fmt, ...)
{
 char arr[256];

 va_list argptr;
 va_start(argptr, fmt);
 #ifndef LINUX
 vsprintf_s(arr, 256, fmt, argptr);
 #else
 vsnprintf(arr, 256, fmt, argptr);
 #endif
 va_end(argptr);

 FILE *f;

 if (!LOGinit)
 {
  int i=0;
  int OK=0;
  do
  {
   sprintf(LOGname, "GS%04d.log", i);
   f=fopen(LOGname, "r");
   if (f)
   {
	fclose(f);
	i++;
   }
   else
   {
	f=fopen(LOGname, "w");
	OK=1;
   }
  } while (!OK);
  LOGinit=1;
 }
 else f=fopen(LOGname, "a");

 fprintf(f, "%s\n", arr);
 fclose(f);
}
