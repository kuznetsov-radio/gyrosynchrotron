#include <Windows.h>
#include "SliceWin.h"

BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
 SYSTEM_INFO si;
 GetSystemInfo(&si);
 ProcessorNumber=si.dwNumberOfProcessors;

 return TRUE;
}

