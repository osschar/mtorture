#ifndef MatriplexCommon_H
#define MatriplexCommon_H

#include <cmath>
#include <cstring>
#include <stdexcept>

// Use intrinsics version of code when available
#define  MPLEX_USE_INTRINSICS

namespace Matriplex
{
   typedef int idx_t;

   inline void align_check(const char* pref, void *adr)
   {
      printf("%s 0x%llx  -  modulo 64 = %d\n", pref, adr, (long long)adr%64);
   }
}

#endif
