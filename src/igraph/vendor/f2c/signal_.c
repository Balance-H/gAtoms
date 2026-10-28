#include "f2c.h"
#include "signal1.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

 ftnint
#ifdef KR_headers
signal_(sigp, proc) integer *sigp; sig_pf proc;
#else
signal_(integer *sigp, sig_pf proc)
#endif
{
	int sig;
	intptr_t prev;
	sig = (int)*sigp;

	/* signal() returns the previous handler as a function pointer.  Convert it
	   in two steps -- through void*, then through intptr_t -- so the pointer is
	   widened to its own width before being narrowed to f2c's legacy 32-bit
	   ftnint return type.  Behaviour preserving: the original code also narrowed
	   straight to ftnint, so no truncation is added or removed.  No igraph code
	   calls signal_ or uses its result. */
	prev = (intptr_t)(void *)signal(sig, proc);
	return (ftnint)prev;
	}
#ifdef __cplusplus
}
#endif
