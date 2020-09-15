#if defined(_MSC_VER)
# pragma warning(push,1)
#endif

#include <cm3p/kwiml/abi.h>
#include <cm3p/kwiml/int.h>

/* Define to `int' if <sys/types.h> does not define. */
#define ssize_t KWIML_INT_intptr_t

/* sizeof(int *) */
#define SIZEOF_INT_P KWIML_ABI_SIZEOF_DATA_PTR

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1
