#if defined(_MSC_VER)
#  include "jconfig_windows.h"
#elif defined(__APPLE__)
#  include "jconfig_mac.h"
#elif defined(__GNUC__)
#  include "jconfig_gcc.h"
#endif
