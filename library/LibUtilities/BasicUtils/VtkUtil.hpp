#include <vtkVersionMacros.h>

#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0 && VTK_PATCH_VERSION == 0
#include <vtkMathConfigure.h>
#undef VTK_HAS_ISNAN
#undef VTK_HAS_ISINF
#endif
