#include "stdint.h"
#include "stdio.h"
#include "scotch.h"

int main() {
    float ver = SCOTCH_VERSION;
    int patch = SCOTCH_PATCHLEVEL;
    printf("%2.1f.%d", ver, patch);
    return 0;
}
