#include "stdint.h"
#include "stdio.h"
#include "scotch.h"

int main() {
    int ver = SCOTCH_VERSION;
    printf("%d", ver);
    return 0;
}
