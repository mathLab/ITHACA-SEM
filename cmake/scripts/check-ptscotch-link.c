#include "stdint.h"
#include "stdio.h"
#include "scotch.h"

int main() {
    // Include a simple call to SCOTCH_graphAlloc to test for linkage on macOS.
    SCOTCH_Graph *graph = SCOTCH_graphAlloc();
    SCOTCH_graphFree(graph);

    return 0;
}
