#include "Halide.h"

// Use an extern stage to do a sort
extern "C" int argsort(
    halide_buffer_t* data,
    halide_buffer_t* order
);