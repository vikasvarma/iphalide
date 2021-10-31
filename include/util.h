#include "Halide.h"

// Use an extern stage to do a sort
extern "C" int argsort_buffer(
    halide_buffer_t *in,
    halide_buffer_t *out
);

extern "C" class SearchSpace;