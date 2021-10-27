#include "util.h"
#include "Halide.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

extern "C" int argsort_buffer(
    halide_buffer_t *in,
    halide_buffer_t *out
){
    if (in->is_bounds_query()) {
        in->dim[0].min = out->dim[0].min;
        in->dim[0].extent = out->dim[0].extent/2;

    } else {
        float *in_start = (float *)in->host;
        float *out_start = (float *)out->host;

        // initialize original index locations:
        std::vector<uint16_t> index(in->dim[0].extent);
        std::iota(index.begin(), index.end(), 0);

        // sort the values:
        std::stable_sort(index.begin(), index.end(), 
            [in_start](size_t a, size_t b) {
                return *(in_start + a) < *(in_start + b);
            }
        );

        // store index and sorted value together:
        for (auto ptr = 0; ptr < in->dim[0].extent; ptr++){
            *(out_start + 2*ptr) = index[ptr];
            *(out_start + 2*ptr + 1) = *(in_start + index[ptr]);
        }
    }
    return 0;
}