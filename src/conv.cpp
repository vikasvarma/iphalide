#include "Halide.h"
#include <stdio.h>

using namespace Halide;

class conv : public Generator<conv> 
{
    public:
        Input<Buffer<uint8_t>> input{"input",2};
        Input<Buffer<float>> kernel{"kernel",2};
        Output<Buffer<uint8_t>> output{"output",2};

        void generate()
        {
            // Algorithm:
            Func clamped = BoundaryConditions::repeat_edge(input);
            RDom r(kernel);
            output(x,y) = cast<uint8_t>(
                sum(kernel(r.x,r.y)*clamped(x+r.x,y+r.y))
            );
        }

        void schedule()
        {
            if (auto_schedule) {
                input.set_estimates({{0,4096},{0,3072}});
                kernel.set_estimates({{0,11},{0,11}});
                output.set_estimates({{0,4096},{0,3072}});
            } else {
                output.compute_root();
            }
        }
    
    private:
        Var x{"x"}, y{"y"};
};

HALIDE_REGISTER_GENERATOR(conv, conv_generator)