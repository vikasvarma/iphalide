#include "Halide.h"
#include <math.h>

using namespace Halide;

/**
 * @brief 2D Inverse Discrete Cosine Transformation
 * 
 * TODO
 * 
 */
class idct2d : public Generator<idct2d>
{
    public:
        Input<Buffer<double>> coeff{"coeff",2};
        Output<Buffer<uint8_t>> image{"image",2};

        void generate()
        {   
            auto M = coeff.width();
            auto N = coeff.height();

            Expr PI (M_PI);
            Expr sqrt2 (sqrt(2));
            Expr sqrtM (sqrt(M));
            Expr sqrtN (sqrt(N));

            RDom r(coeff);
            image(x,y) = cast<uint8_t>(round(sum( 
                select(r.x > 0, sqrt2/sqrtM, cast<double>(1/sqrtM)) *
                select(r.y > 0, sqrt2/sqrtN, cast<double>(1/sqrtN)) *
                cos((PI*r.x*(2*x+1))/(2*M)) *
                cos((PI*r.y*(2*y+1))/(2*N)) *
                coeff(r.x,r.y)
            )));
        }

        void schedule()
        {
            if (auto_schedule){
                coeff.set_estimates({{0,5},{0,4}});
                image.set_estimates({{0,5},{0,4}});
            } else {
                image.compute_root();
            }
        }

    private:
        Func I{"I"}, idct{"idct"}, dct{"dct"};
        Var x{"x"}, y{"y"};
        
};

HALIDE_REGISTER_GENERATOR(idct2d, idct2d_generator)