#include "Halide.h"
#include <math.h>

using namespace Halide;

/**
 * @brief 2D Discrete Cosine Transformation
 * 
 * TODO
 * 
 */
class dct2d : public Generator<dct2d>
{
    public:
        Input<Buffer<uint8_t>> image{"image",2};
        Output<Buffer<double>> coeff{"coeff",2};

        void generate()
        {
            auto M = image.width();
            auto N = image.height();
            Expr PI (M_PI);
            Expr sqrt2 (sqrt(2));
            Expr sqrtM (sqrt(M));
            Expr sqrtN (sqrt(N));

            I(x,y) = cast<double>(image(x,y));
            RDom r(image);
            dct(x,y) = sum(
                cos((PI*x*(2*r.x+1))/(2*M)) *
                cos((PI*y*(2*r.y+1))/(2*N)) * 
                I(r.x,r.y)
            );

            // Divide first row and col by sqrt(2):
            RDom row(0,1,0,N);
            RDom col(0,M,0,1);
            dct(row.x,row.y) = dct(row.x,row.y)/sqrt2;
            dct(col.x,col.y) = dct(col.x,col.y)/sqrt2;
            
            // Compute final dct:
            coeff(x,y) = (sqrt2/sqrtM) * (sqrt2/sqrtN) * dct(x,y);
        }

        void schedule()
        {
            if (auto_schedule){
                image.set_estimates({{0,5},{0,4}});
                coeff.set_estimates({{0,5},{0,4}});
            } else {
                coeff.compute_root();
            }
        }

    private:
        Func I{"I"}, dct{"dct"};
        Var x{"x"}, y{"y"}, a{"a"}, b{"b"};
};

HALIDE_REGISTER_GENERATOR(dct2d, dct2d_generator)