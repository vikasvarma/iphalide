#include "Halide.h"
#include "util.h"
#include <string.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Halide;
using namespace std;

/**
 * @brief BLOCK MATCHING & 3D DENOISING
 * 
 * TODO
 * 
 */
class bm3d : public Generator<bm3d>
{
    public:
        Input<Buffer<uint8_t>> image{"image",2};
        Input<uint8_t> B{"block_size"};
        Input<uint8_t> S{"block_stride"};
        Input<uint8_t> W{"search_window"};
        Output<Buffer<double>> output{"output",4};
        Output<Buffer<uint16_t>> order{"order",4};

        // Image width and height: `W` `H`
        Expr M, N;

        void generate()
        {
            // Populate state variables for image dimensions:
            M = image.dim(0).extent();
            N = image.dim(1).extent();

            // Replicate edges:
            Func img = BoundaryConditions::repeat_edge(image);

            // Compute 2D dct block lookup table:
            Func LUT = blockdct(img);
            Var wx{"wx"}, wy{"Wy"};
            
            //Compute block distances:
            Func scores{"scores"};
            RDom b(0,B,0,B,"block");
            scores(x,y,wx,wy) = sum(b, 
                pow(
                    LUT(x*S, y*S, b.x, b.y) - 
                    LUT(x*S + wx, y*S + wy, b.x, b.y),
                    2
                )
            )/pow(B,2);
            
            // Sort buffer and return sort order:
            Func sorted{"sorted"};
            std::vector<ExternFuncArgument> args;
            args.push_back(scores);
            sorted.define_extern("argsort", {scores}, Float(64), 4);

            // Realize output:
            output(x,y,wx,wy) = scores(x,y,wx,wy);
            order(x,y,wx,wy) = cast<uint16_t>(sorted(x,y,wx,wy));
        }

        void schedule()
        {
            if (auto_schedule) {
                image.set_estimates({{0,12},{0,12}});
                output.set_estimates({{0,4},{0,4},{0,6},{0,6}});
                order.set_estimates({{0,4},{0,4},{0,6},{0,6}});
            } else {
                output.compute_root();
                order.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, z{"z"};

        Func dct2d(Func I)
        {
            /**
             * @brief Construct a LUT of DCT coefficients of all blocks in the 
             * image of specified block size `B` in each spatial dimension 
             * spanning `W` and `H` pixels.
             * 
             * Returns LUT(x,y,p,q) = dct<x,y>(p,q) where <x,y> is the top-left 
             * pixel coordinate of the block and <p,q> is the index of the 
             * (p,q) coefficient in the B-by-B coefficient matrix of the block.
             */

            // TODO - Separate both X,Y coefficient calculation (performance).
            Func LUT {"LUT"};
            Expr PI (M_PI);
            RDom r(0,B,0,B,"block");
            Var bx{"bx"}, by{"by"};
            LUT(x,y,bx,by) = select(bx==0 || by==0, Expr(1/sqrt(2)), 1) *
                             Expr(2.0f/B) * 
                             sum(
                                cos((PI*bx*(2*r.x+1))/(Expr(2*B))) *
                                cos((PI*by*(2*r.y+1))/(Expr(2*B))) *
                                I(x+r.x, y+r.y)
                             );
            return LUT;
        }

        Func blockdct(Func I)
        {
            // Generate DCT for each block:
            // TODO - Separate both X,Y coefficient calculation (performance).
            Func LUT {"LUT"};
            Expr PI(M_PI);
            RDom b(0,B,0,B);
            Var bx{"bx"}, by{"by"};
            LUT(x,y,bx,by) = Expr(2.0f/B) * sum(
                cos((PI*bx*(2*b.x+1))/(Expr(2*B))) *
                cos((PI*by*(2*b.y+1))/(Expr(2*B))) *
                I(x+b.x, y+b.y)
            );

            // Divide first row and col of each coefficient block by sqrt(2):
            RDom br(0,1,0,B);
            RDom bc(0,B,0,1);
            LUT(x,y,br.x,br.y) = LUT(x,y,br.x,br.y)/Expr(sqrt(2));
            LUT(x,y,bc.x,bc.y) = LUT(x,y,bc.x,bc.y)/Expr(sqrt(2));

            return LUT;
        }


}; // END - BM3D

HALIDE_REGISTER_GENERATOR(bm3d, bm3d_generator)