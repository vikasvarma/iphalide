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
        Input<uint8_t> MAX{"max_blocks"};
        Output<Buffer<uint16_t>> matches{"matches",3};
        Output<Buffer<double>> estimate{"estimate",5};

        // Image width and height: `W` `H`
        Expr M, N;
        Expr thr;

        void generate()
        {
            // Populate state variables for image dimensions:
            M = image.dim(0).extent();
            N = image.dim(1).extent();

            thr = Expr(2.7*25); // lambda * sigma

            // Replicate edges:
            Func img = BoundaryConditions::repeat_edge(image);

            // Compute 2D dct block lookup table:
            Func LUT = dct2d(img);
            Var wx{"wx"}, wy{"Wy"};
            
            //Compute block distances:
            Func scores{"scores"};
            RDom b(0,B,0,B,"block");
            scores(x,y,wx,wy) = sum(b, 
                pow(
                    LUT(x*S, y*S, b.x, b.y) - 
                    LUT(x*S + wx, y*S + wy, b.x, b.y),
                    2
                ),
                "score-sum"
            )/pow(B,2);
            
            // Sort buffer and return sort order:
            Func sorted{"sorted"};
            std::vector<ExternFuncArgument> args;
            args.push_back(scores);
            sorted.define_extern("argsort", {scores}, UInt(16), 4);

            Var mb,bx,by;

            // Group similar blocks and perform 1-D DCT along z-dim:
            Func blockgroup{"blockgroup"};
            blockgroup(x,y,mb,bx,by) = LUT(
                x*S + sorted(x,y,mb%W,mb/W)%W, 
                y*S + sorted(x,y,mb%W,mb/W)/W,
                bx, by
            );

            // Pick only the top N candidates for z-dim dct:
            RDom rz(0,MAX,"rz");
            RDom rb(0,B,"rb");
            Expr PI (M_PI);
            Func bdct3d{"bdct3d"};
            bdct3d(x,y,mb,bx,by) = 
                select(bx==0, Expr(sqrt(1/2.0)), Expr(1.0)) *
                Expr(2.0f/MAX) * sum(
                    cos((PI*mb*(2*rz.x+1))/(Expr(2*MAX))) *
                    LUT(
                        x*S + sorted(x,y,rz.x%W,rz.x/W)%W,
                        y*S + sorted(x,y,rz.x%W,rz.x/W)/W,
                        bx, by
                    ),
                    "z-dct-sum"
                );
            
            // Hard-threshold the DCT coefficients:
            //bdct3d(x,y,mb,bx,by) = select(
            //    bdct3d(x,y,mb,bx,by) < thr,
            //    bdct3d(x,y,mb,bx,by), 0
            //);

            // Perform 3-D inverse DCT:
            Func idctx, idcty, idctz;
            idctz(x,y,mb,bx,by) = 
                round(sum( 
                    select(rz > 0, 
                        sqrt(Expr(2.0f/MAX)), 
                        sqrt(Expr(1.0f/MAX))
                    ) * 
                    cos((PI*rz*(2*mb+1))/(2*MAX)) * 
                    bdct3d(x,y,rz,bx,by),
                    "idctz-sum"
                ));
            
            idcty(x,y,mb,bx,by) = 
                round(sum( 
                    select(rb > 0, 
                        sqrt(Expr(2.0f/B)), 
                        sqrt(Expr(1.0f/B))
                    ) * 
                    cos((PI*rb*(2*by+1))/(2*B)) * 
                    idctz(x,y,mb,bx,rb),
                    "idcty-sum"
                ));
            
            idctx(x,y,mb,bx,by) = 
                round(sum( 
                    select(rb > 0, 
                        sqrt(Expr(2.0f/B)), 
                        sqrt(Expr(1.0f/B))
                    ) * 
                    cos((PI*rb*(2*bx+1))/(2*B)) * 
                    idcty(x,y,mb,rb,by),
                    "idctx-sum"
                ));

            // Assign output:
            estimate(x,y,mb,bx,by) = idctz(x,y,mb,bx,by);
            matches(x,y,z) = sorted(x,y,z%W,z/W);
        }

        void schedule()
        {
            if (auto_schedule) {
                image.set_estimates({{0,12},{0,12}});
                matches.set_estimates({{0,4},{0,4},{0,36}});
                estimate.set_estimates({{0,4},{0,4},{0,10},{0,4},{0,4}});
            } else {
                matches.compute_root();
                estimate.compute_root();
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
            LUT(x,y,0,0) = LUT(x,y,0,0)/Expr(sqrt(2));
            return LUT;
        }


}; // END - BM3D

HALIDE_REGISTER_GENERATOR(bm3d, bm3d_generator)