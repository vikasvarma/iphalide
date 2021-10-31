#include "Halide.h"
#include "util.h"
#include <string.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Halide;
using namespace std;

#define PROFILE = "REGULAR"

struct REGULAR 
{
    struct step1
    {
        int block_size = 8;
        int step_size = 3;
        int search_win = 39;
        int thr = 2500;
        int max_match = 16;
    } STEP1;

    struct step2
    {
        int block_size = 8;
        int step_size = 3;
        int search_win = 39;
        int thr = 400;
        int max_match = 32;
    } STEP2;
};

/**
 * @brief BLOCK MATCHING & 3D DENOISING
 * 
 * TODO
 * 
 */
class bm3d : public Generator<bm3d>
{
    public:
        Input<Buffer<uint8_t>> noisy{"noisy",2};
        Output<Buffer<float>> matches{"matches",2};

        // Image width and height: `W` `H`
        Expr W, H, B;

        void generate()
        {
            // Populate state variables for image dimensions:
            W = noisy.width();
            H = noisy.height();
            B = Expr(4);

            // Edges should be replicated:
            Func I = BoundaryConditions::repeat_edge(noisy);
            Func dct = block_dct2d(I,Expr(4));

            // Hard-threshold the DCT coefficients:
            dct(x,y,p,q) = select(abs(dct(x,y,p,q)) < 0, 0, dct(x,y,p,q));
            
            // Find the matches for a specific block:
            Halide::Runtime::Buffer<float> distance = 
            Halide::Runtime::Buffer<float>(8,8);

            matches(x,y) = distance(x,y);
            matches.trace_stores();
            
            //matches(x,y) = block_matching(dct, x*B, y*B, 8, 4, 10);
        }

        void schedule()
        {
            if (auto_schedule)
            {
                noisy.set_estimates({{0,12},{0,12}});
                matches.set_estimates({{0,3},{0,3}});
            } else {
                matches.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, bx{"bx"}, by{"by"}, z{"z"}, p{"p"}, q{"q"};

        Func block_dct2d(Func I, Expr B)
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
            Func LUT {"LUT"};

            // Generate DCT for each block:
            // TODO - Separate both X,Y coefficient calculation (performance).
            Expr PI (M_PI);
            RDom b(0,B,0,B);
            LUT(x,y,p,q) = Expr(2.0f/B) * sum(
                cos((PI*p*(2*b.x+1))/(Expr(2*B))) *
                cos((PI*q*(2*b.y+1))/(Expr(2*B))) *
                I(x+b.x, y+b.y)
            );

            // Divide first row and col of each coefficient block by sqrt(2):
            RDom br(0,1,0,B);
            RDom bc(0,B,0,1);
            LUT(x,y,br.x,br.y) = LUT(x,y,br.x,br.y)/Expr(sqrt(2));
            LUT(x,y,bc.x,bc.y) = LUT(x,y,bc.x,bc.y)/Expr(sqrt(2));

            return LUT;
        }

        Tuple get_search_window(
            Expr xref,
            Expr yref,
            Expr window_size,
            Expr block_size
        ){
            /** Construct search window reduction domain.
             * @brief Construct a reduction domain that spans the extent of the 
             * search window.
             *      
             *     * * * * * * * * * * xxxx
             *     * * * * * * * * * *    x
             *     * * * * * * * * * *    x
             *     * * * * * @ * * * *    x
             *     * * * * * * * * * *    x
             *     * * * * * * * * * *    x
             *     x                      x
             *     xxxxxxxxxxxxxxxxxxxxxxxx
             * 
             *  * - Top left coordainte of the search block.
             *  @ - Reference block top-left coordaintes.
             *  x - Search extents.
             */

            Expr xmin, ymin, xmax, ymax;
            xmin = xref - floor((window_size-block_size)/Expr(2.0));
            ymin = yref - floor((window_size-block_size)/Expr(2.0));
            xmax = xmin + window_size - 1;
            ymax = ymin + window_size - 1;

            Tuple window = {
                max(0, xmin), max(0, ymin),
                min(W - block_size + 1, xmax), min(H - block_size + 1, ymax)
            };

            return window;
        }

        /**
        Tuple block_matching(
            Func LUT,
            Expr xref,
            Expr yref,
            Expr window_size,
            Expr block_size,
            Expr max_blocks
        ){
            /** Block matching routine for BM3D
             * @brief Search the window space of <xref,yref> block to determine 
             * similar blocks.
             * 
             * Returns the top-left block coordiantes that match with the 
             * reference block.
             

            Func similar{"similar"};
            
            Tuple limits = get_search_window(xref,yref,window_size,block_size);
            RDom block (0,block_size, 0,block_size);

            // Compute distances:
            Func distance(x,y) = 
                sum(
                    pow(
                        LUT(xref,yref,block.x,block.y) - 
                        LUT(
                            window.x,
                            window.y, 
                            block.x,
                            block.y
                        ), 2
                    )
                ) / pow(window_size,2);

            // sorted_dist calls the external C function to return the sorted 
            // argmin order:
            std::vector<ExternFuncArgument> args;
            args.push_back(distance);
            args.push_back(pow(window_size,2));
            args.push_back(max_blocks);
            similar.define_extern("argsort_buffer", args, Float(32), 1);

            return similar;
        }
        */

}; // END - BM3D

HALIDE_REGISTER_GENERATOR(bm3d, bm3d_generator)