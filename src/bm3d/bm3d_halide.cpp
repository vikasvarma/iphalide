/**
 * @file bm3d_halide.cpp
 * @author Vikas Varma (vikas.2110varma@gmail.com)
 * @brief Halide AOT generator modules for BM3D algorithm.
 * @version 0.1
 * @date 2021-10-29
 * 
 * @copyright Copyright (c) 2021
 */

#include "Halide.h"
#include "util.h"
#include <string.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Halide;
using namespace std;

/**-----------------------------------------------------------------------------
 * @brief BLOCKDCT - Generates a LUT of 2D DCT coefficients.
 * 
 * Constructs a LUT of DCT coefficients of all blocks in the 
 * image of specified block size `B` in each spatial dimension 
 * spanning `W` and `H` pixels.
 * 
 * Returns LUT(x,y,p,q) = dct<x,y>(p,q) where <x,y> is the top-left 
 * pixel coordinate of the block and <p,q> is the index of the 
 * (p,q) coefficient in the B-by-B coefficient matrix of the block.
 */
class blockdct : public Generator<blockdct>
{
    public:
        Input<Buffer<uint8_t>> image{"image",2};
        Input<uint8_t> B{"block_size"};
        Output<Buffer<double>> LUT{"LUT",4};

        void generate()
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
            Func I = BoundaryConditions::repeat_edge(image);

            // Generate DCT for each block:
            // TODO - Separate both X,Y coefficient calculation (performance).
            Expr PI(M_PI);
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
        }

        void schedule()
        {
            if (auto_schedule)
            {
                image.set_estimates({{0,12},{0,12}});
                LUT.set_estimates({{0,12},{0,12},{0,4},{0,4}});
            } else {
                LUT.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, p{"p"}, q{"q"};

};

HALIDE_REGISTER_GENERATOR(blockdct, blockdct_generator)

/**-----------------------------------------------------------------------------
 * @brief BLOCKDIST - Identify similar blocks in a search window.
 * 
 * Slides a search window of size `S` in each spatial dimension and computes
 * distances between blocks and the central refrence block in each window to
 * collect similar blocks.
 * 
 * Returns BLOCKS(x,y,n) = blocks<x,y>(n) where <x,y> are coordinate in the 
 * search window super pixel space. <x,y> are <x*S,y*S> in pixel coordiantes.
 * `n` is the index of the matched block to reference in the search window.
 */
class blockdist : public Generator<blockdist>
{
    Input<Buffer<double>> LUT{"LUT",4};
    Input<uint8_t> S{"window_size"};
    Output<Buffer<double>> blocks{"blocks",4};

    public:
        void generate()
        {
            W = LUT.dim(0).extent();
            H = LUT.dim(1).extent();
            B = LUT.dim(2).extent();
            blocks(x,y,m,n) = distance(x,y,m,n);
        }

        void schedule()
        {
            if (auto_schedule)
            {
                LUT.set_estimates({{0,12},{0,12},{0,4},{0,4}});
                blocks.set_estimates({{0,2},{0,2},{0,8},{0,8}});
            } else {
                blocks.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, m{"m"}, n{"n"};
        Expr W, H, B; // image & block extents

        Tuple search_limits(Tuple origin)
        {
            Expr xmin ( max(0,origin[0]) );
            Expr ymin ( max(0,origin[1]) );
            Tuple window = {
                xmin, ymin,
                min(W-1, xmin + S - 1), min(H-1, ymin + S - 1)
            };
            return window;
        }

        Expr distance(
            Expr orgx,
            Expr orgy,
            Expr srcx,
            Expr srcy
        ){
            // Origin of the search window:
            Tuple origin = {orgx * S, orgy * S};
            Expr offset  = ceil(cast<double>(S)/Expr(2.0)) - 
                           ceil(cast<double>(B)/Expr(2.0));

            // Find the center reference block:
            Expr xref = cast<int>(origin[0] + offset);
            Expr yref = cast<int>(origin[1] + offset);

            // Construct a search window extent:
            Tuple slim = search_limits(origin);

            // Compute distances:
            RDom b (0,B,0,B,"block_ind");
            Expr d = sum(pow(
                LUT(xref,yref,b.x,b.y) - 
                LUT(
                    min(slim[0] + srcx, slim[2]),
                    min(slim[1] + srcy, slim[3]), 
                    b.x, b.y
                ), 2
            )) / pow(B,2);

            return d;
        }
};

HALIDE_REGISTER_GENERATOR(blockdist, blockdist_generator)


/**-----------------------------------------------------------------------------
 * @brief ARGSORT - Sorted indices of the block data.
 * 
 */
class argsort : public Generator<argsort>
{
    Input<Buffer<double>> data{"data",4};
    Output<Buffer<uint16_t>> output{"output",4};

    public:
        void generate()
        {
            W = data.dim(0).extent();
            H = data.dim(1).extent();
            S = cast<double>(data.dim(3).extent());
            output(x,y,m,n) = sortmin(data,x,y);
        }

        void schedule()
        {
            if (auto_schedule)
            {
                data.set_estimates({{0,2},{0,2},{0,8},{0,8}});
                output.set_estimates({{0,2},{0,2},{0,8},{0,8}});
            } else {
                output.compute_root();
            }
        }

    private:
        Var x{"x"}, y{"y"}, m{"m"}, n{"n"};
        Expr W,H,S;

        Tuple sortmin(
            Func data,
            Expr xind,
            Expr yind
        ){
            // Extract the data from specified <x,y> indices:
            Func values{"values"};
            values(x,y) = data(xind,yind,x,y);

            // Call external C function to return the sorted argmin order:
            Func sorted{"sorted"};
            std::vector<ExternFuncArgument> args;
            args.push_back(values);
            sorted.define_extern("argsort_buffer", args, Float(64), 1);

            // Hard-coded tuples:
            //sorted(xind,yind,m,n) = sorted( (n * cast<int>(S)) + m );

            return Tuple(sorted(x));
        }
};

HALIDE_REGISTER_GENERATOR(argsort, argsort_generator)