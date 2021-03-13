/* Copyright (c) 2020-2021
   Andreas Paffenholz
   Technische Universität Berlin, Germany
   https://www2.mathematik.tu-darmstadt.de/~paffenholz

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/permutations.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"

#include "polymake/polytope/row_permutations.h"

namespace polymake { namespace polytope {


    /*
        Determines an order of a subset of rows of a matrix whose columns are subject to a permutation
        such that these rows are in lex max order with the given permutation
        we get a matrix m
        a permutation on the columns
        a subset of the rows

        returns the subset of rows reordered such that lex larger rows come first
    */
    Array<Int> lexicographic_maximal_row_order(const Matrix<Integer>& m, const Array<Int>& row_list, const Array<Int>& col_perm ) {

        std::vector<Int> row_list_converted;
        for ( const auto& r: row_list ) {
            row_list_converted.push_back(r);
        }
        return lexicographic_maximal_row_order(m, row_list_converted, col_perm);
    }


    /*
        given a vector v, a permutation perm of its entries 
        and a subdivision of the indices (after application of the permutation)
        within which we may further permute the entries 
        (given as list of start indices of each block, and a final entry with the total length, 
        so that block ranges are [block[i]..block[i+1]-1])

        return a new permutation and a new set of blocks such that those produce a lex max vector 
        with the given restrictions
    */
    ListReturn lexicographic_maximal_permutation(const Vector<Integer>& v, const Array<Int>& perm, const Array<Int>& blocks) {

        std::vector<Int> blocks_as_vector(blocks.begin(), blocks.end());
        RestrictedPermutation rp_in(perm,blocks_as_vector);
        RestrictedPermutation rp = lexicographic_maximal_permutation(v,rp_in);

        ListReturn res;
        res << rp.col_perm;
        res << rp.blocks;
        return res;
    }

Function4perl(&lexicographic_maximal_permutation, "lexicographic_maximal_permutation");
Function4perl(&compare_remaining_rows, "compare_remaining_rows");
Function4perl(&lexicographic_maximal_row_order, "lexicographic_maximal_row_order");
Function4perl(&compare_row_permutations, "compare_row_permutations");

} }