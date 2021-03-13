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

#ifndef ROW_PERMUTATIONS_H
#define ROW_PERMUTATIONS_H

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

#include <numeric>


//#define DEBUG

namespace polymake { namespace polytope {



class RestrictedPermutation {

    public:
    Array<Int> col_perm;
    std::vector<Int> blocks;

    RestrictedPermutation() {
        col_perm.resize(0);
        blocks.resize(0);
    }

    RestrictedPermutation(const Array<Int>& col_perm_in, const std::vector<Int>& blocks_in) {
        col_perm.resize(col_perm_in.size());
        col_perm = col_perm_in;
        blocks = blocks_in;
    }

    RestrictedPermutation(Int n) {
        col_perm.resize(n);
        std::iota (std::begin(col_perm), std::end(col_perm), 0);

        blocks.resize(2);
        blocks[0] = 0;
        blocks[1] = n;
    }
};

namespace {
    
class RowPermutation {

    public:
    Int row;
    RestrictedPermutation restricted_perm;

    RowPermutation() {}

    RowPermutation(Int r, Array<Int>& col_perm, std::vector<Int>& blocks) {
        row = r;
        restricted_perm = RestrictedPermutation(col_perm, blocks);
    }

    RowPermutation(Int r, RestrictedPermutation& rp) {
        row = r;
        restricted_perm = rp;
    }

    void replace(Int r, Array<Int>& col_perm, std::vector<Int>& blocks) {
        row = r;
        restricted_perm = RestrictedPermutation(col_perm, blocks);
    }
};


    /* 
        given two pairs of a row index and a permutation of the columns, check which one is lex larger

        -1 : first smaller than second
        0: equal
        1: first larger than second
    */
    int compare_row_permutations ( const Matrix<Integer>& m,
                                   const Int x, const Array<Int>& x_perm,
                                   const Int y, const Array<Int>& y_perm ){
        Int n = m.cols();
        for ( Int i = 0; i < n; ++i ) {
            if ( m(x,x_perm[i]) < m(y,y_perm[i]) ) {
                return -1;
            }
            if ( m(x,x_perm[i]) > m(y,y_perm[i]) ) {
                return 1;
            }
        }

        return 0;
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
    RestrictedPermutation lexicographic_maximal_permutation(const Vector<Integer>& v, const RestrictedPermutation& rp) {
        std::vector<Int> new_blocks;
        Array<Int> new_perm(rp.col_perm);

        new_blocks.push_back(0);
        Int n_blocks = rp.blocks.size();
        for(Int i = 0; i < n_blocks-1; ++i ) {
            Int a = rp.blocks[i];
            while ( a < rp.blocks[i+1] ) {
                Integer max = v[new_perm[a]];
                Int max_index = a;
                for ( Int b = a+1; b < rp.blocks[i+1]; ++b) {
                    if ( v[new_perm[b]] > max ) {
                        max = v[new_perm[b]];
                        max_index = b;
                    }
                }
                Int temp = new_perm[a];
                new_perm[a] = new_perm[max_index];
                new_perm[max_index] = temp;
                for ( Int c = max_index+1; c < rp.blocks[i+1]; ++c ) {
                    if ( v[new_perm[c]] == max ) {
                        ++a;
                        temp = new_perm[a];
                        new_perm[a] = new_perm[c];
                        new_perm[c] = temp;
                    }
                }
                ++a;
                new_blocks.push_back(a);
            }
        }

        RestrictedPermutation res(new_perm,new_blocks);
        return res;
    }


    /*
        Determines an order of a subset of rows of a matrix whose columns are subject to a permutation
        such that these rows are in lex max order with the given permutation
        we get a matrix m
        a permutation on the columns
        a subset of the rows

        returns the subset of rows reordered such that lex larger rows come first
    */
    Array<Int> lexicographic_maximal_row_order(const Matrix<Integer>& m, const std::vector<Int>& row_list, const Array<Int>& col_perm ) {
        Int n = row_list.size();
        Array<Int> lexicographic_maximal_row_list(row_list);

        for ( Int i = 0; i < n; ++i ) 
            for ( Int j = 0; j < n-1-i; ++j ) 
                if ( -1 == compare_row_permutations(m, lexicographic_maximal_row_list[j], col_perm, lexicographic_maximal_row_list[j+1], col_perm) ) 
                      std::swap(lexicographic_maximal_row_list[j],lexicographic_maximal_row_list[j+1]);

        return lexicographic_maximal_row_list;
    }

        /* 
        given a matrix and
        two pairs of a subset of the rows and a permutation on the columns (same size subsets)
        check which is lex larger

        return -1 if first is lex smaller than second, 1 if larger, 0 if equal.
    */ 
    int compare_remaining_rows ( const Matrix<Integer>& pairing_matrix, 
                                 const Array<Int>& row_order_a, const Array<Int>& col_perm_a,
                                 const Array<Int>& row_order_b, const Array<Int>& col_perm_b ) {
        
        Int n = row_order_a.size();
        for ( Int i = 0; i < n; ++i ) {
            int cmp = compare_row_permutations(pairing_matrix, row_order_a[i], col_perm_a, row_order_b[i], col_perm_b);
            if ( cmp != 0 ) {
                return cmp;
            }
        }
        return 0;
    }

}

} }

#endif