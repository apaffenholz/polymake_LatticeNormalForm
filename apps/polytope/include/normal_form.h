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

#ifndef NORMAL_FORM_H
#define NORMAL_FORM_H

#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/permutations.h"
#include "polymake/integer_linalg.h"

#include "polymake/polytope/row_permutations.h"
#include "polymake/polytope/lexicographic_comparison.h"

#include <numeric>
#include <algorithm>

//#define DEBUG

namespace polymake { namespace polytope {


class NormalFormTransformation {
    public:
    Matrix<Integer> normal_form;
    Matrix<Integer> transformation;
};

class AffineNormalFormTransformation {
    public:
    Matrix<Integer> normal_form;
    Matrix<Integer> transformation;
    Vector<Integer> translation;
};

class PartialPermutation {

    public: 
    RestrictedPermutation restricted_perm;
    std::vector<Int> row_perm;
    std::vector<Int> remaining_rows;

    PartialPermutation() {
        restricted_perm = RestrictedPermutation();
        row_perm.resize(0);
        remaining_rows.resize(0);
    }

    PartialPermutation(Int n, Int f) {
        restricted_perm = RestrictedPermutation(n);
        row_perm.resize(0);
        remaining_rows.resize(f);
        std::iota (std::begin(remaining_rows), std::end(remaining_rows), 0);
    } 

  friend std::ostream& operator<< (std::ostream & os, const PartialPermutation& pp) {
    os << "col_perm:  ";
    for ( const auto& c: pp.restricted_perm.col_perm ) {
        os << c << " ";
    } 
    os << endl;
    os << "row_perm:  ";
    for ( const auto& c: pp.row_perm ) {
        os << c << " ";
    } 
    os << endl << "blocks:   ";
    for ( const auto& c: pp.restricted_perm.blocks) {
        os << c << " ";
    }
    os << endl;
    return os;
  }
};


/* 
    Given a partial permutation of the pairing matrix
    compute a list of all rows that we can add to the permutation, together with the appropriate column permutations and blocks
*/
Array<RowPermutation> lexicographic_maximal_rows(const Matrix<Integer>& pm, const PartialPermutation& partial_permutation) {

    // temp variables the keep a representative of the currently best permutation
    RowPermutation current_max_perm;
    std::vector<RowPermutation> list_of_max_rows;

    // the current restricted permutation
    // all new permutations must be a refinement of this one
    RestrictedPermutation restricted_perm = partial_permutation.restricted_perm;

    // test all remaining rows
    for ( const auto& i: partial_permutation.remaining_rows ) {

        // bring it into lex max order with the given restriction
        RestrictedPermutation new_partial_permutation = lexicographic_maximal_permutation(pm[i], restricted_perm);

        // check whether it is 
        // - better than the current, then empty the list and put this one in
        // - worse than the current, then discard it
        // - leads to the same partial permutation, then add it to the list
        int replace = list_of_max_rows.empty() ? 1 : compare_row_permutations(pm, i, new_partial_permutation.col_perm, current_max_perm.row, current_max_perm.restricted_perm.col_perm );

        if ( replace == 1 ) {
            current_max_perm.restricted_perm = std::move(new_partial_permutation);
            current_max_perm.row = i;
            list_of_max_rows.clear();
            list_of_max_rows.push_back(current_max_perm);
        }

        if ( replace == 0 ) {
            RowPermutation rp(i,new_partial_permutation);
            list_of_max_rows.push_back(rp);
        }
    }

    return Array<RowPermutation>(list_of_max_rows);
}

/* 
    Given a partial permutation of the pairing matrix
    compute all partial permutations with one additional row
*/
std::vector<PartialPermutation> next_row(const Matrix<Integer>& pm, const PartialPermutation& partial_permutation) {

    Array<RowPermutation> next_rows = lexicographic_maximal_rows(pm,partial_permutation);

    std::vector<PartialPermutation> list_extended_partial_permutations;
    for (const auto& l: next_rows ) {
        // create new partial permutation
        PartialPermutation pp;
        
        // move the new restricted permutation
        pp.restricted_perm = std::move(l.restricted_perm);
        
        // get the old row permutatation and add the new row
        pp.row_perm = partial_permutation.row_perm;
        pp.row_perm.push_back(l.row);

        // remove the row from the set of rows that we still need to consider
        pp.remaining_rows = partial_permutation.remaining_rows;
        pp.remaining_rows.erase(std::remove(pp.remaining_rows.begin(), pp.remaining_rows.end(), l.row), pp.remaining_rows.end());

        // and add the partial permutation to the list
        list_extended_partial_permutations.push_back(pp);
    }

    return list_extended_partial_permutations;
}

/*
    compute a list of all column permutations that lead to a maximal paring matrix
*/
Array<Array<Int> > all_max_permutations(const Matrix<Integer>& pm) {

    Int n = pm.cols();  // number of vertices
    Int f = pm.rows(); // number of facets;

    // initialize the first permutation: identity permutation and just one block
    std::vector<PartialPermutation> list_of_permutations(1);
    PartialPermutation initial_perm(n,f);
    list_of_permutations[0] = initial_perm;

    // keep track of the number of rows we still need to work on
    Int n_remaining_rows = f;
    // if the blocksize is n, then the column permutation is fixed
    // and we can only reorder the rows
    Int blocksize = 1;

    for ( Int i = 0; i < f; ++i ) {

        // no further column permutations are possible
        // we fill the rest of the row permutation below with a simpler method
        // as this is just reordering the rows
        if ( blocksize == n ) {
            break;
        }
        --n_remaining_rows;

        std::vector<PartialPermutation> new_list_of_permutations;
        RowPermutation current_max_perm;

        for (const auto& partial_permutation: list_of_permutations) {

            std::vector<PartialPermutation> row_extensions = next_row(pm, partial_permutation);
            Int n_row_extensions = row_extensions.size();

            if ( n_row_extensions == 0 ) {
                continue;
            }

            // all computed permutations lead to the same partial permutation (but the remaining rows, the column permutation and the blocks leading to it may differ)
            // hence, we can decide using just the first one in the list of row_extensions
            // whether we add, replace, or discard
            PartialPermutation l = row_extensions[0];
            int replace = new_list_of_permutations.empty() ? 1 : compare_row_permutations(pm, l.row_perm.back(), l.restricted_perm.col_perm, current_max_perm.row, current_max_perm.restricted_perm.col_perm );

            // discard
            if ( replace == -1 ) {
                continue;
            }

            // replace
            // here we just make the list empty again
            if ( replace == 1 ) {
                current_max_perm.replace(l.row_perm.back(), l.restricted_perm.col_perm, l.restricted_perm.blocks);
                blocksize = l.restricted_perm.blocks.size()-1;
                new_list_of_permutations.clear();
            }
            // now we insert all row_extensions
            new_list_of_permutations.insert(std::end(new_list_of_permutations), std::begin(row_extensions), std::end(row_extensions));
        }
        // move the newly constructed list to our global store for the next interation
        list_of_permutations = std::move(new_list_of_permutations);
    }

    // if the column permutation is fixed but there are rows left
    // we just order the rows in lex max order
    // we need to do this although we don't care about the row permutation in the end
    // as not all sets of remaining rows for the partial permutations we have computed so far
    // necessaryly lead to the same lex max order
    if ( n_remaining_rows > 0 ) {
        std::vector<PartialPermutation> new_list_of_permutations;

        Array<Int> comparison_row_perm  = list_of_permutations[0].restricted_perm.col_perm;
        Array<Int> comparison_row_order = lexicographic_maximal_row_order(pm, list_of_permutations[0].remaining_rows, comparison_row_perm);

        new_list_of_permutations.push_back(list_of_permutations[0]);

        for ( unsigned int i = 1; i < list_of_permutations.size(); ++i ) {
            Array<Int> rem_row_order = lexicographic_maximal_row_order(pm, list_of_permutations[i].remaining_rows, list_of_permutations[i].restricted_perm.col_perm );
            Int cmp = compare_remaining_rows(pm, comparison_row_order, comparison_row_perm, rem_row_order, list_of_permutations[i].restricted_perm.col_perm );

            if ( cmp == - 1) {
                continue;
            }
            if ( cmp == 1 ) {
                new_list_of_permutations.clear();
                comparison_row_order = rem_row_order;
                comparison_row_perm  = list_of_permutations[i].restricted_perm.col_perm;
            }
            new_list_of_permutations.push_back(list_of_permutations[i]);
        }
        list_of_permutations = std::move(new_list_of_permutations);
    }

    Array<Array<Int> > perms(list_of_permutations.size());
    Int i = 0;
    for ( const auto& partial_permutation: list_of_permutations ) {
        perms[i++] = partial_permutation.restricted_perm.col_perm;
    }

    return perms;
}

Array<Int> palp_permutation(const Matrix<Integer>& pm ) {

    Array<Integer> col_max(pm.cols());
    Array<Integer> col_sum(pm.cols());
    for ( Int i = 0; i < pm.cols(); ++i ) {
        col_sum[i] = accumulate(pm.col(i),operations::add());
        col_max[i] = accumulate(pm.col(i),operations::max());
    }

    Array<Int> col_perm(pm.cols());
    std::iota (std::begin(col_perm), std::end(col_perm), 0);
    for ( Int i = 0; i < pm.cols(); ++i ) {
        Int k = i;
        for ( Int j = i+1; j < pm.cols(); ++j ) {
            if ( col_max[col_perm[j]] < col_max[col_perm[k]] 
                   || col_max[col_perm[j]] == col_max[col_perm[k]] && col_sum[col_perm[j]] < col_sum[col_perm[k]] ) {
                k = j;
            }
        }
        if ( k != i ) {
            Int temp = col_perm[i];
            col_perm[i] = col_perm[k];
            col_perm[k] = temp;
        }
    }    
    return col_perm;
}

NormalFormTransformation normal_form_from_vertices_pairing(const Matrix<Integer>& v, const Matrix<Integer>& pm,  bool apply_palp_perm = false ) {

    const Array<Array<Int> > vperms = all_max_permutations(pm);
    Matrix<Integer> hnf_min(0,0);
    Matrix<Integer> min_transformation;

    Array<Int> palp_col_perm;
    if ( apply_palp_perm ) {
        Array<Int> one_perm = vperms[0];
        palp_col_perm = palp_permutation(permuted_cols(pm,one_perm));
    }

    for ( const auto& perm: vperms ) {

        Matrix<Integer> v_permuted(permuted_rows(v,perm));

        if ( apply_palp_perm ) {
            v_permuted = permuted_rows(v_permuted,palp_col_perm);
        }

        HermiteNormalForm<Integer> hnf_with_companions = hermite_normal_form(v_permuted);
        Matrix<Integer> hnf(hnf_with_companions.hnf);
        if ( !hnf_min.rows() || lexicographic_compare(hnf_min, hnf) == 1 ) {
            hnf_min = hnf;
            min_transformation = std::move(hnf_with_companions.companion);
        }
    }
    
    NormalFormTransformation ret;
    ret.normal_form = std::move(hnf_min);
    ret.transformation = std::move(min_transformation);
    return ret;
}

AffineNormalFormTransformation affine_normal_form_from_vertices_pairing(const Matrix<Integer>& v, const Matrix<Integer>& pm,  bool apply_palp_perm = false ) {

    // the main step: get all permutations of the pairing matrix leading to the maximal one
    const Array<Array<Int> > vperms = all_max_permutations(pm);

    // these variables will hold the minimal HNF and the transformation leading to it
    Matrix<Integer> hnf_min(0,0);
    Matrix<Integer> min_transformation;

    // we do affine normal form, so one vertex is shifted into the origin
    // r will in the end hold its index
    Int r = -1;

    // in the original palp version there was an additional permutation of the pairing matrix in the end
    // if apply_palp_perm is true, then we compute it and use it later
    Array<Int> palp_col_perm;
    if ( apply_palp_perm ) {
        Array<Int> one_perm = vperms[0];
        palp_col_perm = palp_permutation(permuted_cols(pm,one_perm));
    }

    // move one vertex into the origin
    for ( Int v_index = 0; v_index < v.rows(); ++v_index ) {
        Matrix<Integer> v_shifted = v - repeat_row(v[v_index],v.rows());

        // now loop over all permutations
        for ( const auto& perm: vperms ) {

            Matrix<Integer> v_permuted = permuted_rows(v_shifted,perm);
            if ( apply_palp_perm ) {
                v_permuted = permuted_rows(v_permuted,palp_col_perm);
            }

            HermiteNormalForm<Integer> hnf_with_companions = hermite_normal_form(v_permuted);
            Matrix<Integer> hnf(hnf_with_companions.hnf);

            // we have the HNF for this permutation
            // if it is the first or smaller than the current min, then
            // we replace the global vars for the min
            if ( !hnf_min.rows() || lexicographic_compare(hnf_min, hnf) == 1 ) {
                hnf_min = std::move(hnf);
                min_transformation = std::move(hnf_with_companions.companion);
                r = v_index;
            }
        }
    }

    // we return the triple of HNF, transformation, and index of vertex shifted into the origin
    AffineNormalFormTransformation ret;
    ret.normal_form = std::move(hnf_min);
    ret.transformation = std::move(min_transformation);
    ret.translation = -v[r];

    return ret;
}

} }

#endif