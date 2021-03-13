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
#include "polymake/integer_linalg.h"

#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/common/lattice_tools.h"
#include "polymake/graph/compare.h"

#include "polymake/polytope/row_permutations.h"
#include "polymake/polytope/lexicographic_comparison.h"
#include "polymake/polytope/normal_form.h"

#include <numeric>
#include <algorithm>

//#define DEBUG

namespace polymake { namespace polytope {

ListReturn affine_normal_form_from_vertices_and_pairing_matrix(const Matrix<Integer>& v, const Matrix<Integer>& pm,  bool apply_palp_perm = false ) {
    AffineNormalFormTransformation ANFT = affine_normal_form_from_vertices_pairing(v,pm,apply_palp_perm);
    ListReturn ret;
    ret << ANFT.normal_form;
    ret << ANFT.transformation;
    ret << ANFT.translation;
    return ret;
}

ListReturn normal_form_from_vertices_and_pairing_matrix(const Matrix<Integer>& v, const Matrix<Integer>& pm,  bool apply_palp_perm = false ) {
    NormalFormTransformation NFT = normal_form_from_vertices_pairing(v,pm,apply_palp_perm);
    ListReturn ret;
    ret << NFT.normal_form;
    ret << NFT.transformation;
    return ret;
}

Function4perl(&affine_normal_form_from_vertices_and_pairing_matrix, "affine_normal_form_from_vertices_and_pairing_matrix");
Function4perl(&normal_form_from_vertices_and_pairing_matrix, "normal_form_from_vertices_and_pairing_matrix");
} }