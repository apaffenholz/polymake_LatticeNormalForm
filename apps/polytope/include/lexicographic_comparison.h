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

#ifndef LEXICOGRAPHIC_COMPARISON_H
#define LEXICOGRAPHIC_COMPARISON_H

#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"

namespace polymake { namespace polytope {

namespace {

template <typename E>
int lexicographic_compare(const GenericVector<E>& v1, const GenericVector<E>& v2) {
    //auto it1=entire(v1.top());
    auto it2=entire(v2.top());
    for (auto it1 = entire(v1.top()); !it1.at_end(); ++it1, ++it2) {
        if ( *it1 < *it2 )
            return -1;
        
        if ( *it1 > *it2 )
            return 1;
    }
    return 0;
}

template <typename E>
int lexicographic_compare(const Matrix<E>& m1, const Matrix<E>& m2) {

    for ( Int i = 0; i < m1.rows(); ++i ) {
        int cmp = lexicographic_compare(m1[i],m2[i]);
        if ( cmp != 0 )
            return cmp;
    }
    return 0;
}

} // end anonymous namespace

} } // end of namespaces polymake::polytope

#endif