
compare_object('1_OK', polytope_in_normal_form(load('1'),apply_palp_permutation=>false) );
compare_object('2_OK', polytope_in_affine_normal_form(load('1'),apply_palp_permutation=>false) );
compare_object('3_OK', polytope_in_affine_normal_form(load('3'),apply_palp_permutation=>false) );
check_boolean( '4', lattice_equivalent(load('4-1'), load('4-2')) );