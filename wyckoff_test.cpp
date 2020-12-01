#include "wyckoff.hpp"
//#include "symop_4d.hpp"
// is there a way to include all the headers in the avdv-factorgroup folder??

bool test_wyckoff_construction()
{
    Eigen::Vector3d basis_vector1{0, 0, 1};
    Eigen::Vector3d basis_vector2{0, 0, 0};
    Eigen::Vector3d basis_vector3{1, 0, 0};
    Eigen::Vector3d offset{0.25, 0.25, 0};

    Subspace subspace1 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);
    Eigen::Matrix3d expected_basis_matrix;
    expected_basis_matrix.col(0) = basis_vector1;
    expected_basis_matrix.col(1) = basis_vector2;
    expected_basis_matrix.col(2) = basis_vector3;

    return (expected_basis_matrix == subspace1.basis_col_matrix() && offset == subspace1.offset());
}

Eigen::Matrix3d make_xyz_rotation_matrix(double degrees)
{
    Eigen::AngleAxisd rotation_generator(degrees * M_PI / 180.0, Eigen::Vector3d(1, 1, 1).normalized());
    return rotation_generator.matrix();
}
Eigen::Matrix3d make_y_rotation_matrix(double degrees)
{
    Eigen::AngleAxisd rotation_generator(degrees * M_PI / 180.0, Eigen::Vector3d(0, 1, 0));
    return rotation_generator.matrix();
}

PeriodicGroup make_identity_group(double tol)
{
    Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    PeriodicGroup identity_group({identity_matrix}, symop_comparator, symop_multiplier);
    return identity_group;
}

PeriodicGroup make_rotation_60_group(double tol)
{
    Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    PeriodicGroup rotation_60_group({make_z_rotation_matrix(60)}, symop_comparator, symop_multiplier);
    return rotation_60_group;
}
PeriodicGroup make_rotation_group(double degrees, double tol)
{
    Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    PeriodicGroup rotation_60_group({make_z_rotation_matrix(degrees)}, symop_comparator, symop_multiplier);
    return rotation_60_group;
}
PeriodicGroup make_rotation_60_group_diagonal(double tol)
{
    Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    Eigen::Matrix3d rotation_generator_matrix = make_xyz_rotation_matrix(60);
    PeriodicGroup rotation_60_group({rotation_generator_matrix}, symop_comparator, symop_multiplier);
    return rotation_60_group;
}

bool test_find_coset(double tol)
{
    /*    make a rotation group
     *    make identity group
     *    find coset of identity group in rotations group
     *    should be length rotation group -1, and not contain identity, and contain all other elements of rotation
     * group.*/

    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    PeriodicGroup identity_group = make_identity_group(tol);
    std::vector<SymOp> coset = find_coset(rotation_60_group, identity_group);

    // check if subgroup is not in coset:

    BinarySymOpPeriodicCompare_f symop_comparator = rotation_60_group.get_comparator();
    auto identity_compare = [identity_group, symop_comparator](SymOp input_symop) {
        for (SymOp group_op : identity_group.operations())
        {
            if (symop_comparator(group_op, input_symop))
            {
                return true;
            }
        }
        return false;
    };
    if (std::find_if(coset.begin(), coset.end(), identity_compare) != coset.end())
    {
        return false;
    }

    // Check if every coset element is in larger group(factor group)
    for (SymOp coset_symop : coset)
    {
        auto coset_compare = [coset_symop, symop_comparator](SymOp input_symop) {
            return symop_comparator(coset_symop, input_symop);
        };
        if (std::find_if(rotation_60_group.operations().begin(), rotation_60_group.operations().end(), coset_compare) ==
            rotation_60_group.operations().end())
        {
            return false;
        }
    }

    return true;
}

bool test_make_SymOp4dMatrix()
{
    // make symop
     // make symop 4d, with translation
     // check all elements
     //
    Eigen::Matrix3d test_matrix;
    test_matrix << 1, 0, 0, 0, -1, 0, 0, 0, 1;
    Eigen::Vector3d test_translation;
    test_translation << .25, .25, .25;
    SymOp test_input_symop(test_matrix, test_translation);
    //    Eigen::Matrix4d symop_4d = make_symop_4dmatrix(test_input_symop);
    Eigen::Matrix4d symop_4d = test_input_symop.get_combined_matrix();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (test_matrix(i, j) != symop_4d(i, j))
            {
                return false;
            }
        }
    }
    for (int i = 0; i < 3; i++)
    {
        if (test_translation(i) != symop_4d(i, 3))
        {
            return false;
        }
    }
    for (int i = 0; i < 3; i++)
    {
        if (i < 3 && symop_4d(3, i) != 0)
        {
            return false;
        }
    }
    if (symop_4d(3, 3) != 1)
    {
        return false;
    }

    return true;
}

bool test_make_reynolds_operator(double tol)
{
    /*make roation group
     * find average
     * average should be
     */
    Eigen::Matrix4d expected_average;
    expected_average << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    /*if the translation was not 000, would this move the axis from the 00z?
     */
    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    Eigen::Matrix4d average_matrix = make_4d_reynolds_operator(rotation_60_group);
    return average_matrix.isApprox(expected_average, tol);
}

bool test_make_3d_reynolds_operator(double tol)
{
    Eigen::Matrix3d expected_average;
    expected_average << 0, 0, 0, 0, 0, 0, 0, 0, 1;
    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    Eigen::Matrix3d average_matrix = make_3d_reynolds_operator(rotation_60_group);
    return average_matrix.isApprox(expected_average, tol);
}

bool test_find_invariant_subspace(double tol)
{
    Eigen::Matrix3d expected_basis_col_matrix;
    expected_basis_col_matrix << 0, 0, 0, 0, 0, 0, 0, 0, 1;
    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    Subspace invariant_subspace = find_invariant_subspace(rotation_60_group);
    std::cout << invariant_subspace.formula() << std::endl;
    return expected_basis_col_matrix.isApprox(invariant_subspace.basis_col_matrix());
}

bool test_find_diagonal_axis_invariant_subspace(double tol)
{
    // TODO: can we test a rotation axis group with diagonal axis??

    Eigen::Matrix3d expected_basis_col_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d expected_basis_col = Eigen::Vector3d::Ones().normalized();
    expected_basis_col_matrix.col(0) = expected_basis_col;
    PeriodicGroup rotation_60_group_diagonal = make_rotation_60_group_diagonal(tol);
    Subspace invariant_subspace = find_invariant_subspace(rotation_60_group_diagonal);
    std::cout << invariant_subspace.formula() << std::endl;
    std::cout << invariant_subspace.formula();
    // TODO change comparison so the order of the vectors don't matter??
    return expected_basis_col_matrix.isApprox(invariant_subspace.basis_col_matrix());
}

std::vector<Subspace> load_c6_wyckoff_positions()
{
    // TODO finish to be used below to check rotation group wyckoff positions. two type, the 00z axis and the 6 xyz
    // types
    std::vector<Subspace> wyckoff_positions;
    return wyckoff_positions;
};

bool test_subspace_compare_same(double tol)
{
    // TODO: test subspace comparison function

    Eigen::Vector3d basis_vector1{0, 0, 1};
    Eigen::Vector3d basis_vector2{0, 0, 0};
    Eigen::Vector3d basis_vector3{1, 0, 0};
    Eigen::Vector3d offset{0.25, 0.25, 0};

    Subspace subspace1 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);

    return subspaces_are_equal(subspace1, subspace1, tol);
}

bool test_subspace_compare_different(double tol)
{
    // TODO: test subspace comparison function

    Eigen::Vector3d basis_vector1{0, 0, 1};
    Eigen::Vector3d basis_vector2{0, 0, 0};
    Eigen::Vector3d basis_vector3{1, 0, 0};
    Eigen::Vector3d offset{0.25, 0.25, 0};

    Subspace subspace1 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);
    Subspace subspace2 = Subspace({basis_vector2, basis_vector2, basis_vector3}, offset);

    return !subspaces_are_equal(subspace1, subspace2, tol);
}
bool test_subspace_compare_diff_order(double tol)
{
    // TODO: test subspace comparison function

    Eigen::Vector3d basis_vector1{0, 0, 1};
    Eigen::Vector3d basis_vector2{0, 0, 0};
    Eigen::Vector3d basis_vector3{1, 0, 0};
    Eigen::Vector3d offset{0.25, 0.25, 0};

    Subspace subspace1 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);
    Subspace subspace2 = Subspace({basis_vector2, basis_vector1, basis_vector3}, offset);

    return subspaces_are_equal(subspace1, subspace2, tol);
}
bool test_subspace_compare_equivalent_basis(double tol)
{
    // TODO: test subspace comparison function

    Eigen::Vector3d basis_vector1{0, 0, 1};
    Eigen::Vector3d basis_vector2{0, 0, 0};
    Eigen::Vector3d basis_vector3{1, 0, 0};
    Eigen::Vector3d basis_vector4{0.5, 0, 1};
    Eigen::Vector3d offset{0.25, 0.25, 0};

    Subspace subspace1 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);
    Subspace subspace2 = Subspace({basis_vector1, basis_vector2, basis_vector3}, offset);

    return subspaces_are_equal(subspace1, subspace2, tol);
}

bool test_find_symmetrically_equivalent_wyckoff_positions(double tol)
{
    // test for known subgroup without translation

    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    PeriodicGroup rotation_180_group = make_rotation_group(180, tol);
    PeriodicGroup identity_group = make_identity_group(tol);
    std::vector<SymOp> coset_identity = find_coset(rotation_60_group, identity_group);
    std::vector<SymOp> coset_180 = find_coset(rotation_60_group, rotation_180_group);

    Subspace wyckoff_prototype = find_invariant_subspace(identity_group);
    Subspace wyckoff_prototype_180 = find_invariant_subspace(rotation_180_group);
    std::vector<Subspace> wyckoff_e =
        find_symmetrically_equivalent_wyckoff_positions(coset_identity, wyckoff_prototype);
    std::vector<Subspace> wyckoff_180 =
        find_symmetrically_equivalent_wyckoff_positions(coset_180, wyckoff_prototype_180);

    if (wyckoff_e.size() != rotation_60_group.operations().size())
    {
        return false;
    }
    std::vector<Subspace> expected_wyckoff_positions = load_c6_wyckoff_positions();

    std::cout << "identity group wyckoff positions" << std::endl;
    for (auto position : wyckoff_e)
    {
        std::cout << position.formula() << std::endl;
    }

    std::cout << "size of 180 coset:   " << coset_180.size() << std::endl;
    std::cout << "180 group wyckoff positions" << std::endl;
    for (auto position : wyckoff_180)
    {
        std::cout << position.formula() << std::endl;
    }

    /* TODO: check that they are the list of wycoff positions expected (check against list from bilbao with find if)
     */
    return true;
}

bool test_wyckoff_position_compare(double tol)
{
    /*TODO:Compares vectors of symmetrically equivalent wyckoff positions
     * within tolerance and irrespective of order of positions.
     */
    return false;
}

WyckoffList load_expected_cubic_point_group_wyckoff_positions()
{ // TODO: function to generated expected cubic point group
    std::vector<std::vector<Subspace>> point_group_wyckoff_list;
    return point_group_wyckoff_list;
}
bool test_cubic_point_group_wyckoff_positions(double tol)
{
    // TODO:test for simple factor group (ie cubic)
    /* load simple cubic structure from POSCAR,
     * generate factor group(same as point group)
     * generate wyckoff_positions
     * compare number first to expected number of wyckoff positions
     * compare to expected wycoff positions
     * */
    return false;
}
/* TODO:
 * test for known subgroup with translation
 * test for harder factor group (hexagonal and diamond)
 */
/*
bool test_sympp_4d_construction()
{
    Eigen::Matrix3d cart_matrix;
    cart_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    Eigen::Vector3d translation;
    translation << 0.5, 0.5, 0.5;
    Eigen::Matrix4d eigen_4d_matrix;
    eigen_4d_matrix << 1, 0, 0, 0.5, 0, 1, 0, 0.5, 0, 0, 1, 0.5, 0, 0, 0, 1;
    SymOp sym_op(cart_matrix, translation);

    Symop_4d symop_4d_matrix_0(cart_matrix, translation);
    Symop_4d symop_4d_matrix_1(eigen_4d_matrix);
    Symop_4d symop_4d_matrix_2(sym_op);

    if (symop_4d_matrix_0.get_matrix().isApprox(eigen_4d_matrix) &&
        symop_4d_matrix_1.get_matrix().isApprox(eigen_4d_matrix) &&
        symop_4d_matrix_2.get_matrix().isApprox(eigen_4d_matrix))
    {
        return true;
    }
    std::cout << symop_4d_matrix_0.get_matrix() << std::endl;
    std::cout << symop_4d_matrix_1.get_matrix() << std::endl;
    std::cout << symop_4d_matrix_2.get_matrix() << std::endl;

    return false;
}
*/
bool test_projector_operator()
{
    Eigen::Vector3d a;a<<1,1,1;
    Eigen::Vector3d b; b<<0,0,1;
    return (b==projector_operator(a,b));

}

bool test_convert_subspace_to_vectors(double tol)
{

    Eigen::Vector3d a;a<<1,0,0;
    Eigen::Vector3d b;b<<0,1,0;
    Eigen::Vector3d c;c<<0,0,1;
    Eigen::Vector3d d;d<<0.5,0.5,0;

    Subspace test_subspace({a,b,c},d);
    std::vector<Eigen::Vector3d> expected_vector({a,b,c,d});
    std::vector<Eigen::Vector3d> converted_vector = convert_subspace_to_vector(test_subspace, tol);
//needs like comparing each vector to input vectors a ,b, c, and d
    return false;
}

//bool 

int main()
{
    double tol = 1e-6;
    /* all the tests above*/
    // EXPECT_TRUE(//test_wyckoff_construction, "Construction of Wyckoff class");
    EXPECT_TRUE(test_make_SymOp4dMatrix(), "Convert SymOp to 4d Matrix");
    EXPECT_TRUE(test_find_coset(tol), "Generate Expected Coset");
    EXPECT_TRUE(test_make_reynolds_operator(tol), "Generate 4d Reynolds operator");
    EXPECT_TRUE(test_make_3d_reynolds_operator(tol), "Generate 3d Reynolds operator");
    EXPECT_TRUE(test_subspace_compare_same(tol), "Subspace compare same subspaces works");
    EXPECT_TRUE(test_subspace_compare_different(tol), "Subspace compare different subspaces works");
    EXPECT_TRUE(test_find_invariant_subspace(tol), "Find the expected invariant subspace of subgroup");
    EXPECT_TRUE(test_find_symmetrically_equivalent_wyckoff_positions(tol),
                "Find and print symmetrically equivalent wyckoff positions");
    EXPECT_TRUE(test_find_diagonal_axis_invariant_subspace(tol),
                "Find the expected invariant subspace of subgroup with diagnol rotation axis");
    EXPECT_TRUE(test_wyckoff_position_compare(tol), "Wyckoff Positions Comparison");
    EXPECT_TRUE(test_cubic_point_group_wyckoff_positions(tol), "Generated expected cubic point group");
//    EXPECT_TRUE(test_sympp_4d_construction(), "Test symop4d class");
    return 0;
}
