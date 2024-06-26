using Test
using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals


@testset "Test suite for print_polynomials" begin
    A = rand(10, 3)

    sets = AVI.construct_SetsOandG(A)

    append!(sets.G_coefficient_vectors, [reshape(   [   -0.77, 
                                                        0.44, 
                                                        1.0, 
                                                        0.0 ], 4, 1)])

    append!(sets.G_coefficient_vectors, [Matrix([       [0.16 -0.11 0.20 -0.19 0.17];
                                                        [-0.98 -0.40 -0.53 0.23 -0.03];
                                                        [-0.07 0.18 -0.48 -0.34 -1.00];
                                                        [1.0 0.0 0.0 0.0 0.0];
                                                        [0.0 1.0 0.0 0.0 0.0];
                                                        [0.0 0.0 1.0 0.0 0.0];
                                                        [0.0 0.0 0.0 1.0 0.0];
                                                        [0.0 0.0 0.0 0.0 1.0]  ])])

    sets.O_terms = [    [0 1 0];
                        [0 0 0];
                        [0 0 1] ]

    sets.leading_terms =  [     [0  2  1  1  0  0];
                                [1  0  1  0  1  0];
                                [0  0  0  1  1  2]  ]

    sets.O_degree_indices = [2]

    polys = [   "x_{2} + 0.44x_{1} - 0.77", 
                "x_{1}^{2} - 0.07x_{3} - 0.98x_{1} + 0.16",
                "x_{1}x_{2} + 0.18x_{3} - 0.4x_{1} - 0.11",
                "x_{1}x_{3} - 0.48x_{3} - 0.53x_{1} + 0.2",
                "x_{2}x_{3} - 0.34x_{3} + 0.23x_{1} - 0.19",
                "x_{3}^{2} - x_{3} - 0.03x_{1} + 0.17"
            ]

    constructed_polys = AVI.print_polynomials(sets; ret=true)

    @test all(polys .== constructed_polys)
end;


@testset "Test suite for convert_term_to_latex" begin
    terms = [   [0 1 2 10 0];
                [0 0 0 7 1];
                [0 1 0 0 0];
                [0 0 1 0 0] ]

    converted_terms = ["", "x_{1}x_{3}", "x_{1}^{2}x_{4}", "x_{1}^{10}x_{2}^{7}", "x_{2}"]

    for i in 1:size(terms, 2)
        term = terms[:, i]
        @test AVI.convert_term_to_latex(term) == converted_terms[i]
    end
end;
