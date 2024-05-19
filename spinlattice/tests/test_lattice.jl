using Test
using Mocking

# Include the module containing the `get_lattice_data` function
test_folder = @__DIR__
parent_folder = dirname(test_folder)
include(joinpath(parent_folder, "lattice.jl"))

# result = get_lattice_data("dimelar-hexagonal-lattice", "dimer-hexagonal", [2, 2], "periodic")


@testset verbose = false "TestGetLatticeData" begin

    @testset "2x2_periodic" begin
        expected_output = (
            bonds=[
                [0, 1], [1, 2], [2, 3], [2, 9], [3, 4], [3, 12], [4, 5], [5, 6],
                [6, 7], [6, 13], [7, 0], [7, 8], [8, 9], [9, 10], [10, 11], [10, 1],
                [11, 12], [11, 4], [12, 13], [13, 14], [14, 15], [14, 5], [15, 8], [15, 0]
            ],
            bond_types=[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            coordinates=[
                [0.0000000000, 0.0000000000], [0.3333333333, 0.0000000000], [0.5000000000, 0.2886751346], [0.8333333333, 0.2886751346],
                [1.0000000000, 0.0000000000], [1.3333333333, 0.0000000000], [1.5000000000, 0.2886751346], [1.8333333333, 0.2886751346],
                [0.0000000000, 0.5773502692], [0.3333333333, 0.5773502692], [0.5000000000, 0.8660254038], [0.8333333333, 0.8660254038],
                [1.0000000000, 0.5773502692], [1.3333333333, 0.5773502692], [1.5000000000, 0.8660254038], [1.8333333333, 0.8660254038]
            ]
        )

        result = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [2, 2], "periodic")
        @test result == expected_output
    end

    @testset "1x1_open" begin
        expected_output = (
            bonds=[[0, 1], [1, 2], [2, 3]],
            bond_types=[1, 0, 0],
            coordinates=[
                [0.0000000000, 0.0000000000], [0.3333333333, 0.0000000000], [0.5000000000, 0.2886751346], [0.8333333333, 0.2886751346]
            ]
        )
        result = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [1, 1], "open")
        @test result !== nothing
        @test result == expected_output
    end

    @testset "100x100_periodic" begin
        result = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [100, 100], "periodic")
        @test result !== nothing
        @test length(result.bonds) == 100 * 100 * 6
    end

    @testset "test_name_with_space" begin
        result_with_space = get_lattice_data("dimer hexagonal lattice", "dimer hexagonal", [10, 10], "periodic")
        result_with_hyphen = get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [10, 10], "periodic")
        @test result_with_space == result_with_hyphen
    end

    @testset "Exception: invalid_boundary" begin
        @test_throws ArgumentError get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [1, 1], "invalid_boundary")
    end

    @testset "Exception: incorrect_size_list_length" begin
        @test_throws ArgumentError get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", [0], "open")
    end

    @testset "Exception: size_list_contains_non_integer" begin
        @test_throws MethodError get_lattice_data("dimer-hexagonal-lattice", "dimer-hexagonal", Any[1, "a"], "open") #No method matching get_lattice_data(::String, ::String, ::Any[], ::String)
    end

    @testset "Exception: invalid_lattice_name" begin
        @test_throws ArgumentError get_lattice_data("invalid_lattice_name", "dimer-hexagonal", [1, 1], "periodic")
    end

    @testset "Exception: invalid_cell_name" begin
        @test_throws ArgumentError get_lattice_data("dimer-hexagonal-lattice", "invalid_cell_name", [1, 1], "open")
    end
end
