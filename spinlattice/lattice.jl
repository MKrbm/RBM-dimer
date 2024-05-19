using Printf

function get_lattice_data(lattice_name::String, cell_name::String, size::Vector{Int}, boundary::String)::NamedTuple{(:bonds, :bond_types, :coordinates), Tuple{Vector{Vector{Int}}, Vector{Int}, Vector{Vector{Float64}}}}
    # Check if the size has exactly two integers
    if length(size) != 2
        throw(ArgumentError("Size must be a vector of two integers."))
    end

    # Check if the boundary condition is valid
    if !(boundary in ["periodic", "open"])
        throw(ArgumentError("Boundary must be 'periodic' or 'open'."))
    end

    # If lattice_name or cell_name contains space, replace space with _
    lattice_name = replace(lattice_name, " " => "-")
    cell_name = replace(cell_name, " " => "-")

    # Path to the executable
    script_dir = @__DIR__
    executable_path = joinpath(script_dir, "build", "main")

    if !isfile(executable_path)
        throw(ArgumentError("Executable not found. You must compile the C++ code first under the build directory. Current working directory: " * pwd()))
    end

    args = [
        "--lattice_name", lattice_name,
        "--cell_name", cell_name,
        "--L1", string(size[1]),
        "--L2", string(size[2]),
        "--boundary", boundary
    ]

    cmd = Cmd(`$executable_path $args`)
    io = IOBuffer()
    ie = IOBuffer()
    try
        run(pipeline(cmd, stdout=io, stderr=ie))
    catch e
        throw(ArgumentError(String(take!(ie))))
    end
    output = String(take!(io))

    bonds = Vector{Vector{Int}}()
    bond_types = Vector{Int}()
    coordinates = Vector{Vector{Float64}}()

    lines = split(output, "\n")
    parsing_bonds = false
    parsing_coordinates = false

    for line in lines
        if occursin("Bond Types and Bonds:", line)
            parsing_bonds = true
            parsing_coordinates = false
            continue
        elseif occursin("Coordinates:", line)
            parsing_bonds = false
            parsing_coordinates = true
            continue
        end

        if parsing_bonds
            if occursin("Bond Type:", line)
                bond_type_str = split(line, "Bond Type: ")[2] |> x -> split(x, ", Bonds:")[1]
                try
                    bond_type = parse(Int, bond_type_str)
                    push!(bond_types, bond_type)
                catch
                    error("Failed to convert bond type to integer: $bond_type_str")
                end
                bond_list_str = split(split(line, "{")[2], "}")[1] |> x -> split(x, ", ")
                bond_list = [parse(Int, b) for b in bond_list_str]
                push!(bonds, bond_list)
            end
        elseif parsing_coordinates
            if occursin("{", line) && occursin("}", line)
                coord_list_str = split(split(line, "{")[2], "}")[1] |> x -> split(x, ", ")
                coord_list = [parse(Float64, c) for c in coord_list_str]
                push!(coordinates, coord_list)
            end
        end
    end

    return (bonds=bonds, bond_types=bond_types, coordinates=coordinates)
end

# Example usage:
# lattice_data = get_lattice_data("square lattice", "simple2d", [10, 10], "periodic")
# println(lattice_data)

