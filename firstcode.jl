# This code will find the eigenvalues and eigenvectors of a matrix

import LinearAlgebra.eigen

function matrix_construct(x, y)
    mat = zeros(x, y)
    for i in eachindex(mat)
        println("Input value for index $i")
        input = parse(Float32, readline())
        setindex!(mat, [input], [i])
    end
return mat
end

function size_input()
    println("Enter number of rows:")
    rows = parse(Int32, readline())
    println("Enter number of columns:")
    cols = parse(Int32, readline())
    return (rows, cols)
end

size = size_input()

complete_matrix = matrix_construct(size[1], size[2])
eigenvals = eigen(complete_matrix).values
eigenvectors = eigen(complete_matrix).vectors

println("Given matrix is :\n")
display(complete_matrix)
println("with eigenvalues :\n")
display(eigenvals)
println("and eigenvectors :\n")
display(eigenvectors)