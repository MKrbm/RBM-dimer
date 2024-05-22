using Pkg
Pkg.activate(".")
using PyCall

np = pyimport("numpy")

a = np.array([1, 2, 3])

typeof(a)

using LinearAlgebra, SparseArrays


rows = [1,3,4,2,1,3,1,4,1,5]
cols = [1,1,1,2,3,3,4,4,5,5]
vals = [5,-2,-4,5,-3,-1,-2,-10,7,9]

A = sparse(rows, cols, vals, 5, 5)

using Random

function generate_random_sparse_matrix(m, n, density)
    rows = []
    cols = []
    vals = Float64[]
    for i in 1:(m * n * density)
        push!(rows, rand(1:m))
        push!(cols, rand(1:n))
        push!(vals, rand(0.1 : 0.6) * 20 - 10)
    end
    println("typeof", typeof(vals))
    return sparse(rows, cols, vals, m, n)
end

A = generate_random_sparse_matrix(30, 100, 0.2)
c = randn(30)
qra = qr(A)
cc = A * (qra \ c);
sum((cc-c).^2)


function foo1()
    print("Hello1")
    return false
end

function foo2()
    print("Hello2")
    return true
end

function foo3()
    print("Hello3")
    return true
end

foo1() && foo3() # second function is not called

foo2() || foo3() # second function is not called


using Random    
a = randn(10);

println(a)
map(x -> -x^2, a)
println(a)

f(y) = x -> x^2 + y

g = f(3)
g(4)

println(a)
map!((x, y) -> x^2 - y^2, a, a, a)
println(a)

