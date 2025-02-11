using Plots

function print_coordinates()
    println("Showing example of discrete mesh...")
    height = 9
    for i in Base.OneTo(height)
        if i % 2 == 0
            if i == height
                println("|-h-|-h-|-h-|-h-|-h-|")
            else
                println("k---|---|---|---|---|")
            end
        else
            if i == height
                println("|-h-|-h-|-h-|-h-|-h-|")
            elseif i == height / 3
                println("|---|---|\033[5mj,n\033[0m|---|---|")
            else
                println("|---|---|---|---|---|")
            end
        end
    end
    println("This 2D mesh solves a 1 spatial + 1 temporal dimensional PDE; the j index is the space coordinate, and the n index is the time coordinate.")
end

print_coordinates()

# Let's first consider the bounds of our domain. I'll start with a simple 1D domain. I'll make it a tuple.

x0 = 0.0
xJ = 2.0
x = (x0, xJ)

# We now need to discretize this domain into finite steps of size h.

h = 0.01
xsize = length(range(start = x0, stop = xJ, step = h))
X = [i for i in range(start = x0, stop = xJ, step = h)]

# Now we repeat this for the time domain.

t0 = 0.0
tN = 50.0
t = (t0, tN)




k = 0.01
tsize = length(range(start = t0, stop = tN, step = k))
T = [i for i in range(start = t0, stop = tN, step = k)]

# Note, this method is only stable IFF k^2/h^2 <= 1/2

# We wish to solve the wave equation in 1D, which is U_tt = c^2 * U_xx, with Dirichlet boundary conditions (edges must be a constant value, like 0)
# Using finite differences, we can convert this into a simple recursion algorithm. Let U(x,t) be the function we're solving for. Since we
# discretized our domain into a mesh, then the function is also defined over this mesh; so let U(x,t) = U[x, t], where U is an 
# array with (xJ - x0) / h elements on its first axis and (tN - t0) / k on its second axis.

U = Array{Float64}(undef, tsize, xsize)
# print(eachindex(U))

# We'll also need to define c.

c = 1.0

# These are functions to help iterate through entire axes of an array. Julia doesn't directly allow me to return the axes of an array.

function return_col(matr, n)
    init = Vector{eltype(matr)}(undef, size(matr, 1))
    for i in Base.OneTo(size(matr, 1))
        init[i] = matr[i, n]
    end
    return init
end

function return_row(matr, n)
    init = Vector{eltype(matr)}(undef, size(matr, 2))
    for i in Base.OneTo(size(matr, 2))
        init[i] = matr[n, i]
    end
    return init
end

function return_axis(matr, axisnum, element)
    if axisnum > length(size(matr))
        throw("Requested axis larger than highest axis number in this matrix.")
    else
        init = Vector{eltype(matr)}(undef, size(matr, axisnum))
        for i in Base.OneTo(size(matr, axisnum))
            init[i] = matr[element, i]
        end
        return init
    end
end

# We'll need some initial values (U0) and boundary conditions as well. The boundaries will be 0 (Dirichlet) and the inital value will be 
# some function over space.

function gaussian(x, mean, var, amplitude)
    return amplitude * exp(-(((x - mean) / var)^2))
end

function sinewave(x, frequency, amplitude, phase)
    return amplitude * sin(frequency * x + phase)
end

function step(x, start, stop; amplitude = 1.0)
    if start < x <= stop
        return amplitude
    else
        return 0.0
    end
end

# U0 = 0 .* X .+ 0.2
# U0 = gaussian.(X, 1 / 2, 1/30, 0.7)
# U0 = step.(X, 0.3, 0.4; amplitude=0.2)
U0 = sinewave.(X, 5, 0.1, 0)

function assign_IV(matr, func)
    for i in Base.OneTo(size(matr, 2))
        matr[1, i] = func[i]
    end
end

assign_IV(U, U0)

bc_x0 = 0
bc_xJ = 0

# Now we tell the program how to update each row, which is the same as solving the diffeq. Basically, we're going to code in a way
# to update each space coordinate over time based on the finite difference method. Here, "n" is a step in time and "j" is a step in space.
# We have an r term the defines the ratio of the time & space discretization, and apparently this needs to be below a certain value (1/2 ?).
# Since this is a wave equation, we see that c^2 is included. Also, the second derivative in time leads to a n - 1 index, which will
# throw an error at the start of the simulation. So, we need temporal initial conditions as well.

v = 0

function fin_diff(matr, n, j)
    r = ((k * c) / h)^2
    if j == 1

        return 1/2 * (r * matr[n, j+1] + 2 * (1 - r) * matr[n, j] + r * bc_x0- k^2 * v)

    elseif j == size(matr, 2)

        return 1/2 * (r * bc_xJ + 2 * (1 - r) * matr[n, j] + r * matr[n, j - 1]- k^2 * v)

    else

        return 1/2 * (r * matr[n, j+1] + 2 * (1 - r) * matr[n, j] + r * matr[n, j - 1] - k^2 * v)

    end

end

# We then apply that formula to the entire matrix, which should give us a set of 1D lists that show how our initial gaussian
# evolves in time over each row.

for n in range(1, tsize - 1)
    for j in range(1, xsize)
        U[n + 1, j] = fin_diff(U, n, j)
    end
end

# And finally we make a function to convert this matrix into a set of real points (x, t) based on the actual domain values.
function return_points(matr, row)
    init = [(0.0, 0.0) for i in Base.OneTo(size(matr, 2))]
    for i in Base.OneTo(size(matr, 2))
        init[i] = (X[i], return_row(matr, row)[i])
    end
    return init
end

#Lastly, we animate it.

p1 = plot()
for i in Base.OneTo(tsize)
    global p1 = plot(return_points(U, i), title="Wave Equation via FDM", legend=false, ylim=(- maximum(U0),maximum(U0)), color="red", annotations=((0.7, 0.9), text("Timer : $(round(i * k, digits=2))", :left, 10)))
    display(p1)
    print("\r Loading Animation: $(i) / $(tsize)")
end