# 2D Wave Equation Solver
using Plots

function print_coordinates()
    println("Showing example of discrete mesh...")
    height = 9
    for i in Base.OneTo(height)
        if i % 2 == 0
            if i == height - 1
                println("k---|---|---|---|---|  |⟋")
            elseif i == height - 3
                println("k---|---|---|---|---|  |⟋| g")
            else
                println("k---|---|---|---|---|  |⟋|")
            end
        else
            if i == height
                println("|-h-|-h-|-h-|-h-|-h-|g⟋ ")
            elseif i == height - 2
                println("|---|---|---|---|---|  |⋰|⟋")
            elseif i == height / 3
                println("|---|---|\033[5mnjl\033[0m|---|---|  |⟋|")
            else
                println("|---|---|---|---|---|  |⟋|")
            end
        end
    end
    println("This 3D mesh is used to solve the 2D wave equation.")
    println("The space domains include are discretized by h & g, while the")
    println("time domain is discretized by k. The coorindates of each")
    println("element in the mesh is given as a tuple (n,j,l)")
end

print_coordinates()

#Space domain and discretization

x0 = 0.0
xJ = 1.0
x = (x0, xJ)

y0 = 0.0
yL = 1.0
y = (y0, yL)

h = 0.01
xsize = length(range(start = x0, stop = xJ, step = h))
X = [i for i in range(start = x0, stop = xJ, step = h)]

q = 0.01
ysize = length(range(start = y0, stop = yL, step = q))
Y = [i for i in range(start = y0, stop = yL, step = q)]

# Now we repeat this for the time domain.

t0 = 0.0
tN = 0.025
t = (t0, tN)

k = 0.00001
tsize = length(range(start = t0, stop = tN, step = k))
T = [i for i in range(start = t0, stop = tN, step = k)]


U = Array{Float64}(undef, tsize, xsize, ysize)

# We'll also need to define c.

c = 1.0

r = ((k * c) / (h*q))^2

if r > 0.5
    throw("Warning! Numerical instability likely, r value ($(round(r, digits=2))) is above 0.5.")
end

# These are functions to help iterate through entire axes of an array. Julia doesn't directly allow me to return the axes of an array.

function return_row(matr, n)
    init = Vector{eltype(matr)}(undef, size(matr, 2))
    for i in Base.OneTo(size(matr, 2))
        init[i] = matr[n, i]
    end
    return init
end

# We'll need some initial values (U0) and boundary conditions as well. The boundaries will be 0 (Dirichlet) and the inital value will be 
# some function over space.

xx = [i for i in range(x[1], x[end], step=h)]
yy = [i for i in range(x[1], x[end], step=q)]

function gaussian(x,y)
    return exp(-(10*x)^2 - (10*y)^2)

end

function sine(x, y)
    return sin(x * pi * x[end]) + sin(y* pi * y[end])
end

zz = @. gaussian(xx' - (x[end] + x[1]) / 2, yy - (y[end] + y[1]) / 2) * 5
# zz = @. sine(xx', yy)

function edge_set(matr)
    for i in size(matr, 2)
        matr[:, 1, :] = zeros(tsize, ysize)
        matr[:, xsize, :] = zeros(tsize, ysize)
        matr[:, :, 1] = zeros(tsize, xsize)
        matr[:, :, ysize] = zeros(tsize, xsize)
    end
end

U[1, :, :] = zz

edge_set(U)

bc_edge = 0

function fin_diff(matr, n, j, l)
    if (j == 1 || l == 1 || j == size(matr, 2) || l == size(matr, 3))
        #println("(n),(n), (j), $(l)")
        #display(return_slice(U, n))
        return bc_edge
    elseif n == 1

        #println("(n),(n), (j), $(l)")
        return (r/2) * (matr[n, j-1, l] + matr[n, j+1, l] - 4 *matr[n, j, l] + matr[n, j, l - 1] + matr[n, j, l + 1]) + matr[n, j, l]

    elseif n != 1


        #println("(n),(n), (j), $(l)")
        return (r) * (matr[n, j-1, l] + matr[n, j+1, l] - 4 *matr[n, j, l] + matr[n, j, l - 1] + matr[n, j, l + 1]) + 2 * matr[n, j, l] - matr[n - 1, j, l]
    
    else

        return bc_edge
        
    end

end

function FDM_Solver(matr)
    println("Solving wave equation...")
    for n in range(1, size(matr, 1) - 1) # time
        for j in range(1, size(matr, 2)) # space (x)
            for l in range(1, size(matr, 3)) # space (y)
                matr[n + 1, j, l] = fin_diff(matr, n, j, l)
            end
        end
    end
    println("Done! Plotting results soon...")
end


FDM_Solver(U)

height_limit = ceil(maximum(zz))

p1 = surface(size=(900, 720), clims=(-height_limit, height_limit))

for i in Base.OneTo(tsize)
    global p1 = surface(xx, yy, reshape(U[i, :, :], length(U[i, :, :])), size=(900, 720), zlim=(-height_limit, height_limit), clims=(-height_limit,height_limit))
    display(p1)
    sleep(0.01)
end

# anim = @animate for i in Base.OneTo(tsize)
#     surface(xx, yy, reshape(U[i, :, :], length(U[i, :, :])), size=(900, 720), zlim=(-5, 5))
# end

# gif(anim, "C:\\Users\\sgtar\\Desktop\\Chase Programming\\Output\\gaussian_wave.gif", fps=30)
