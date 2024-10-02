using GeneralizedGenerated

# Function to create a typed function from an experssion and tuple of arguments
function make_function(typed_vars::Tuple, expr::Expr,name::Symbol)
    args = [Expr(:(::), var, typ) for (var, typ) in typed_vars] #Converts tuple of (::Symbol,::Datatype) into expr that represents x::Real
    
    #= function_expr = Expr(:function, Expr(:call, name, args...), expr) #Builds expr for the function f(args...) = expr - you can name the function with name
    
    eval(function_expr) # Evaluates the function generator to combine the lhs and rhs of the expr
    
    eval(function_expr.args[1].args[1]) # converts :(f(args...)=expr) into f(args...) and sends the function to global scope to be called from anywhere =#
    name=mk_function(args,[],expr)
end

#Builds a multi-threading loop around a 0D function
function create_multithreaded_func(func, result::Symbol, vars::Vector) 
    return quote #Effectively builds a more complex version of an expr
        Threads.@threads for i in 1:length($result) #Multithreading with the output vector interpoalted in
            $result[i] = $func($(vars...),i) #Interpolates the funtion call and variables, they have to be different as now need an x[i] next to all vector components and no typing
        end
    end
end

function main()
    grid = [0.0,1.0,2.0] # Builds a grid aka the z-grid of the slab
    las = :(1/2 * exp(-$grid[i])) #Defines an equation to represent a laser, interpolates the grid co-ordinates and leaves the [i] to be repalced in the multi-threading

    typed_vars = ((:x, Real), (:y, Real),(:i,Int)) #Create a nested tuple with the symbols and types for the various variables

    expr = :((x+y)/$las) #Expr for the scalar function - is assuming scalar inputs other than laser
    make_function(typed_vars,expr,:scal) # Builds the function called scal(x,y) -> wouldn't work unless i is defined in global scope in this case 

    vars = [:(x[i]),:(y[i])] # Builds the variables for the multi-threaded loop
    loop_expr = create_multithreaded_func(:scal, :du,vars) # Defines the expression for the function becoming multi-threaded

    parallel_vars = ((:x,Vector{<:Real}),(:y,Vector{<:Real}),(:du,Vector{<:Real})) # Defines the variables for the parallel version of the function, the index-wise variables are now their respective vectors
    make_function(parallel_vars, loop_expr,:threaded)#Builds the function called threaded(x,y,du) which replaces each value of du with the answer to scal(x,y), i is now available to be inserted for las


    x_arr = collect(range(1,10,length=3))  #A grid of x vals
    y_arr = collect(range(0,1,length=3)) #A grid of y vals
    du = zeros(length(x_arr)) # The output vector

    # Call the generated function
    threaded(x_arr,y_arr,du) # The calculation

    println(du) #Prints du as this is an in-place operation
end 
main()