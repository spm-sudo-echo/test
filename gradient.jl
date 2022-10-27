#this code will check and calculate the gradients 
# the gradient will be used in the calculations for the distances in the main code
# gardient makes code heavy

using Flux
using ForwardDiff
#single parameter function
#=
f(x) = 4x^2 + 3x + 2;
df(x) = gradient(f, x)[1]; # df/dx = 8x + 3
df(1)   
=#
#double parameter function
#=
c = [1.1 2.2;3.1 4.3]

#f(x, y) = x.* x + 4*y
 x = [2, 1];
 y = [2, 0];
 #f(x, y) = sum((x .- y).^2);
 f(x, y) = x.^2 .+ 4*y
  gs = ForwardDiff.gardient(() -> f(x, y),Flux.params(x, y)) 
=#



#=
f(x,y) = 4x^2 + 3x + 4y;
a= gradient(f,[1,1])
a[1]
#df(y) = gradient(f,y)[1];
=#
#f(1,1)
#x = zeros(length(2),2)

#--------------------------------------------------------------------------------------------
#Using direct vector
#=
x = [2.0;2.0]
f(x) = 4x[1]^2 + 3x[1] + 2x[2]

#df = ForwardDiff.gradient(f, [x[1],x[2]])
df = zeros(length(2),1)
#println("$df")
df=ForwardDiff.gradient(f, [x[1],x[2]])

g = sum(df)
=#
#--------------------------------------------------------------------------------------------
# Passing with the function
#=


function define(x::Array{Float64})
    f(x) = 4x[1]^2 + 3x[1] + 2x[2]
   
    return f(x)
end
define(x)

#--------------------------------------------------------------------------------------------
#x = [20.0,2.0]
# I am calling the functions and changing the x values randomly like in our code
# here one element is called 
#this code is working

function  scale(c::Float64)
    x = zeros(length(2),2)

    #x = c * rand(2,2) # it makes a 2x1 matric
    x = c * ones(2,2) # it makes a 2x1 matric
    check = grad(x) # check is an array of two elements check[1] and check[2] containing the gradient w.r.t x and y resp
    
     return check
end

function grad(x::Array{Float64})
     
    f(x) = 4x[1]^2 + 3x[1] + 2x[2]
    
    df = zeros(length(2),1)
    #println("$df")
    df=ForwardDiff.gradient(f, [x[1],x[2]])
    
    g = sum(df)    # gives the sum of all derivates and hence a gradient value 
    return df
end

#check=grad(x)  
result1=scale(10.0) # arguiment here is a scaling factor

=#
#--------------------------------------------------------------------------------------------
#x = [20.0,2.0]
# I am calling the functions and changing the x values randomly like in our code
# here mutiple element are called using the dot function

function  scale(c::Float64)
    #x = zeros(length(2),2)
    
    #x = c * rand(2,2) # it makes a 2x1 matric
    #x = c * ones(2,2) # it makes a 2x1 matric
    x= [rand(2) for _ in 1:5]
    check = grad.(x) # check is an array of two elements check[1] and check[2] containing the gradient w.r.t x and y resp
    
     return check
end

function grad(x::Array{Float64})
     
    f(x) = 4x[1]^2 + 3x[1] + 10*(sin.(x[2]))
    
    df = zeros(length(2),1)
    #println("$df")
    df=ForwardDiff.gradient(f, [x[1],x[2]])
    
    g = sum(df)    # gives the sum of all derivates and hence a gradient value 
    return df
end

#check=grad(x)  
result=scale(10.0) # arguiment here is a scaling factor

#--------------------------------------------------------------------------------------------
# this is above code but called multiple times
#=
function  scale(c::Float64, N::Int64)
    x = zeros(N,2)

   #x = c * rand(N,2) # it makes a 2x1Nmatric
 
    #x = c * ones(N,2) # it makes a 2xN matric
         
    x = [1.1 2.2;3.1 4.3]

    check = grad!(x) # check is an array of two elements check[1] and check[2] containing the gradient w.r.t x and y resp
    #f = 4x[1]^2 + 3x[1] + 2x[2]

    
     return check
end
#result1=scale(10.0,2)
=#

#x = [1.1 ;4.3]

#=
function grad!(x::AbsractArray{Float64,2},y::AbstractArray{Float64,2})
     
    f(x) = 4x[1]^2 + 3x[1] + 2x[2]
    #f = zeros(2,1)
    #f(x) = 4*x[:,1].*x[:,1]+ 3*x[:,1] + 2*x[:,2]
    #f(x) = 4*x[:,1].*x[:,1]+ 3*x[:,1] 
    #g(x)=2*x[:,2]
    #f = 4*x;
  
  
    k1= y[1,1]
    k2= y[1,2]

    k3= y[2,1]
    k4= y[2,2]
    #println("$df")
    #df=ForwardDiff.gradient(f, [1.0, 1.0])
    #df=ForwardDiff.derivative(f, 1.0)
   # df=ForwardDiff.gradient(f, [x[1],x[2]])
   df=ForwardDiff.gradient(f, [y[1,1],y[1,2], y[2,1],y[2,2]])
    #g = sum(df)    # gives the sum of all derivates and hence a gradient value 
    return k1, k2, k3, k4
end
=#

#=
x = z[:,1]
y=  z[:,2]
f= zeros(2,1)


x = zeros(2,2)
#f = zeros(2,1)
x = [1.1 2.2;3.1 4.3]

function cal(x::Array{Float64,2})
    for i in 1:2:3 
        p = 1*x[i]*x[i] + 1*x[i+1]
    end
 return p
end

df = ForwardDiff.gradient(cal, [1.0,2.0])
  


#=
function val!(x::Array{Float64,2})
    for i in 1:2
    
        f = 4*x[i] + 4*x[i+1]
        df = ForwardDiff.gradient(f, [x[i],x[i+1]])
        #println("f[$i]")
     end
     return df  
end

val!(x)
=#
#f(x) = 4x.*x + 3x+ 2y
#f(x) = 4x[:,1]^2 + 3x[:,1] + 2x[:,2]
#ForwardDiff.gradient(f,x::AbstractArray)
#check=grad!(x,x)  
#result1=scale(10.0,2) # arguiment here is a scaling factor
=#
#-------------------------------------------------------------


