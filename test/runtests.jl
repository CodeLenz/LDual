using Test
using LDual

@testset "sum" begin

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x+1.0)+(2.0+x)+(x+x)

    # Gets the derivative
    df = f(z).dual

    # test
    @test isapprox(df,4.0)


end

@testset "subtraction" begin


    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x-1.0)+(2.0-x)+(x-x)

    # Gets the derivative
    df = f(z).dual

    @test isapprox(df,0.0)

end

@testset "multiplication" begin


  z = LDual.Dual(2.0, 1.0)

  # Tests the derivative 
  f(x) = (x*2.0)+(3.0*x)+(x*x)

  # Gets the derivative
  df = f(z).dual

  # Tests it
  @test isapprox(df,(5+(2*z.real)))

end


@testset "division" begin
    
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x/2.0)+(3.0/x)+(x/x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox(df,(0.5-(3/(z.real^2))))

end


@testset "sin" begin
    
    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*sin(x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox(df,(3*cos(z.real)))

end

@testset "cos" begin
    

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*cos(x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox(df,(-3*sin(z.real)))

end

@testset "tan" begin

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*tan(x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox( df - (3/(cos(z.real)^2)), 0.0; atol=1E-12)

end

@testset "exp" begin
   
    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*exp(x)

    # Gets the derivative
    df = f(z).dual

    # Evaluates the derivative analytically
    df_analytic = 3*exp(z.real)

    # Tests it
    @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "tanh" begin

    # Initializes the input
     z = LDual.Dual(2.0, 1.0)

     # Tests the derivative 
     f(x) = 3*tanh(x)
 
     # Gets the derivative
     df = f(z).dual
 
     # Evaluates the derivative analytically
     df_analytic = 3*(1-(tanh(z.real)^2))
 
     # Tests it
     @test isapprox(df-df_analytic,0,atol=1E-12)
 
end

@testset "sinh" begin

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*sinh(x)

    # Gets the derivative
    df = f(z).dual

    # Evaluates the derivative analytically
    df_analytic = 3*cosh(z.real)

    # Tests it
    @test isapprox(df-df_analytic,0,atol=1E-12)

end


@testset "cosh" begin

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*sinh(x)

    # Gets the derivative
    df = f(z).dual

    # Evaluates the derivative analytically
    df_analytic = 3*cosh(z.real)

    # Tests it
    @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "log" begin

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*log(x)

    # Gets the derivative
    df = f(z).dual

    # Evaluates the derivative analytically
    df_analytic = 3/z.real

    # Tests it
    @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "log10" begin
    
     # Initializes the input
     z = LDual.Dual(2.0, 1.0)

     # Tests the derivative 
     f(x) = 3*log10(x)
 
     # Gets the derivative
     df = f(z).dual
 
     # Evaluates the derivative analytically
     df_analytic = 3/(log(10)*z.real)
 
     # Tests it
     @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "log any base" begin

     # Initializes the input
     z = LDual.Dual(2.0, 1.0)

     # Defines the base of the logarithm
     base_log = 2.0
 
     # Tests the derivative 
     f(x) = 3*log(base_log, x)
 
     # Gets the derivative
     df = f(z).dual
 
     # Evaluates the derivative analytically
     df_analytic = 3/(log(base_log)*z.real)
 
     # Tests it
     @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "sqrt" begin

     # Initializes the input
     z = LDual.Dual(2.0, 1.0)

     # Tests the derivative 
     f(x) = 3*sqrt(x)
 
     # Gets the derivative
     df = f(z).dual
 
     # Evaluates the derivative analytically
     df_analytic = 3*0.5/sqrt(z.real)
 
     # Tests it
     @test isapprox(df-df_analytic,0,atol=1E-12)

end

@testset "abs" begin

    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*abs(x)

    # Gets the derivative
    df = f(z).dual

    # Evaluates the derivative analytically
    df_analytic = 3*sign(z.real)

    # Tests it
    @test isapprox(df-df_analytic,0,atol=1E-12)

end

# Defines a @testset to test the derivative of x^y w.r.t. x

@testset "x^y" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = x^3.0

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*(z.real^2)

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of y^x w.r.t. x

@testset "y^x" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3^x

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = log(3)*(3^z.real)

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of x^x w.r.t. x

@testset "x^x" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = x^x

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = (z.real^z.real)*(log(z.real)+1)

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of x*A

@testset "x*Array" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Initializes a dual matrix and a number matrix

    B_dual = [z z^2; z^3 z^4]

    B = [1 2; 3 4]

    # Tests the derivative 

    f(x, A) = (3*A)+(A*2)+(x*B)+(B*x)+(A*x)+(x*A)

    # Gets the derivative

    df = LDual.get_dualComponents(f(z, B_dual))

    # Evaluates the derivative analytically

    dA_analytic = [1.0 2*z.real; 3*(z.real^2) 4*(z.real^3)]

    df_analytic = ((5+(2*z.real))*dA_analytic)+(2*B)+(2*LDual.
     get_realComponents(B_dual))

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of A*X, where X is dual

@testset "X*A" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Initializes a dual matrix and a number matrix

    B_dual = [z z^2; z^3 z^4]

    B = [1 0; 0 1]

    # Tests the derivative 

    f(A) = (B*A)+(A*B)+(A*A)

    # Gets the derivative

    df = LDual.get_dualComponents(f(B_dual))

    # Evaluates the derivative analytically

    dA_analytic = [1.0 2*z.real; 3*(z.real^2) 4*(z.real^3)]

    df_analytic = (B*dA_analytic)+(dA_analytic*B)+(dA_analytic*LDual.get_realComponents(B_dual))+(
     LDual.get_realComponents(B_dual)*dA_analytic)

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of the dot product

@testset "dot" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Initializes a dual matrix and a number matrix

    b_dual = [z; z^2; z^3]

    b = [1; 2; 3]

    a = [2; 3; 4]

    # Tests the derivative 

    f(x) = dot(b,x)+dot(x,a)+dot(x,x)

    # Gets the derivative

    df = f(b_dual).dual

    # Evaluates the derivative analytically

    df_analytic = 3+(10*z.real)+(21*(z.real^2))+(2*z.real)+(4*(z.real^3))+(6*(z.real^5))

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end

# Defines a @testset to test the derivative of the norm p=2

@testset "norm p=2" begin

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Initializes a dual matrix and a number matrix

    b_dual = [z; z^2; z^3]

    # Tests the derivative 

    f(x) = norm(x)

    # Gets the derivative

    df = f(b_dual).dual

    # Evaluates the derivative analytically

    df_analytic = (((2*z.real)+(4*(z.real^3))+(6*(z.real^5)))/(2*sqrt((z.real 
     ^2)+(z.real^4)+(z.real^6))))

    # Tests it

    @test isapprox(df-df_analytic,0,atol=1E-5)

end