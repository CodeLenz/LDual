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