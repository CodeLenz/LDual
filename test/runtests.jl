using Test
using LDual

@testset "Sum" begin

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x+1.0)+(2.0+x)+(x+x)

    # Gets the derivative
    df = f(z).dual

    # test
    @test isapprox(df,4.0)


end

@testset "Subtract" begin


    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x-1.0)+(2.0-x)+(x-x)

    # Gets the derivative
    df = f(z).dual

    @test isapprox(df,0.0)

end

@testset "Multiplication" begin


  z = LDual.Dual(2.0, 1.0)

  # Tests the derivative 
  f(x) = (x*2.0)+(3.0*x)+(x*x)

  # Gets the derivative
  df = f(z).dual

  # Tests it
  @test isapprox(df,(5+(2*z.real)))

end


@testset "Division" begin
    
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative
    f(x) = (x/2.0)+(3.0/x)+(x/x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox(df,(0.5-(3/(z.real^2))))

end


@testset "Sine" begin
    
    # Initializes the input
    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 
    f(x) = 3*sin(x)

    # Gets the derivative
    df = f(z).dual

    # Tests it
    @test isapprox(df,(3*cos(z.real)))

end