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
