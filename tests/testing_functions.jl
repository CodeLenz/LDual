# Tests for LDual

using LinearAlgebra

# Defines a function to test the derivative of the sum of a function

function test_sum()

    # Initializes the input

    z = Dual(2.0, 1.0)

    # Tests the derivative of f(x) = x+1

    f(x) = (x+1.0)+(2.0+x)

    # Gets the derivative

    df = f(z)

    # Tests it

    if df!=2.0

        println("The derivative of the sum is wrong. The value of the ",
         "analytical derivative is ", 2, " and the dual derivative is ",
         df, ".\n")

        return false

    else 

        println("The derivative of the sum is right.\n")

        return true 

    end

end