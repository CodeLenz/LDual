# Tests for LDual

using LinearAlgebra

include("LDual.jl")

# Defines a function to test the derivative of sum

function test_sum()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative

    f(x) = (x+1.0)+(2.0+x)+(x+x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=4.0

        println("The derivative of the sum is wrong. The value of the ",
         "analytical derivative is ", 4, " and the dual derivative is ",
         df, ".\n")

        return false

    else 

        println("The derivative of the sum is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of subtraction

function test_subtraction()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative

    f(x) = (x-1.0)+(2.0-x)+(x-x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=0.0

        println("The derivative of the subtraction is wrong. The value",
         " of the analytical derivative is ", 0, " and the dual deriva",
         "tive is ", df, ".\n")

        return false

    else 

        println("The derivative of the subtraction is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of multiplication

function test_multiplication()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = (x*2.0)+(3.0*x)+(x*x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=(5+(2*z.real))

        println("The derivative of the multiplication is wrong. The va",
         "lue of the analytical derivative is ", (5+(2*z.real)), " and",
         " the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the multiplication is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of division

function test_division()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative

    f(x) = (x/2.0)+(3.0/x)+(x/x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=(0.5-(3/(z.real^2)))

        println("The derivative of the division is wrong. The value of",
         " the analytical derivative is ", (0.5-(3/(z.real^2))), " and",
         " the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the division is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of sine function

function test_sine()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*sin(x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=(3*cos(z.real))

        println("The derivative of the sine is wrong. The value of the",
         " analytical derivative is ", (3*cos(z.real)), " and the dual",
         " derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the sine is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of cosine function

function test_cosine()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*cos(x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if df!=(-3*sin(z.real))

        println("The derivative of the cosine is wrong. The value of t",
         "he analytical derivative is ", (-3*sin(z.real)), " and the d",
         "ual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the cosine is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the tangent function

function test_tangent()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*tan(x)

    # Gets the derivative

    df = f(z).dual

    # Tests it

    if abs(df-(3/(cos(z.real)^2)))>1E-5

        println("The derivative of the tangent is wrong. The value of ",
         "the analytical derivative is ", (3/(cos(z.real)^2)), " and t",
         "he dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the tangent is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the exponential function

function test_exponential()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*exp(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*exp(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the exponential is wrong. The value",
         " of the analytical derivative is ", df_analytic, " and the d",
         "ual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the exponential is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the hyperbolic tangent 

function test_tanh()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*tanh(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*(1-(tanh(z.real)^2))

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the hyperbolic tangent is wrong. Th",
         "e value of the analytical derivative is ", df_analytic, " an",
         "d the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the hyperbolic tangent is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the hyperbolic sine 

function test_sinh()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*sinh(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*cosh(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the hyperbolic sine is wrong. The v",
         "alue of the analytical derivative is ", df_analytic, " and t",
         "e dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the hyperbolic sine is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the hyperbolic cosine 

function test_cosh()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*sinh(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*sinh(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the hyperbolic cosine is wrong. The",
         " value of the analytical derivative is ", df_analytic, " and",
         " the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the hyperbolic sine is right.\n")

        return true 

    end

end