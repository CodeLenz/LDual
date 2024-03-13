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

    df_analytic = 3*cosh(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the hyperbolic cosine is wrong. The",
         " value of the analytical derivative is ", df_analytic, " and",
         " the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the hyperbolic cosine is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the natural logarithm 

function test_naturalLogarithm()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*log(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3/z.real

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the natural logarithm is wrong. The",
         " value of the analytical derivative is ", df_analytic, " and",
         " the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the natural logarithm is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the logarithm of base 10 

function test_logarithm10()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*log10(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3/(log(10)*z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the logarithm of base 10 is wrong. ",
         "The value of the analytical derivative is ", df_analytic, " ",
         "and the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the logarithm of base 10 is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the logarithm of any base 

function test_logarithmAny()

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

    if abs(df-df_analytic)>1E-5

        println("The derivative of the logarithm of base ", base_log, 
         " is wrong. The value of the analytical derivative is ",
         df_analytic, " and the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the logarithm of base ", base_log,
         " is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the square root

function test_squareRoot()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*sqrt(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*0.5/sqrt(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the square root is wrong. The value",
         " of the analytical derivative is ", df_analytic, " and the d",
         "ual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the square root is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the absolute value

function test_absoluteValue()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3*abs(x)

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*sign(z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of the absolute value is wrong. The va",
         "lue of the analytical derivative is ", df_analytic, " and th",
         "e dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the absolute value is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of x^y w.r.t. x

function test_xRaisedY()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = x^3.0

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = 3*(z.real^2)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of x^y is wrong. The value of the anal",
         "ytical derivative is ", df_analytic, " and the dual derivati",
         "ve is ", df, ".\n")

        return false

    else 

        println("The derivative of x^y is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of y^x w.r.t. x

function test_yRaisedX()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = 3^x

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = log(3)*(3^z.real)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of y^x is wrong. The value of the anal",
         "ytical derivative is ", df_analytic, " and the dual derivati",
         "ve is ", df, ".\n")

        return false

    else 

        println("The derivative of y^x is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of x^x w.r.t. x

function test_xRaisedX()

    # Initializes the input

    z = LDual.Dual(2.0, 1.0)

    # Tests the derivative 

    f(x) = x^x

    # Gets the derivative

    df = f(z).dual

    # Evaluates the derivative analytically

    df_analytic = (z.real^z.real)*(log(z.real)+1)

    # Tests it

    if abs(df-df_analytic)>1E-5

        println("The derivative of x^x is wrong. The value of the anal",
         "ytical derivative is ", df_analytic, " and the dual derivati",
         "ve is ", df, ".\n")

        return false

    else 

        println("The derivative of x^x is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of x*A

function test_xMultArray()

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

    if norm(df-df_analytic)>1E-5

        println("The derivative of dual multiplied by array is wrong. ",
         "The value of the analytical derivative is ", df_analytic, " ",
         "and the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of dual multiplied by array is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of A*X, where X is dual

function test_XArrayMultArray()

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

    if norm(df-df_analytic)>1E-5

        println("The derivative of dual array multiplied by array is w",
         "rong. The value of the analytical derivative is ",
         df_analytic, " and the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of dual array multiplied by array is r",
         "ight.\n")

        return true 

    end

end

# Defines a function to test the derivative of the dot product

function test_dotX()

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

    if norm(df-df_analytic)>1E-5

        println("The derivative of the dual dot product is wrong. The ",
         "value of the analytical derivative is ", df_analytic, " and ",
         "the dual derivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the dual dot product is right.\n")

        return true 

    end

end

# Defines a function to test the derivative of the norm p=2

function test_normP2()

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

    if norm(df-df_analytic)>1E-5

        println("The derivative of the norm is wrong. The value of the",
         " analytical derivative is ", df_analytic, " and the dual der",
         "ivative is ", df, ".\n")

        return false

    else 

        println("The derivative of the norm is right.\n")

        return true 

    end

end