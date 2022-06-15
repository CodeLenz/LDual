#
# Implementação de diferenciação automática - Forward Differences
# Eduardo Lenz 11/05/2017
#
# Para usar:: using FD
#
# x = Dual(1.0, 1.0)
# cos(x)^2 --> vai dar um numero dual, com a segunda parte sendo a derivada, ie 2*cos(x.real)*(-sin(x.real))
#
# Ultima revisão: 18/05/2017 (FUNÇÕES HIPERBÓLICAS)
# Ultima revisão: 16/09/2019 (Julia 1.*)
# Ultima revisão: 15/06/2022 (Arrumando para colocar no git) 
# 
# using LDual
#
# f(x) = cos(2*x) + 3*x
# x = Dual(1.0,1.0)
#
# Se rodarmos f(x) obteremos  Dual(f(x), df(x))
#
#
module LDual

    # Exporta os metodos deste módulo
    export Dual, +, *, -,/,sin,cos, exp, log, sqrt, abs

    # exporta as definições basicas de um e zero para o tipo
    # isto estende automaticamente os comandos zeros e ones
    export one, zero

    # Metodo de conversão Int e Float para dual
    export convert

    # Exporta as operações sobre arrays
    export transpose, dot, Rand, norm

    # Define o tipo básico - Número Dual
    # real -> parte real
    # dual -> parte dual
    struct  Dual #<:Number
	     real::Float64
	     dual::Float64
    end

    #################################################
    #       DEFINIÇÕES DAS IDENTIDADES
    #################################################

    # Identidade multiplicativa
    import Base:one
    function one(T::Type{Dual})
        return Dual(1.0,0.0)
     end

    # Identidade aditiva
    import Base:zero
    function zero(T::Type{Dual})
        return Dual(0.0,0.0)
     end

    import Base:zero
    function zero(a::Dual)
        return Dual(0.0,0.0)
     end

    # Converte um número "normal" para Dual
    import Base:convert
    function convert(::Type{Dual},x::Number)
       Dual(x,0.0)
    end


    ############################################################
    #       DEFINIÇÕES DAS OPERAÇÕES BÁSICAS
    #
    # Lembrando que a parte dual é sempre a derivada
    # da operação em relação a entrada e a parte
    # primal (real) é a própria operação * dual (liga/desliga)
    ############################################################

    # SOMA
    import Base:+
    function +(x::Dual, y::Dual)
	      p = x.real + y.real
	      d = x.dual + y.dual
         Dual( p , d )
    end

    # PRODUTO
    import Base:*
    function *(x::Dual, y::Dual)
	    p = x.real*y.real
	    d = x.real*y.dual + x.dual*y.real
	    Dual( p , d )
    end


    # NEGATE
    import Base:-
    function -(x::Dual)
	    p = -x.real
	    d = -x.dual
	    Dual( p , d )
    end

    # SUBTRAI
    import Base:-
    function -(x::Dual,y::Dual)
	    y1 = -y
	    return x + y1
    end


    #X/Y
    import Base:/
    function /(x::Dual,y::Dual)

	   # Evitamos y.real nulo
	   @assert y.real!=0

	   # Calcula a parte 1/y
	   um   =  1.0/y.real
	   dois = -1*y.dual/(y.real^2)

	   # Devolve o valor
	   return x*Dual(um,dois)

    end

    # SENO
    import Base:sin
    function sin(x::Dual)
	    p = sin(x.real)
	    d = cos(x.real)*x.dual
	    Dual( p , d )
    end

    # COSSENO
    import Base:cos
    function cos(x::Dual)
	    p = cos(x.real)
	    d = -sin(x.real)*x.dual
	    Dual( p , d )
    end


    # EXP
    import Base:exp
    function exp(x::Dual)
	    p = exp(x.real)
	    d = p*x.dual
	    Dual( p , d )
    end

    # TANH - em função de exp
    import Base:tanh
    function tanh(x::Dual)

           e2x = exp(2*x)
           return (e2x-1.0)/(e2x+1.0)

    end

    # SINH - em função de exp
    import Base:sinh
    function sinh(x::Dual)

           return (1/2)*(exp(x)-exp(-x))

    end

    # COSH - em função de exp
    import Base:cosh
    function cosh(x::Dual)

           return (1/2)*(exp(x)+exp(-x))

    end


    # LOG (LN..)
    import Base:log
    function log(x::Dual)

	    # Evita divisao por zero
	    @assert x.real != 0

	    p = log(x.real)
	    d = x.dual/x.real
	    Dual( p , d )
    end

    # SQRT
    import Base:sqrt
    function sqrt(x::Dual)

	    # Evita divisao por zero
	    @assert x.real != 0

	    p = sqrt(x.real)
	    d = x.dual/(2*sqrt(x.real))
	    Dual( p , d )
    end

    # ABS
    import Base:abs
    function abs(x::Dual)

        # Evita divisao por zero
        @assert x.real != 0

        p = abs(x.real)
        d = x.dual*(x.real/abs(x.real))

        Dual( p , d)
    end

    #################################################
    #               CORNER CASES
    #################################################

    # FACILITA A MULTIPLICACAO POR UMA CTE NAO DUAL,
    function *(x, y::Dual)

	    # Trasforma x em dual
	    x1 = Dual(x,0.0)

	    # Multiplica
	    x1*y
    end

    function *(x::Dual, y)

	    # Trasforma y em dual
	    y1 = Dual(y,0.0)

	    # Multiplica
	    x*y1
    end

    # FACILITA A SOMA POR UMA CTE NAO DUAL
    function +(x,y::Dual)

	   # Transforma x em dual
	   x1 = Dual(x,0.0)

	   # Soma
	   x1 + y

    end

    function +(x::Dual,y)

	   # Transforma y em dual
	   y1 = Dual(y,0.0)

	   # Soma
	   x + y1

    end

    # FACILITA A SUBTR POR UMA CTE NAO DUAL
    function -(x,y::Dual)

	   # Transforma x em dual
	   x1 = Dual(x,0.0)

	   # Subtrai
	   x1 - y

    end

    function -(x::Dual,y)

	   # Transforma y em dual
	   y1 = Dual(y,0.0)

	   # Subtrai
	   x - y1

    end


    #################################################
    #       DEFINIÇÕES VETORIAS - Em teste
    #################################################

    # TRANSPOSE
    import LinearAlgebra:transpose
    function transpose(A::Array{Dual})

       # Verifica a dimensão de A
       dims = size(A)

       # Vamos definir transposta para vetor e matriz
       if length(dims)==1
          return reshape(A, 1, length(A))
       elseif length(dims)==2
          return  permutedims(A, (2, 1))
       else
           error("Transpose dual -> somente arrays com dimensão 1 ou 2")
       end


    end


    # PRODUTO INTERNO
    import LinearAlgebra:dot
    function dot(A::Vector{Dual},B::Vector{Dual})
	    soma = Dual(0.0, 0.0)
	    for i=1:length(A)
		    soma += A[i]*B[i]
	    end
        return soma
     end


    # Produto de um escalar por um array
    function *(x::Float64,A::Array{Dual})


           # Aloca a saida
           V = zeros(A)

           # Transforma o x em um numero dual
           x1 = Dual(x,0.0)

           # Aplica o produto em cada uma das posições
           for i in eachindex(A)
             V[i] = x1*A[i]
           end

           return V

    end


    # Produto de um dual por um array
    function *(x::Dual,A::Array{Dual})


           # Aloca a saida
           V = zeros(A)

           # Aplica o produto em cada uma das posições
           for i in eachindex(A)
             V[i] = x*A[i]
           end

           return V

    end


    # RAND PARA ARRAYS - Aqui estou com uma dificuldade pois está dando
    # conflito com
    #   rand(T::Type, d1::Integer, dims::Integer...) at random.jl:232
    #   rand(T::Type{FD.dual}, dims...) at /home/lenz/Dropbox/dif_automatica.jl:245
    # VOU USAR Rand
    #import Base.rand
    function Rand(T::Type{Dual}, dims...)

         # Aloca array
         V = Array{T}(undef,dims)

         # Inicializa todas as posições como um e dual nulo
         for i in eachindex(V)
            V[i] = Dual(rand(), 0.0)
         end

         # Retorna o V completo
         return V

    end


    #
    # calcula a norma 2 - Só a 2 !
    #
    import LinearAlgebra:norm
    function norm(A::Vector{Dual},p=2)

       # Verifica se p é 2 (norma Euclidiana)
       p==2 || throw("LDual::norm nonly p=2 is implemented")

       # Converte para vetor
       a = vec(A)

       # Faz o produto interno
       prod = dot(a,a)

       # Retorna a raiz
       sqrt(prod)


    end




end #LDual
