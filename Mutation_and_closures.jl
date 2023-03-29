function main(A,D,d)
    a = 3
    
    B = A
    @show B === A
    
    B[2] = 20
    @show B,A

    B = A[1:2]
    @show B === A
    
    B2 = copy(B)
    @show B2,B
    @show B2 === B

    B[2] = 15.0
    @show B2,B,A
    @show B == A
    @show B === A

    B3 = copy(B)
    @show B3 === B2
    
    @show "setup",a,A,B,B2,B3,typeof(B)
    @show eltype(B)

    sharing!(B3)
    @show "shared",B3
    
    test = test_outer(A,B,B2,B3,a)
   
    # Similar to test = test_outer(A,B,B2,B3,a)
    # Uncomment to run
    #test = t -> test_inner(t,A,B,B2,B3,a)

    # Equivalent to test = t -> test_inner(t,A,B,B2,B3,a)
    # Uncomment to run
    #test(t) = test_inner(t,A,B,B2,B3,a)
    
    C = test(12.0)
    @show "after first test",A,a,B,B2,B3,C,typeof(C)
    
    #A = Array(1:3)
    A2 = A
    A2[2] = 25
    a = 30
    B = A
    B[1] = -1

    # Also try these
    # B2 = A.*20
    # B = [0,1,2]
    # B[1] = -10

    @show "before second test",A,a,A2,B,B2,B3,C,typeof(C)
    C = test(13)  
    @show "after second test",A,a,A2,B,B2,B3,C,typeof(C)
    

    push!(A2,42)
    @show "push",A,A2

    D[2] = 13
    d = 13
    E[2] = 13
    e = 13

    return d
end

function test_outer(A,B,B2,B3,a)
    @show "closure A,a,B,B2,B3",A,a,B,B2,B3
    function test_inner(t)
        B[2] = A[2] 
        B2 .= A*10
        B3 = A*a
        #a = 5
        @show "closed A,a,B,B2,B3",A,a,B,B2,B3
        return A*t
    end
end

function test_outer_lim(A::Array{T,1},B::Array{T,1},B2::Array{T,1},B3::Array{T,1},a::T) where T <: Real
    @show "closure A,a,B,B2,B3",A,a,B,B2,B3
    function test_inner(t)
        B[2] = A[2] 
        B2 .= A*10
        B3 = A*a
        #a = 5
        @show "closed A,a,B,B2,B3",A,a,B,B2,B3
        return A*t
    end
end

function test_inner(t,A,B,B2,B3,a)
    @show "B=",B
    B[2] = A[2] 
    B2 .= A*11
    B3 = A.*a
    @show "closed A,a,B,B2,B3",A,a,B,B2,B3
    return A.*t
end


function sharing!(B)
    @show "incoming value of B",B
    B .+= 1
    @show "pointer-like sharing",B
    # Also try this
   # sharing2!(B)
    return nothing
end

function sharing2!(C)
    @show "incoming value of B",C
    C .+= 1
    @show "pointer-like sharing",C
    return nothing
end


D = Array(1:2)
d = 2

E = collect(1:2)
e = 2

a = 2

f::Float64 = main(collect(1:2),D,d)
@show "globals",a,D,d,E,e,f
@show A


