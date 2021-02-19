using Images, ImageMagick, LinearAlgebra
include("matrixsQ.jl")

"""
    DTC1(x)
Applies DTC type-II to vectors.
# Argument
- `x::Array{Float64,1}`
- `x::Ntuple`
"""
function DTC1(x)
    N=length(x)
    X=zeros(N)
    for k=0:N-1
        for n=0:N-1
            X[k+1] = X[k+1]+x[n+1]*cos(π/N*(n+1/2)*k)
        end
    end
    return X
end

"""
    invDTC1(x)
Inverse DTC type-II to vectors.
# Argument
- `x::Array{Float64,1}`
- `x::Ntuple`
"""
function invDTC1(X)
    N=length(X)
    x=(0.5*X[1]).*ones(N)
    for k=0:N-1
        for n=1:N-1
            x[k+1] = x[k+1]+X[n+1]*cos((π/N)*(k+1/2)*n)
        end
    end
    return (2.0/N).*x
end

"""
    DTC2(A)
Applies DTC type-II to arrays.
# Argument
- `x::Array{Float64,2}`
"""
function DTC2(x)
    a,b=size(x)
    M=zeros(a,b)
    for j=1:b
        M[:,j]=DTC1(x[:,j])
    end
    for i=1:a
        M[i,:]=DTC1(M[i,:])
    end
    return M
end

"""
    invDTC2(A)
Inverse DTC type-II to arrays.
# Argument
- `x::Array{Float64,2}`
"""
function invDTC2(X)
    a,b=size(X)
    N=zeros(a,b)
    for i=1:a
        N[i,:]=invDTC1(X[i,:])
    end
    for j=1:b
        N[:,j]=invDTC1(N[:,j])
    end
    return N
end

"""
    divsixteen(M)
Cut the edge elements of the M to transform the size multiple of 16.
# Argument
- `M::Array{,2}`
"""
function divsixteen(M)
    N = Vector{Array{Float64,2}}(undef,3)
    l,c=size(M[1,:,:])
    m=Int(floor(l/16))
    n=Int(floor(c/16))
    m=16*m
    n=16*n
    x=l-m
    y=c-n
    if x==0
        up=1
        down=0
    else
        up=Int(floor(x/2))+1
        down=Int(floor(x/2))
    end
    if y==0
        left=1
        right=0
    else
        left=Int(floor(y/2))+1
        right=Int(floor(y/2))
    end
    for k=1:3
        N[k]=M[k,up:m+down,left:n+right]
    end
    O=zeros(3,size(N[1])[1],size(N[1])[2])
    for p=1:3
        O[p,:,:]=N[p]
    end
    return O,x,y
end

"""
    ad_divsixteen(M)
Execute divsixteen function and print advices.
# Argument
- `M::Array{,2}`
"""
function ad_divsixteen(M)
    O,x,y=divsixteen(M)
    l,c=size(M[1,:,:])
    if x!=0
        if y!=0
            println("Attention, have been removed $x rows e $y columns.")
            println("Will not be counted $(3*(x*c+y*l-x*y)) data.")
        else
            println("Attention, have been removed $x rows.")
            println("Will not be counted $(3*x*c) data.")
        end
    else
        if y!=0
            println("Attention, have been removed $y columns.")
            println("Will not be counted $(3*y*l) data.")
        end
    end
    return O
end

"""
    RGB_to_YUV(A)
Transform RGB (red, green, blue) into YUV (luminance and chrominances).
# Argument
- `A::Array{Float64,3}`: size(A)=(3,any,any)
"""
function RGB_to_YUV(RGB) #transforma RGB em YUV (luminosidade e cromancias)
    Y=0.299*RGB[1,:,:]+0.587*RGB[2,:,:]+0.114*RGB[3,:,:]
    U=-0.1687*RGB[1,:,:]-0.3313*RGB[2,:,:]+0.5*RGB[3,:,:]
    V=0.5*RGB[1,:,:]-0.4187*RGB[2,:,:]-0.0813*RGB[3,:,:]
    return Y, U, V
end

"""
    YUV_to_RGB(Y,U,V)
Transform YUV (luminance and chrominances) into RGB (red, green, blue).
# Arguments
- `Y::Array{,3}`: Luminance
- `U::Array{,3}`: Chrominance 1
- `V::Array{,3}`: Chrominance 2
"""
function YUV_to_RGB(Y,U,V) #transforma YUV em RGB
    R=Y + 1.402 * V
    G=Y - 0.34414 * U -0.71414 * V
    B=Y + 1.772 * U
    return R, G, B
end

"""
    cut(A)
Remove all rows and columns with odd index.

Reduces array size to 1/4.
# Argument
- `A::Array{,2}`
"""
function cut(A) #reduz o tamanho das matrizes de cromancias
    m,n=size(A)
    m=Int(floor(m/2))
    B=zeros(m,Int(n))
    for i=1:m
        B[i,:]=A[2*i,:]
    end
    n=Int(floor(n/2))
    C=zeros(m,n)
    for j=1:n
        C[:,j]=B[:,2*j]
    end
    return C
end

"""
    expand(A)
Duplicates all rows and columns.

Increases the size of the matrix four times.
# Argument
- `A::Array{,2}`
"""
function expand(A) #retorna as matrizes de cromancia para o tamanho comum
    m,n=size(A)
    B=zeros(2*m,n)
    for i=1:m
        B[2*i,:]=A[i,:]
        B[2*i-1,:]=A[i,:]
    end
    C=zeros(2*m,2*n)
    for i=1:n
        C[:,2*i]=B[:,i]
        C[:,2*i-1]=B[:,i]
    end
    return C
end

"""
    block(A)
Makes a array subdivide into 8x8 blocks, stored in a vector.

Ignores final parts not divisible by 8.
# Argument
- `A::Array{,2}`
"""
function block(A) #molda os blocos 8x8
    (m,n)=size(A)
    mblocos = floor(Int,m/8)
    nblocos = floor(Int,n/8)
    B = Vector{Array{Float64,2}}(undef,mblocos*nblocos)
    k=1
    for i=1:8:m-7
        for j=1:8:n-7
            B[k] = A[i:i+7,j:j+7]
            k=k+1
        end
    end
    return B, mblocos, nblocos
end

"""
    frame(A,x,y)
Assembles a square array from a vector where each element is a 8x8 array.

Inverse of `block`.
# Argument
- `A::Array{Array{,2},1}`
- `x::Int64` number of blocks in rows
- `y::Int64` number of blocks in columns
"""
function frame(Gb,r,s) #une os bloquinhos 8x8
    (x,)=size(Gb)
    if x!=r*s
        r=2*r
        s=2*s
    end
    B=zeros(8*r,8*s)
    aux=1
    for h=0:r-1
        for k=0:s-1
            for i=1:8
                for j=1:8
                    B[i+8*h,j+8*k]=Gb[aux][i,j]
                end
            end
            aux=aux+1
        end
    end
    return B
end

"""
    quantization(A,Q)
Performs element-by-element division of A by Q and round results.
# Argument
- `A::Array{Float64,2}` size(A)=(8,8)
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function quantization(M,Q)
    for i=1:8
        for j=1:8
            M[i,j]=round(M[i,j]/Q[i,j])
        end
    end
    return M
end

"""
    unquantized(A,Q)
Performs element-by-element multiplication of A by Q.
# Argument
- `A::Array{Float64,2}` size(A)=(8,8)
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function unquantized(N,Q)
    for i=1:8
        for j=1:8
            N[i,j]=N[i,j]*Q[i,j]
        end
    end
    return N
end

"""
    compress(image,Q)
Performs the previous functions to compress an image.

Returns compressed arrays.
# Argument
- `image::String`: archive name .jpeg or .png
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function compress(imagem,Q) #resulta nos dados que serão armazenados
    img=channelview(load(imagem))
    Y,U,V=RGB_to_YUV(ad_divsixteen(img))
    Yb,r,s=block(Y)
    Ub,r,s=block(cut(U))
    Vb,r,s=block(cut(V))
    Yb=255*Yb
    Ub=255*Ub
    Vb=255*Vb
    m=length(Ub)
    n=length(Yb)

    for i=1:m
        Ub[i]=quantization(DTC2(Ub[i]),Q)
        Vb[i]=quantization(DTC2(Vb[i]),Q)
    end
    for i=1:n
        Yb[i]=quantization(DTC2(Yb[i]),Q)
    end

    return Yb, Ub, Vb, r, s
end

"""
    decompress(Y,U,V,x,y,Q)
Decompress the arrays.

Return full arrays.
# Arguments
- `Y::Array`: Compressed luminance
- `U::Array`: Compressed chrominance 1
- `V::Array`: Compressed chrominance 2
- `x::Int64` number of blocks in rows to use in function frame step
- `y::Int64` number of blocks in columns to use in function frame step
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function decompress(Yb,Ub,Vb,r,s,Q)
    m=length(Ub)
    n=length(Yb)

    for i=1:m
        Ub[i]=invDTC2(unquantized(Ub[i],Q))
        Vb[i]=invDTC2(unquantized(Vb[i],Q))
    end
    for i=1:n
        Yb[i]=invDTC2(unquantized(Yb[i],Q))
    end
    U=expand(frame(Ub,r,s))
    V=expand(frame(Vb,r,s))
    Y=frame(Yb,r,s)
    return Y, U, V
end

"""
    jpeg_compress(image,Q)
Performs compress and decompress functions to returns compressed image.
# Argument
- `image::String`: archive name .jpeg or .png
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function jpeg_compress(figura,Q)
    Y,Ub,Vb,r,s=compress(figura,Q)
    Y,U,V=decompress(Y,Ub,Vb,r,s,Q)
    R,G,B=YUV_to_RGB(Y,U,V)
    return colorview(RGB,R/255,G/255,B/255)
end
