include("jpeg_compress.jl")
include("quality.jl")

"""
    nonull(A)
Counts how many non-null elements are in an array.
# Argument
- `A::Array{Float64,2}`
"""
function nonull(A)
    n=0
    a,=size(A)
    for i=1:a
        for j=1:8
            for k=1:8
                if A[i][j,k] != 0
                    n=n+1
                end
            end
        end
    end
    return n
end

"""
    jpeg_results(image,Q)
Performs as final function and returns numerical data related to compression.
# Argument
- `image::String`: archive name .jpeg or .png
- `Q::Array{Int64,2}` size(Q)=(8,8) is the quantization matrix.
"""
function jpeg_results(figura,Q) #pega os dados e mostra a imagem
    Y,Ub,Vb,r,s=compress(figura,Q)
    F=nonull(Y)+nonull(Ub)+nonull(Vb)
    F=2*F+65
    Y,U,V=decompress(Y,Ub,Vb,r,s,Q)
    R,G,B=YUV_to_RGB(Y/255,U/255,V/255)
    m,n=size(R)
    I=3*m*n
    println("Initially was used $I data.")
    println("When compressing are used $F data.")
    t=F/I*100
    t=round(t,digits=3)
    println("Storage equal $t% of the original.")
    imgf=colorview(RGB,R,G,B)
    imgO=divsixteen(channelview(load(figura)))[1]
    MSR=quality(imgO,channelview(imgf))
    println("Loss of quality rate is $MSR. (0 = without, 1 = maximum)")
    return imgf
end
