using Images, ImageMagick, LinearAlgebra
include("quality.jl")

function svdE_results(i,p)

    #compressao
    img=channelview(load(i))
    imgY,PimgU,PimgV=RGB_to_YUV(img)
    imgU=cut(PimgU)
    imgV=cut(PimgV)
    UY,SY,VY=svd(imgY)
    UU,SU,VU=svd(imgU)
    UV,SV,VV=svd(imgV)
    k=Int(ceil(length(SY)*p/100))
    h=Int(ceil(length(SV)*p/100))
    Y=UY[:,1:k]*diagm(0=>SY[1:k])*transpose(VY[:,1:k])
    U=UU[:,1:h]*diagm(0=>SU[1:h])*transpose(VU[:,1:h])
    V=UV[:,1:h]*diagm(0=>SV[1:h])*transpose(VV[:,1:h])

    #contagem
    vav,a,b=size(img)
        println("Anteriormente precisavamos armazenar $(3*a*b) dados.")
    c=k*(a+b+1)+2*h*(a+b+1)
        println("Atualmente precisamos armazenar $c dados.")
    d=3*a*b -c
    if d>0
        println("A comprenssão resulta em $(d) dados economizados.")
    else
        println("A comprenssão resulta em $(-d) dados gastos a mais.")
    end
        println("Armazenamos aproximadamente $(round((c*100)/(3*a*b),digits=3))% da quantia de dados iniciais.")
    #descompressao
    UEx=expand(U)
    VEx=expand(V)
    imgR,imgG,imgB=YUV_to_RGB(Y,UEx,VEx)
    imgCom=colorview(RGB,imgR,imgG,imgB)
    MSR=quality(img,channelview(imgCom))
    println("Qualidade perdida da imagem é $MSR (0=sem perda, 1=máximo).")
    return imgCom
end
