using Images, ImageMagick, LinearAlgebra
include("quality.jl")

function svd_results(i,p)
    img=channelview(load(i)) #carrega a imagem
    UR,SR,VR=svd(img[1,:,:]) #calcula SVD da imagem vermelha
    UG,SG,VG=svd(img[2,:,:]) #calcula SVD da imagem verde
    UB,SB,VB=svd(img[3,:,:]) #calcula SVD da imagem azul
    k=Int(ceil(length(SR)*p/100)) #transforma p em %; k definido de acordo com quantos valores singulares existem
    imgR=UR[:,1:k]*diagm(0=>SR[1:k])*transpose(VR[:,1:k]) # criação da imagem nova vermelha
    imgG=UG[:,1:k]*diagm(0=>SG[1:k])*transpose(VG[:,1:k]) # criação da imagem nova verde
    imgB=UB[:,1:k]*diagm(0=>SB[1:k])*transpose(VB[:,1:k]) # criação da imagem nova azul
    vav,a,b=size(img)
        println("Anteriormente precisavamos armazenar $(3*a*b) dados.")
    c=3*k*(a+b+1)
        println("Atualmente precisamos armazenar $c dados.")
    d=3*a*b -c
    if d>0
        println("A comprenssão resulta em $(d) dados economizados.")
    else
        println("A comprenssão resulta em $(-d) dados gastos a mais.")
    end
        println("Armazenamos aproximadamente $(round((c*100)/(3*a*b),digits=3))% da quantia de dados iniciais.")
    imgCom=colorview(RGB,imgR,imgG,imgB)
    MSR=quality(img,channelview(imgCom))
    println("Qualidade perdida da imagem é $MSR (0=sem perda, 1=máximo).")
    return imgCom
end
