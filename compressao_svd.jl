using Images, ImageMagick, LinearAlgebra

function compress_svd(i,p) #i é o nome da imagem, p é um valor relacionado com a quantia de valores singulares da imagem 0<p<100 %
    img=channelview(load(i)) #carrega a imagem
    UR,SR,VR=svd(img[1,:,:]) #calcula SVD da imagem vermelha
    UG,SG,VG=svd(img[2,:,:]) #calcula SVD da imagem verde
    UB,SB,VB=svd(img[3,:,:]) #calcula SVD da imagem azul

    k=Int(ceil(length(SR)*p/100)) #transforma p em %; k definido de acordo com quantos valores singulares existem
    imgR=UR[:,1:k]*diagm(SR[1:k])*transpose(VR[:,1:k]) # criação da imagem nova vermelha
    imgG=UG[:,1:k]*diagm(SG[1:k])*transpose(VG[:,1:k]) # criação da imagem nova verde
    imgB=UB[:,1:k]*diagm(SB[1:k])*transpose(VB[:,1:k]) # criação da imagem nova azul

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

    return colorview(RGB,imgR,imgG,imgB)
end

function svd_results(image,percent)
    imgCom=compress_svd(image,percent)
    MSR=quality(channelview(load(image)),channelview(imgCom))
    println("Qualidade perdida da imagem é $MSR (0=sem perda, 1=máximo).")
end
