using Images, ImageMagick, LinearAlgebra

function svd_compress(i,p) #i é o nome da imagem, p é um valor relacionado com a quantia de valores singulares da imagem 0<p<100 %
    img=channelview(load(i)) #carrega a imagem
    UR,SR,VR=svd(img[1,:,:]) #calcula SVD da imagem vermelha
    UG,SG,VG=svd(img[2,:,:]) #calcula SVD da imagem verde
    UB,SB,VB=svd(img[3,:,:]) #calcula SVD da imagem azul
    k=Int(ceil(length(SR)*p/100)) #transforma p em %; k definido de acordo com quantos valores singulares existem
    imgR=UR[:,1:k]*diagm(0=>SR[1:k])*transpose(VR[:,1:k]) # criação da imagem nova vermelha
    imgG=UG[:,1:k]*diagm(0=>SG[1:k])*transpose(VG[:,1:k]) # criação da imagem nova verde
    imgB=UB[:,1:k]*diagm(0=>SB[1:k])*transpose(VB[:,1:k]) # criação da imagem nova azul
    return colorview(RGB,imgR,imgG,imgB)
end
