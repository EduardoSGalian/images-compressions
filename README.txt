Compressão de imagens

 Via SVD
  Arquivo "svd_compress.jl"
   svd_compress("nome da imagem","porcentagem de valores singulares a ser mantido")
  Arquivo "svd_results.jl", requer "quality.jl"
   svd_results("nome da imagem","porcentagem de valores singulares a ser mantido")
  
 Via SVD+
  Arquivo "svdE_results.jl", requer "quality.jl"
   svd_results("nome da imagem","porcentagem de valores singulares a ser mantido")
 
 Via JPG
  Arquivo "jpeg_compress.jl", requer "matrixsQ.jl"
   jpeg_compress("nome da imagem","matriz Q")
  Arquivo "jpeg_results.jl", requer "jpeg_compress.jl"
   jpeg_results("nome da imagem","matriz Q")
   
