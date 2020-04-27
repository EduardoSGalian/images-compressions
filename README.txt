Uso "compressao_svd.jl"
 -Inclua em seu terminal
 -Use a função: compress_svd("nome da imagem","porcentagem de valores singulares a ser mantido")
 	+Uso da função svd_results() da mesma forma, porém precisa da função quality() presente no arquivo "compresao_jpeg.jl"

Uso "compressao_jpeg.jl"
 -Inclua em seu terminal, juntamente com "Qmatrixs.jl"
 -Apenas para ver a imagem comprimida use a função: final("nome da imagem")
 -Para análises use: jpeg_results("nome da imagem")
 -A tabela de quantização pode ser alterada alterando a matriz 8x8 Q, no arquivo Qmatrixs.jl existe 4 sugestões.