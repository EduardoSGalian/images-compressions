"""
    quality(O,C)
Use the Mean Square Error metric to calculate the noise.

0 is true image and 1 is maximum result.
# Argument
- `O::Array{Float64,3}`: original image in tensor form
- `C::Array{Float64,3}`: compressed image in tensor form
"""
function quality(IMA,GES)
    q,m,n=size(IMA)
    q,mb,nb=size(GES)
    if m==mb && n==nb
        e=m*n
        t=zeros(3)
        r=zeros(3)
        for k=1:3
            for i=1:m
                for j=1:n
                    t[k]=t[k] + (IMA[k,i,j]-GES[k,i,j])^2
                end
            end
        end
        for h=1:3
            r[h]=t[h]/e
        end
        R=0
        for o=1:3
        R=R+r[o]/3
        end
        return R
    else
        println("The images have different sizes.")
    end
end
