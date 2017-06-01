
export getLoop

function getLoop(;loc=[0.,0.,0.], r=0.5, dir="x")
    d = [-r, -r, r, r, -r]
    if dir == "x"
        L = [zeros(5)+loc[1] d[end:-1:1]+loc[2] d+loc[2]]
    elseif dir == "y"
        L = [d+loc[1] zeros(5)+loc[2] d[end:-1:1]+loc[2]]
    elseif dir == "z"
        L = [d[end:-1:1]+loc[1] d+loc[2] zeros(5)+loc[2]]
    end
    return L
end
