"""
    Return an offset vector `p`, where p[i] is the prior 
    probability of any hypothesis such that sum_{j=1}^m γj = i,
    according to Womack (2015)'s proposal
"""
function ph0(m, ζ)
    p = OffsetArray([zeros(m); 1.0], 0:m)
    for l = m-1:-1:0
        for j = 1:m-l
            p[l] += ζ * p[l + j] * binomial(l + j, l)
        end
    end
    p /= sum(p)
    for l = 0:m-1
        p[l] /= binomial(m, l)
    end
    return p    
end

function γcode(γ)
    code = 1
    @inbounds for j = 2:length(γ)
        code += γ[j] * 2^(j - 2)
    end
    code
end

function γvector(J, code)
    γ = zeros(Int, J)
    code -= 1
    for j = J:-1:2
        if (code ÷ 2^(j-2)) > 0
            code -= 2^(j-2)
            γ[j] = 1
        end
    end
    return γ
end

γstr(J, code) = string(γvector(J, code)[2:J]...)