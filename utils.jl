@inline second(v) = v[eachindex(v)[2]]
@inline safe_acos(x::T) where T = x≈1 ? zero(x) : x≈-1 ? T(π) : acos(x)
