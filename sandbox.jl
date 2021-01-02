using Random
using StatsBase
using BenchmarkTools

# Initializations
T = 1000;                           # sample size
s = zeros(Int, 2);                  # anchors
c = [1 + (i - 1) % 4 for i in 1:T]; # vector of cluster indicators
s̄_clu = zeros(Int, T);              # s̄ (cluster indicators)
s̄_idx = zeros(Int, T);              # s̄ (indices)
σ = collect(1:T);

# Update s
function update_s!(s, T)
    StatsBase.sample!(1:T, s, replace=false)
end
update_s!(s, T);

# Update s̄
function update_s̄!(s̄_clu, s̄_idx, s, c)
    card_s̄ = 0
    T = length(c)
    resize!(s̄_clu, T)
    resize!(s̄_idx, T)
    c1, c2 = c[s[1]], c[s[2]]
    for i ∈ 1:T
        i ∈ s && continue
        (c[i] == c1) && (card_s̄ += 1; s̄_clu[card_s̄] = 1; s̄_idx[card_s̄] = i)
        (c[i] == c2) && (card_s̄ += 1; s̄_clu[card_s̄] = 2; s̄_idx[card_s̄] = i)
    end
    card_s̄ += 2
    s̄_clu[card_s̄ - 1] = c1; s̄_idx[card_s̄ - 1] = s[1]
    s̄_clu[card_s̄ - 0] = c2; s̄_idx[card_s̄ - 0] = s[2]
    resize!(s̄_clu, card_s̄)
    resize!(s̄_idx, card_s̄)
    return card_s̄
end
card_s̄ = update_s̄!(s̄_clu, s̄_idx, s, c);

# Update σ
function update_σ!(σ, card_s̄)
    resize!(σ, card_s̄ - 2);
    randperm!(σ);
    resize!(σ, card_s̄);
    σ[card_s̄ - 1] = card_s̄ - 1;
    σ[card_s̄ - 0] = card_s̄ - 0;
    return σ
end
update_σ!(σ, card_s̄);

# Tests
@btime update_s!($s, $T);
@btime update_s̄!($s̄_clu, $s̄_idx, $s, $c);
@btime update_σ!($σ, $card_s̄);
