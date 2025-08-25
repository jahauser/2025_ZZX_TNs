function ZZ_corrs(ρ::MPS)
    ρ /= doubledtrace(ρ)
    N = length(ρ)÷2
    sites = siteinds(ρ)
    output_matrix = Matrix{ComplexF64}(I, N, N)

    reduced_ρs = Vector{ITensor}()
    for i in 1:N
        bell_i = bell(sites[2i-1:2i])
        push!(reduced_ρs, bell_i[1]*bell_i[2]*ρ[2i-1]*ρ[2i])
    end

    left_stack = ITensor(1)
    for i in 1:N-1
        Z_i = M_bra(sites[2i-1:2i], PauliZ, 1)
        LS_Z = left_stack * Z_i[1]*Z_i[2]*ρ[2i-1]*ρ[2i]

        right_stack = ITensor(1)
        for j in N:-1:i+1
            Z_j = M_bra(sites[2j-1:2j], PauliZ, 1)
            Z_RS = Z_j[1]*Z_j[2]*ρ[2j-1]*ρ[2j] * right_stack
            
            output_matrix[i,j] = (reduce(*, [LS_Z, reduced_ρs[i+1:j-1]..., Z_RS])[])^2 # squared because kappa ea!?
            output_matrix[j,i] = output_matrix[i,j]

            right_stack = reduced_ρs[j] * right_stack
        end

        left_stack = left_stack * reduced_ρs[i]
    end

    return output_matrix
end

function κEA(ρ::MPS; ref=false)
    L = length(ρ)÷2 - ref
    ZZs = ZZ_corrs(ρ)
    return sum(ZZs[1:L,1:L] .^ 2) / L
end

function ZZZZ_corrs(ρ::MPS)
    ρ /= doubledtrace(ρ)
    ρ = combined(ρ)
    return correlation_matrix(ρ, "ZZ", "ZZ")
end

function κ2(ρ::MPS; ref=false)
    L = length(ρ)÷2 - ref

    corrs = ZZZZ_corrs(ρ)
    return sum(corrs[1:L,1:L])/L
end

function Ic2(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)

    S2_QR = -log2(inner(ρ, ρ))

    a, b = id_mps(sites[N-1], sites[N])
    
    ρ_Q = MPS(ρ[1:N-2])
    ρ_Q[N-2] *= ρ[N-1]*ρ[N]*a*b

    S2_Q = -log2(inner(ρ_Q, ρ_Q))
    return S2_Q - S2_QR
end

function SR(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)

    A = ITensor(1)
    for i in 1:2:N-2
        B, C = id_mps(sites[i],sites[i+1])
        A = A * B * C * ρ[i] * ρ[i+1]
    end
    ρ_R = A * ρ[N-1] * ρ[N]
    M = Matrix(ρ_R, inds(ρ_R))
    λs = abs.(eigvals(M))
    return sum([-λ*log2(λ) for λ in λs if abs(λ) > 1e-16])
end

function Svn(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    U, Λ, V = svd(ψ[b], linkind(ψ, b))

    ps = [Λ[i,i]^2 for i in 1:dim(Λ,1) if Λ[i,i] > 0]
    return -sum([p*log2(p) for p in ps])
end

function pure_SR(ψ::MPS)
    L = length(siteinds(ψ))-1
    return Svn(ψ, L)
end

function pure_κEA(ψ::MPS; ref=false)
    L = length(ψ) - ref
    corr = correlation_matrix(ψ, "Z", "Z")
    return sum(corr[1:L,1:L] .^ 2) / L
end

function boundary_ZZ_corrs(ρ::MPS)
    N = length(ρ)÷2
    sites = siteinds(ρ)
    # println(inner(MPS(sites, _ -> "+"), ρ))
    # ρ /= inner(MPS(sites, _ -> "+"), ρ)
    output_matrix = Matrix{ComplexF64}(I, N, N)

    denom = inner(MPS(sites, _ -> "+"), ρ)
    # println(denom)

    reduced_ρs = Vector{ITensor}()
    for i in 1:N
        plus_i = MPS(sites[2i-1:2i], _ -> "+")
        push!(reduced_ρs, plus_i[1]*plus_i[2]*ρ[2i-1]*ρ[2i])
    end

    left_stack = ITensor(1)
    for i in 1:N-1
        plus_i = MPS(sites[2i-1:2i], _ -> "+")
        Z_i = apply(op(PauliZ, sites[2i-1]), plus_i)
        # println(expect(Z_i, "X"))
        LS_Z = left_stack * Z_i[1]*Z_i[2]*ρ[2i-1]*ρ[2i]

        right_stack = ITensor(1)
        for j in N:-1:i+1
            plus_j = MPS(sites[2j-1:2j], _ -> "+")
            Z_j = apply(op(PauliZ, sites[2j-1]), plus_j)
            # println(expect(Z_j, "X"))
            Z_RS = Z_j[1]*Z_j[2]*ρ[2j-1]*ρ[2j] * right_stack
            
            # println( (reduce(*, [LS_Z, reduced_ρs[i+1:j-1]..., Z_RS])[])^2)
            output_matrix[i,j] = (reduce(*, [LS_Z, reduced_ρs[i+1:j-1]..., Z_RS])[])/denom
            output_matrix[j,i] = output_matrix[i,j]

            right_stack = reduced_ρs[j] * right_stack
        end

        left_stack = left_stack * reduced_ρs[i]
    end

    return output_matrix
end

function boundary_ZZZZ_corrs(ρ::MPS)
    N = length(ρ)÷2
    sites = siteinds(ρ)
    # ρ /= inner(MPS(sites, _ -> "+"), ρ)
    output_matrix = Matrix{ComplexF64}(I, N, N)
    denom = inner(MPS(sites, _ -> "+"), ρ)

    reduced_ρs = Vector{ITensor}()
    for i in 1:N
        plus_i = MPS(sites[2i-1:2i], _ -> "+")
        push!(reduced_ρs, plus_i[1]*plus_i[2]*ρ[2i-1]*ρ[2i])
    end

    left_stack = ITensor(1)
    for i in 1:N-1
        plus_i = MPS(sites[2i-1:2i], _ -> "+")
        ZZ_i = apply(op(kron(PauliZ,PauliZ), sites[2i-1:2i]), plus_i)
        LS_Z = left_stack * ZZ_i[1]*ZZ_i[2]*ρ[2i-1]*ρ[2i]

        right_stack = ITensor(1)
        for j in N:-1:i+1
            plus_j = MPS(sites[2j-1:2j], _ -> "+")
            ZZ_j = apply(op(kron(PauliZ,PauliZ), sites[2j-1:2j]), plus_j)
            Z_RS = ZZ_j[1]*ZZ_j[2]*ρ[2j-1]*ρ[2j] * right_stack
            
            output_matrix[i,j] = (reduce(*, [LS_Z, reduced_ρs[i+1:j-1]..., Z_RS])[])/denom
            output_matrix[j,i] = output_matrix[i,j]

            right_stack = reduced_ρs[j] * right_stack
        end

        left_stack = left_stack * reduced_ρs[i]
    end

    return output_matrix
end

function boundary_ZZ_susceptibility(ρ::MPS; ref=false)
    L = length(ρ)÷2 - ref

    corrs = boundary_ZZ_corrs(ρ)
    return sum(corrs[1:L,1:L])/L
end

function boundary_ZZZZ_susceptibility(ρ::MPS; ref=false)
    L = length(ρ)÷2 - ref

    corrs = boundary_ZZZZ_corrs(ρ)
    return sum(corrs[1:L,1:L])/L
end

function disorder_EA(ψ::MPS; ref=false)
    L = length(ψ) - ref
    κ = L
    for i in 1:L-1
        for j in i+1:L
            κ += 2expectX(ψ, i, j)^2
        end
    end
    return κ/L
end

function expectX(ψ::MPS, i::Int, j::Int)
    ϕ = apply(op(PauliX, siteinds(ψ)[i]), ψ)
    for site in i+1:j-1
        ϕ = apply(op(PauliX, siteinds(ψ)[site]), ϕ)
    end
    return inner(ϕ, ψ)/norm(ψ)^2
end