module L1Homotopy

export L1_homotopy

"""
    L1_homotopy(A, y, λₑ)

Solve the optimization problem

  minimize  J(x) = ‖A*x - y‖₂² + λₑ‖x‖₁
     x

  subject to  xᵢ > 0 ∀ i

using the homotopy algorithm.  Here we largely follow the notation from the
paper "Fast solution to L1-norm minimization problems when the solution may be
sparse" by David Donoho and Yaakov Tsaig.  The formulation above is labelled
(D_λ) in the paper, but there are several other equivalent formulations,
including two constrained formulations labelled (P_1) and (L_q) in the paper.
L_q is popular for statistical model selection and is called Lasso in the
statistical literature.


# Devdocs

The algorithm works by following the optimal solution from λ=∞ down to λ=λₑ.
It turns out this path is a series of linear segments x0 -> x1 -> ... -> xf
through the space of x's, with `x0=zeros(N)`.
"""
function L1_homotopy(A, y, λₑ)
    # TODO: ATA may be precached, or computed on the fly.
    ATA = Symmetric(A'*A)
    N = size(A,2)
    # Residual correlations of signal with system model,
    #    c(x) ≡ A'*(y - A*x)
    # initialized where x=0
    c = A'*y
    # Initial active set
    λ = maximum(c)
    Iₐ = find(c .== λ)
    ET = eltype(A)
    T = typeof((zero(ET)*one(ET) + zero(ET)) / one(ET))
    xₐ = zeros(Iₐ, T)
    iter = 0
    while λ > λₑ
        # Cache inactive indices
        Iᵢ = setdiff(1:length(y), Iₐ)
        # Update direction, restricted to the active set
        dₐ = Symmetric(ATA[Iₐ,Iₐ]) \ sign.(c[Iₐ])
        # Compute update magnitude, γ. There's two ways this can occur
        #
        # 1) An inactive set correlation residual increases to the current λ,
        #    causing an inactive element to enter the active set
        γ⁺ = Inf
        new_active = 0
        if !isempty(Iᵢ)
            dcᵢ_dγ = ATA[Iᵢ,Iₐ]*dₐ
            γs = (λ .- c[Iᵢ]) ./ (one(T) .- dcᵢ_dγ)
            for (i, γ) ∈ zip(Iᵢ, γs)
                if γ > 0 && γ <= γ⁺
                    # TODO: Multiple additions to active set
                    γ⁺ = γ
                    new_active = i
                end
            end
        end
        # 2) An active set element changes sign (ie, decreases to zero in the
        #    positively constrained formulation)
        # Active element to leave the active set
        γ⁻ = Inf
        old_active_ind_ind = 0
        if !isempty(Iₐ)
            γs = -xₐ ./ dₐ
            for (i, γ) ∈ enumerate(γs)
                if γ > 0 && γ <= γ⁻
                    # TODO: Mutliple removals from active set
                    γ⁻ = γ
                    old_active_ind_ind = i
                end
            end
        end
        γ = min(γ⁺, γ⁻)
        # Incremental update to  c ≡ A'*(y - A*x)
        c .-= γ .* ATA[:,Iₐ]*dₐ
        xₐ .+= γ .* dₐ
        λ = maximum(c)
        if γ⁺ < γ⁻
            push!(Iₐ, new_active)
            push!(xₐ, zero(T))
        else
            splice!(Iₐ, old_active_ind_ind)
            splice!(xₐ, old_active_ind_ind)
        end
        iter += 1
    end
    x = zeros(N)
    x[Iₐ] = xₐ
    x
end


end
