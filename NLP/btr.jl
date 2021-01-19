using LinearAlgebra

struct BasicTrustRegion{T <: Real}
    η1:: T
    η2:: T
    γ1:: T
    γ2:: T
end

function BTRDefaults()
    return BasicTrustRegion(0.01,0.9,0.5,0.5)
end

mutable struct BTRState

    iter::Int
    x::Vector
    xcand::Vector
    g::Vector
    step::Vector
    Δ::Float64
    ρ::Float64

    tol::Float64

    trace
    keepTrace::Bool
    
    function BTRState()
        state = new()
        state.tol = 1e-6
        state.keepTrace = false
        return state
    end
end

function acceptCandidate!(state::BTRState, b::BasicTrustRegion)
    # If the iteration is successful, update the iterate
    if (state.ρ >= b.η1)
        return true
    else
        return false
    end
end

function updateRadius!(state::BTRState, b::BasicTrustRegion)
    if (state.ρ >= b.η2)
        stepnorm = norm(state.step)
        state.Δ = min(1e20,max(4*stepnorm,state.Δ))
    elseif (state.ρ >= b.η1)
        state.Δ *= b.γ2
    else
        state.Δ *= b.γ1
    end
end

function btr(f::Function, g!::Function, H!::Function, Step::Function,
    x0::Vector, state:: BTRState = BTRState(), ApproxH::Bool = false, verbose::Bool = false)
    
    b = BTRDefaults()
    state.iter = 0
    state.x = x0
    n=length(x0)

    tol2 = state.tol*state.tol
    
    state.g = zeros(n)
    # A better initialization procedure should be used with quasi-Newton approximations
    # We could rely on some preconditioner.
    H = zeros(n,n)+I
    
    fx = f(x0)
    g!(x0, state.g)
    state.Δ = 0.1*norm(state.g) # 1.0

    if (ApproxH)
        y = zeros(n)
        gcand = zeros(n)
        # H!(H, y, state.step)
    else
        H!(x0, H)
    end
    
    nmax = 100000
    if (state.keepTrace)
        state.trace= x0'
    end
    
    function model(s::Vector, g::Vector, H::Matrix)
        return dot(s, g)+0.5*dot(s, H*s)
    end
    
    while (dot(state.g,state.g) > tol2 && state.iter < nmax)
        # Compute the step by approximately minimize the model
        state.step = Step(state.g, H, state.Δ)
        state.xcand = state.x+state.step
	println("$(state.iter). $(state.Δ) $(state.step) $(state.g)")

        # Compute the actual reduction over the predicted reduction
        fcand = f(state.xcand)
        state.ρ = (fcand-fx)/(model(state.step, state.g, H))

        if (ApproxH)
            g!(state.xcand, gcand)
            y = gcand-state.g;
            H = H!(H, y, state.step)
        end

        if (acceptCandidate!(state, b))
            state.x = copy(state.xcand)
            if (ApproxH == false)
                g!(state.x, state.g)
                H!(state.x, H)
            else
                state.g = copy(gcand)
            end
            fx = fcand
        end

        if (state.keepTrace)
            state.trace= [state.trace ; state.x']
        end
        
        updateRadius!(state, b)
        state.iter += 1
    end
    
    return state
end

function CauchyStep(g::Vector, H::Matrix, Δ::Float64)
    q = dot(g,H*g)
    normg = norm(g)
    if (q <= 0)
        τ = 1.0
    else
        τ = min((normg*normg*normg)/(q*Δ),1.0)
    end
    return -τ*g*Δ/normg
end

function BFGSUpdate(B::Matrix, y::Vector, s::Vector)
    Bs = B*s
    return B - (Bs*Bs')/dot(s, Bs) + (y*y')/dot(s,y)
end

function stopCG(normg::Float64, normg0::Float64, k::Int, kmax::Int, χ::Float64 = 0.1, θ::Float64 = 0.5)
    if ((k == kmax) || (normg <= normg0*min(χ, normg0^θ)))
        return true
    else
        return false
    end
end

function TruncatedCG(g::Vector, H::Matrix, Δ::Float64)
    n = length(g)
    s = zeros(n)

    normg0 = norm(g)
    v = g
    d = -v
    gv = dot(g,v)
    norm2d = gv
    norm2s = 0
    sMd = 0
    k = 0
    Δ2 = Δ*Δ

    while (stopCG(norm(g), normg0, k, n) == false)
        Hd = H*d
        κ = dot(d,Hd)
 
        # Is the curvature negative in the direction d?
        if (κ <= 0)
            σ = (-sMd+sqrt(sMd*sMd+norm2d*(Δ2-dot(s,s))))/norm2d
            s += σ*d
            break
        end

        α = gv/κ

        # Check is the model minimizer is outside the trust region
        norm2s += α*(2*sMd+α*norm2d)
        if (norm2s >= Δ2)
            σ = (-sMd+sqrt(sMd*sMd+norm2d*(Δ2-dot(s,s))))/norm2d
            s += σ*d
            break
        end

        # The model minimizer is inside the trust region
        s += α*d
        g += α*Hd
        v = g
        newgv = dot(g,v)
        β = newgv/gv
        gv = newgv
        d = -v+β*d
        sMd = β*(sMd+α*norm2d)
        norm2d = gv+β*β*norm2d
        
        k += 1;
    end
    
    return s
end

defaultState = BTRState()
