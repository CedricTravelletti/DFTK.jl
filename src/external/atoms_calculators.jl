# Define AtomsCalculators interface for DFTK.
#
# This interface is inspired by the one used in Molly.jl,
# see https://github.com/JuliaMolSim/Molly.jl/blob/master/src/types.jl
using AtomsBase
using AtomsCalculators


Base.@kwdef struct DFTKParameters
    model_kwargs = (; )
    basis_kwargs = (; )
    scf_kwargs   = (; )
end

mutable struct DFTKCalculator
    parameters::DFTKParameters
end

"""
Updates state between calculations. This can be used to speed up computations 
by interpolating from past results. 

Currently, it re-uses the last density as a starting point for the scf if 
the lattice hasn't changed, and starts from scratch otherwise.
"""
function update_state(algorithm, state, original_system, new_system)
    if bounding_box(original_system) != bounding_box(new_system)
        nothing
    else
        state
    end
end

"""
Construct a [AtomsCalculators](https://github.com/JuliaMolSim/AtomsCalculators.jl)
compatible calculator for DFTK. The `model_kwargs` are passed onto the
[`Model`](@ref) constructor, the `basis_kwargs` to the [`PlaneWaveBasis`](@ref)
constructor, the `scf_kwargs` to [`self_consistent_field`](@ref). At the very
least the DFT `functionals` and the `Ecut` needs to be specified.

By default the calculator preserves the symmetries that are stored inside the
`state` (the basis is re-built, but symmetries are fixed and not re-computed).

## Example
```julia-repl
julia> DFTKCalculator(; model_kwargs=(; functionals=[:lda_x, :lda_c_vwn]),
                        basis_kwargs=(; Ecut=10, kgrid=(2, 2, 2)),
                        scf_kwargs=(; tol=1e-4))
```
"""
function DFTKCalculator(parameters::DFTKParameters)
    DFTKCalculator(parameters, DFTKState())  # Create dummy state if not given.
end

function DFTKCalculator(; verbose=false, model_kwargs, basis_kwargs, scf_kwargs)
    if !verbose
        scf_kwargs = merge(scf_kwargs, (; callback=identity))
    end
    parameters = DFTKParameters(; model_kwargs, basis_kwargs, scf_kwargs)
    DFTKCalculator(parameters)
end

function compute_scf!(system::AbstractSystem, calculator::DFTKCalculator, state::DFTKState)
    parameters = calculator.parameters

    # We re-use the symmetries from the state to avoid issues
    # with accidentally more symmetric structures.
    symmetries = haskey(state.scfres, :basis) ? state.scfres.basis.model.symmetries : true
    model = model_DFT(system; symmetries, parameters.model_kwargs...)
    basis = PlaneWaveBasis(model; parameters.basis_kwargs...)

    ρ = @something state.scfres.ρ guess_density(basis, system)
    self_consistent_field(basis; ρ, state.scfres.ψ, parameters.scf_kwargs...)
end

AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
        ::AtomsCalculators.Energy,
        system::AbstractSystem, calculator::DFTKCalculator,
        parameters, state;
        kwargs...)
    scfres = compute_scf!(system, calculator, state)
    return (; :energy => scfres.energies.total,
              :state => scfres)
end

AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
        system::AbstractSystem, calculator::DFTKCalculator,
        parameters, state;
        kwargs...)
    scfres = compute_scf!(system, calculator, state)
    return (; :forces => compute_forces_cart(scfres),
              :state => scfres)
end
