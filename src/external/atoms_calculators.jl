# Define AtomsCalculators interface for DFTK.
#
# This interface is inspired by the one used in Molly.jl,
# see https://github.com/JuliaMolSim/Molly.jl/blob/master/src/types.jl
using AtomsBase
using AtomsCalculators
import AtomsCalculators: calculator_state, update_state
using Unitful
using UnitfulAtomic


Base.@kwdef struct DFTKParameters
    model_kwargs = (; )
    basis_kwargs = (; )
    scf_kwargs   = (; )
end

empty_state() = (; ψ=nothing, ρ=nothing)

mutable struct DFTKCalculator
    ps::DFTKParameters
    st
end

"""
Updates state between calculations. This can be used to speed up computations 
by interpolating from past results. 

Currently, it re-uses the last density as a starting point for the scf if 
the lattice hasn't changed, and starts from scratch otherwise.
"""
function update_state(algorithm, st, original_system, new_system)
    if bounding_box(original_system) != bounding_box(new_system)
        empty_state()
    else
        st
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
function DFTKCalculator(ps::DFTKParameters)
    DFTKCalculator(ps, empty_state())  # Create dummy state if not given.
end

function DFTKCalculator(; verbose=false, model_kwargs, basis_kwargs, scf_kwargs)
    if !verbose
        scf_kwargs = merge(scf_kwargs, (; callback=identity))
    end
    ps = DFTKParameters(; model_kwargs, basis_kwargs, scf_kwargs)
    DFTKCalculator(ps)
end

function compute_scf!(system::AbstractSystem, calculator::DFTKCalculator, ps::DFTKParameters, st)
    # We re-use the symmetries from the state to avoid issues
    # with accidentally more symmetric structures.
    symmetries = haskey(st.scfres, :basis) ? st.scfres.basis.model.symmetries : true
    model = model_DFT(system; symmetries, ps.model_kwargs...)
    basis = PlaneWaveBasis(model; ps.basis_kwargs...)

    ρ = @something st.scfres.ρ guess_density(basis, system)
    scfres = self_consistent_field(basis; ρ, st.scfres.ψ, ps.scf_kwargs...)
end

AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
        ::AtomsCalculators.Energy, system::AbstractSystem, calculator::DFTKCalculator,
        ps=nothing, st=nothing; kwargs...)
    scfres = compute_scf!(system, calculator, ps, st)
    return (; :energy => scfres.energies.total * u"hartree", :state => scfres)
end

AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
        ::AtomsCalculators.Forces, system::AbstractSystem, calculator::DFTKCalculator,
        ps=nothing, st=nothing; kwargs...)
    scfres = compute_scf!(system, calculator, ps, st)
    return (; :forces => compute_forces_cart(scfres) * u"hartree/bohr", :state => scfres)
end

AtomsCalculators.@generate_interface function AtomsCalculators.calculate(
        ::AtomsCalculators.Virial, system::AbstractSystem, calculator::DFTKCalculator,
        ps=nothing, st=nothing; kwargs...)
    scfres = compute_scf!(system, calculator, ps, st)
    virial = - (compute_stresses_cart(scfres) * scfres.basis.model.unit_cell_volume) * u"hartree"
    return (; :virial => virial, :state => scfres)
end
