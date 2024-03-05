# Experiment effects of strain on bandgap on silicon.
using DFTK
using AtomsBase
using AtomsCalculators
using LazyArtifacts
using Unitful
using UnitfulAtomic
using GeometryOptimization
using ComponentArrays
using ForwardDiff
    
    
lattice = [0.0  5.131570667152971 5.131570667152971;
           5.131570667152971 0.0 5.131570667152971;
           5.131570667152971 5.131570667152971  0.0]
hgh_lda_family = artifact"hgh_lda_hgh"
psp_hgh = joinpath(hgh_lda_family, "si-q4.hgh")

positions = [ones(3)/8, -ones(3)/8]
atoms = fill(ElementPsp(:Si; psp=load_psp(psp_hgh)), 2)
system = periodic_system(lattice, atoms, positions)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-3)
basis_kwargs = (; kgrid = [5, 5, 5], Ecut = 20.0)
scf_kwargs = (; tol = 1e-6)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)


function strain_energy(strain)
    # Only modify strains, not positions.
    positions = ComponentVector(atoms = position(system), strain = strain)
    new_system = update_positions(system, positions)
    
    AtomsCalculators.potential_energy(new_system, calculator)
end

function strain_indirect_band_gap(strain)
    # Only modify strains, not positions.
    positions = ComponentVector(atoms = position(system), strain = strain)
    new_system = update_positions(system, positions)
    
    # Comput energy so we have the scfres.
    AtomsCalculators.potential_energy(new_system, calculator)
    # compute_band_gaps(calculator.state.scfres)[:indirect_bandgap]
end

function f(x)
    new_system = update_positions(system, x * u"bohr")
    
    AtomsCalculators.potential_energy(new_system, calculator)
end

ForwardDiff.gradient(f, zeros(6))

using ForwardDiff: Dual
xin = [2000.0 + Dual(0, (1,0,0)), 20000.0 + Dual(0, (0,1,0)), 80.0]

f([
   [Dual{Float64}(0., 1.), Dual{Float64}(0., 1.), Dual{Float64}(0., 1.)],
   [Dual{Float64}(0., 1.), Dual{Float64}(0., 1.), Dual{Float64}(0., 1.)]
  ])
