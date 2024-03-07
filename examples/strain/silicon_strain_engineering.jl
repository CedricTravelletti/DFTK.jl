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
    
    
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
                   [1 0 1.];  # specified column by column
                   [1 1 0.]];
hgh_lda_family = artifact"hgh_lda_hgh"
psp_hgh = joinpath(hgh_lda_family, "si-q4.hgh")

positions = [ones(3)/8, -ones(3)/8]
atoms = fill(ElementPsp(:Si; psp=load_psp(psp_hgh)), 2)
system = periodic_system(lattice, atoms, positions)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-6)
basis_kwargs = (; kgrid = [6, 6, 6], Ecut = 50.0)
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
    compute_band_gaps(calculator.state.scfres)[:indirect_bandgap]
end
