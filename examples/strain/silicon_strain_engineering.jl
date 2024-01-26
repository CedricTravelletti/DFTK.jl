# Experiment effects of strain on bandgap on silicon.
using DFTK
using AtomsBase
using AtomsCalculators
using LazyArtifacts
using Unitful
using UnitfulAtomic
using GeometryOptimization
    
    
lattice = [0.0  5.131570667152971 5.131570667152971;
           5.131570667152971 0.0 5.131570667152971;
           5.131570667152971 5.131570667152971  0.0]
hgh_lda_family = artifact"hgh_lda_hgh"
psp_hgh = joinpath(hgh_lda_family, "si-q4.hgh")

positions = [ones(3)/8, -ones(3)/8]
atoms = fill(ElementPsp(:Si; psp=load_psp(psp_hgh)), 2)
system = periodic_system(lattice, atoms, positions)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-6)
basis_kwargs = (; kgrid = [5, 5, 5], Ecut = 20.0)
scf_kwargs = (; tol = 1e-6)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)


function strain_energy(strain)
    new_system = apply_voigt_strain(system, strain)
    # Initialize empty state, since we change the cell.
    state = DFTK.DFTKState()
    AtomsCalculators.potential_energy(new_system, calculator; state)
end

function strain_indirect_band_gap(strain)
    new_system = apply_voigt_strain(system, strain)
    # Initialize empty state, since we change the cell.
    state = DFTK.DFTKState()
    # Comput energy so we have the scfres.
    AtomsCalculators.potential_energy(new_system, calculator; state)
    compute_band_gaps(calculator.state.scfres)[:indirect_bandgap]
end
