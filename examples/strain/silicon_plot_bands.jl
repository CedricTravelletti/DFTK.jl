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
hgh_lda_family = artifact"hgh_pbe_hgh"
psp_hgh = joinpath(hgh_lda_family, "si-q12.hgh")

positions = [ones(3)/8, -ones(3)/8]
atoms = fill(ElementPsp(:Si; psp=load_psp(psp_hgh)), 2)
system = periodic_system(lattice, atoms, positions)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:gga_x_pbe, :gga_c_pbe], temperature = 1e-6)
basis_kwargs = (; kgrid = [7, 7, 7], Ecut = 70.0)
scf_kwargs = (; tol = 1e-7)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

AtomsCalculators.potential_energy(system, calculator)
# Construct 2D path through Brillouin zone

plot_bandstructure(calculator.state.scfres; kline_density=10)
