@testitem "Test band gaps calculations" setup=[TestCases] begin
    using DFTK
    using AtomsBase
    using AtomsCalculators
    using Unitful
    using UnitfulAtomic
    silicon = TestCases.silicon

    # Perturb from equilibrium so forces are not 0.
    system = periodic_system(silicon.lattice, silicon.atoms, silicon.positions)

    model_kwargs = (; temperature=1e-6, functionals=[:lda_x, :lda_c_pw])
    basis_kwargs = (; kgrid=[6, 6, 6], Ecut=30.0)
    scf_kwargs = (; tol=1e-6)
    model = model_DFT(system; model_kwargs...)
    basis = PlaneWaveBasis(model; basis_kwargs...)
    
    scfres = self_consistent_field(basis; scf_kwargs...)
    band_gaps = compute_band_gaps(scfres)
    
    # Reference values.
    vi = 4; εMax_valence = 0.258363u"hartree"; εMin_conduction = 0.2815352u"hartree"
    direct_bandgap = 0.0937999u"hartree"
    
    @test band_gaps[:εMax_valence] ≈ εMax_valence rtol=1e-3
    @test band_gaps[:εMin_conduction] ≈ εMin_conduction rtol=1e-3
    @test band_gaps[:direct_bandgap] ≈ direct_bandgap rtol=1e-3
    @test band_gaps[:valence_band_index] == vi
end
