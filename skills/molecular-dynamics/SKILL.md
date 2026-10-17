---
name: molecular-dynamics
description: Runs and analyzes molecular dynamics simulations with OpenMM and MDAnalysis. Sets up protein/small molecule systems, defines force fields, runs energy minimization and production MD, and analyzes trajectories (RMSD, RMSF, contact maps, free energy surfaces). For structural biology, drug binding, and biophysics.
license: MIT
compatibility: Python 3.11+ with OpenMM and MDAnalysis; matplotlib for plots. Defaults to the CUDA platform and falls back to OpenCL then CPU when no GPU verifies, so the GPU path needs a CUDA-enabled OpenMM build and a current driver. PDBFixer and OpenFF are optional. Simulation and analysis run offline.
metadata:
  version: "1.4"
  last-reviewed: "2026-10-17"
  skill-author: Kuan-lin Huang
---

# Molecular Dynamics

## Overview

Molecular dynamics (MD) simulation computationally models the time evolution of molecular systems by integrating Newton's equations of motion. This skill covers two complementary tools:

- **OpenMM** (https://openmm.org/): High-performance MD simulation engine with GPU support, Python API, and flexible force field support
- **MDAnalysis** (https://mdanalysis.org/): Python library for reading, writing, and analyzing MD trajectories from all major simulation packages

**Reviewed versions:** OpenMM 8.6.1 and MDAnalysis 2.10.0. Reference-platform
smoke checks use small synthetic systems; GPU performance and scientific convergence
are not established by them. OpenMM latest API pages identify a development build,
so the installed 8.6.1 API is the executable reference. The platform-selection
script was reviewed against the 8.6.1 sources and exercised without OpenMM
installed and against a stand-in API, including a CUDA plugin that fails to
load. A Claude Science run on 2026-10-17 verified the CUDA platform on an 8 GiB consumer GPU under Windows with `openmm[cuda13]`, but only once the runtime DLL
directory was on the search path (see the Windows note below). Run
`python scripts/select_platform.py` and the installation test below once on the
target machine.

**Installation.** Install the GPU build — the CPU platform ships with it, so
there is no separate CPU-only variant to choose.

```bash
uv venv --python 3.12 .venv-md
uv pip install --python .venv-md/bin/python \
    "openmm[cuda13]==8.6.1" "MDAnalysis==2.10.0" matplotlib pandas
# Windows: .venv-md/Scripts/python.exe ; older CUDA 12 driver: "openmm[cuda12]==8.6.1"
# conda equivalent: conda install -c conda-forge "openmm=8.6.1" "cuda-version=12.9"
#   (win-64 openmm 8.6.1 is built against 12.9 and 13.0 only)
```

The CUDA extra brings its own NVRTC, which compiles kernels at runtime, so no
separate CUDA Toolkit is needed — only the system display driver. Plain
`openmm==8.6.1` gives OpenCL, CPU and Reference but no CUDA platform. Match the
CUDA series to the driver (`nvidia-smi` header, or conda's `__cuda` virtual
package); major releases are not binary compatible, so the pin is not cosmetic.
If OpenMM shares an environment with a CUDA torch build, match majors — a
`+cu132` torch beside `openmm[cuda12]` would install two parallel `nvidia-*` runtime
sets, wasting GiB and obscuring later breakage.

**Windows: the CUDA runtime DLLs are off the search path.** The pip `nvidia-*`
wheels (for `cuda13`, `site-packages/nvidia/cu13/bin/x86_64`) are not on the DLL
search path, so the CUDA plugin fails with error 126, never registers, and
OpenMM offers OpenCL as if nothing were wrong. Importing `select_platform` puts
those directories on the path, which works only if it is imported before the
first `import openmm` in the process; its `notes` list any plugin that failed to
load, so a quiet OpenCL fallback is visible in the report.

**Installing it under Claude Science.** Two obstacles, both measured 2026-10-17.
Host-side pip refuses `--index-url` and the `name[extra]` form outright, so
`openmm[cuda13]` must be installed by invoking pip as a subprocess inside the
sandbox. And a deep workspace path can break the install even when the download
succeeds: torch failed on a 279-character member path against a 260-character
MAX_PATH with long paths off, and needed extraction to a `\\?\`-prefixed
target. Keep the environment root short, or enable `LongPathsEnabled`. A
workspace-resident venv also may not survive a sweep — budget a reinstall.
`pdbfixer` failed to build on this host and was dropped; it is optional.

Validate with `python -m openmm.testInstallation`, which lists every registered
platform, computes forces on each and reports inter-platform differences. A
platform that is listed but errors there is a broken plugin, not a usable
device — rerun it after any driver update. On Windows with the pip CUDA wheels,
plain `python -m openmm.testInstallation` imports openmm without the DLL path
and reports no CUDA platform on a working install; run it from `scripts/` as
`python -c "import select_platform, runpy; runpy.run_module('openmm.testInstallation', run_name='__main__')"`
so the DLL directories are on the path first.

## Platform Selection: GPU by Default, CPU as Fallback

MD is the part of this skill that a GPU changes by orders of magnitude, so the
default is the fastest platform that is **shown to work**, not the fastest one
that is registered. OpenMM registers the CUDA platform whenever its plugin
library loads; that happens with a mismatched driver, a missing device and a
half-installed runtime alike. Registration is not a working device.

[scripts/select_platform.py](scripts/select_platform.py) resolves this. It walks
CUDA → OpenCL → CPU → Reference, and before accepting a platform it builds a
tiny bonded/Lennard-Jones system, computes forces and takes ten steps on it. A
platform that raises is skipped with its error recorded; the next one is tried.
Nothing raises out of the selector itself, including when OpenMM is absent.

```bash
python scripts/select_platform.py   # prints version, registered platforms, choice, notes
```

```python
from select_platform import select_platform, platform_report

choice = select_platform()          # PlatformChoice(name, properties, verified, notes)
print(platform_report(choice))      # e.g. "CUDA, precision=mixed, verified"
simulation = Simulation(modeller.topology, system, integrator,
                        **choice.simulation_kwargs())
```

- **Precision defaults to `mixed`** on CUDA/OpenCL, overriding OpenMM's `single`:
  forces in single, integration in double, which is the right trade on a GeForce
  card whose FP64 is reported at 1/64 of FP32 on Ada, while `single` alone
  degrades energy conservation over long integrations. Override with
  `SCI_PRECISION` or `precision=`. The OpenMM 8 key is `Precision`
  (`CudaPrecision` is a deprecated alias).
- **`SCI_PLATFORM` (or `SCI_DEVICE`) is honored verbatim** and not verified, so a
  forced platform fails loudly at Context creation instead of silently becoming
  something else; `SCI_PLATFORM=CPU` is the escape hatch.
- **Record platform and precision in provenance** — they change trajectories in
  detail, and checkpoints are platform-, hardware- and version-specific.
- **Verification cannot catch** a driver that hard-faults rather than raising, or
  anything size-dependent: a tiny system passing proves the platform runs, not
  that your system fits in VRAM or samples fast enough.

### Constraints on an 8 GB consumer card (8 GiB consumer card)

**One reference point, measured 2026-10-17:** ~23k-atom solvated lysozyme
(1AKI, AMBER14/TIP3P-FB, PME, 2 fs, HBonds constrained), NPT at 300 K/1 bar,
`Precision=mixed` — **more than an order of magnitude faster than the CPU platform**, or
far more sampling overnight. Throughput scales with atom count, cutoff and PME
grid, so measure your own system with `StateDataReporter(..., speed=True)` on an
idle desktop; the same run during a large download reported a measurably lower rate.

- **VRAM is rarely the MD bottleneck**: device state scales with atom count plus
  a PME grid, a small fraction of 8 GB for the 10^4–10^5-atom systems below, so
  throughput is what decides how much you can sample.
- **VRAM does bind** for million-atom systems, fine PME grids, large cutoffs or
  several concurrent Contexts; on Windows the desktop shares the same 8 GB
  (roughly 6 GiB usable), so check free memory with
  `nvidia-smi` first — the selector's tiny-system check cannot predict a failure
  at real Context creation.
- **Overflow does not announce itself.** A torch allocation past the card does
  not raise on Windows; WDDM pages it to host RAM and throughput collapses
  (measured 19–50× on a matmul). Expect a mysteriously slow run rather than an
  out-of-memory error, and suspect memory pressure when ns/day drops.
- **GeForce on Windows runs under WDDM** (no TCC), so a long-running kernel can
  trip the ~2 s display watchdog and reset the context; ordinary steps are far
  shorter, but minimization on a very large system is the realistic exposure,
  and a run dying with a device reset points here.
- **Driving the display from the same GPU** makes `speed=True` timings noisy, so
  treat busy-desktop numbers as a lower bound.

## When to Use This Skill

Use molecular dynamics when:

- **Protein stability analysis**: How does a mutation affect protein dynamics?
- **Drug binding simulations**: Characterize binding mode and residence time of a ligand
- **Conformational sampling**: Explore protein flexibility and conformational changes
- **Protein-protein interaction**: Model interface dynamics and binding energetics
- **RMSD/RMSF analysis**: Quantify structural fluctuations from a reference structure
- **Free energy estimation**: Compute binding free energy or conformational free energy
- **Membrane simulations**: Model proteins in lipid bilayers
- **Intrinsically disordered proteins**: Study IDR conformational ensembles

## Core Workflow: OpenMM Simulation

Examples require system-specific preparation and validation. Review biological
assembly, alternate locations, missing loops, termini, protonation, disulfides,
ligands and cofactors before parameterization. An MD trajectory alone does not
establish binding free energy or residence time.

The functions below form one Python module: execute their import blocks together.
Use a new output prefix for every stage to avoid overwriting earlier results.
`select_platform` is imported from this skill's `scripts/` directory; run from
there or add it to `sys.path`.

### 1. System Preparation

```python
from select_platform import select_platform, platform_report  # before openmm (Windows DLL path)
from openmm.app import *
from openmm import *
from openmm.unit import *

def prepare_system_from_pdb(pdb_file, forcefield_name="amber14-all.xml",
                              water_model="amber14/tip3pfb.xml", ph=7.0):
    """
    Prepare an OpenMM system from a PDB file.

    Args:
        pdb_file: Path to cleaned PDB file (use PDBFixer for raw PDB files)
        forcefield_name: Force field XML file
        water_model: Water model XML file

    Returns:
        modeller, system
    """
    # Load PDB
    pdb = PDBFile(pdb_file)

    # Load force field
    forcefield = ForceField(forcefield_name, water_model)

    # Add hydrogens and solvate
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=ph)

    # TIP3P geometry also serves TIP3P-FB; XML supplies the parameters.
    # 10 Å padding; added salt excludes neutralizing counterions.
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=10*angstroms,
        ionicStrength=0.15*molar, neutralize=True
    )

    print(f"System: {modeller.topology.getNumAtoms()} atoms, "
          f"{modeller.topology.getNumResidues()} residues")

    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,         # Particle Mesh Ewald for long-range electrostatics
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds,           # Constrain bonds involving H (not intermolecular H-bonds)
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )

    return modeller, system
```

### 2. Energy Minimization

```python
from openmm.app import *
from openmm import *
from openmm.unit import *

def minimize_energy(modeller, system, output_pdb="minimized.pdb",
                     max_iterations=1000, tolerance=10.0, platform_name=None):
    """
    Energy minimize the system to remove steric clashes.

    Args:
        modeller: Modeller object with topology and positions
        system: OpenMM System
        output_pdb: Path to save minimized structure
        max_iterations: Maximum minimization steps
        tolerance: Convergence criterion in kJ/mol/nm
        platform_name: Force a platform; None runs verified GPU-first selection

    Returns:
        simulation object with minimized positions
    """
    # Keep the integrator at 2 fs for the subsequent NVT/NPT examples.
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # GPU first, verified by actually running a tiny system, then CPU.
    # A registered GPU plugin does not prove a working device or driver.
    # An explicit platform_name is honored without verification and fails loudly.
    choice = select_platform(platform_name)
    print(f"Platform: {platform_report(choice)}")  # record this with the run
    simulation = Simulation(modeller.topology, system, integrator,
                            **choice.simulation_kwargs())
    simulation.context.setPositions(modeller.positions)

    # Check initial energy
    state = simulation.context.getState(getEnergy=True)
    print(f"Initial energy: {state.getPotentialEnergy()}")

    # Minimize
    simulation.minimizeEnergy(
        tolerance=tolerance*kilojoules_per_mole/nanometer,
        maxIterations=max_iterations
    )

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Minimized energy: {state.getPotentialEnergy()}")

    # Save minimized structure
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    return simulation
```

### 3. NVT Equilibration

```python
from openmm.app import *
from openmm import *
from openmm.unit import *

def run_nvt_equilibration(simulation, n_steps=50000, temperature=300,
                            report_interval=1000, output_prefix="nvt"):
    """
    NVT equilibration: constant N, V, T.
    Equilibrate velocities to target temperature.

    Args:
        simulation: OpenMM Simulation (after minimization)
        n_steps: Number of MD steps (50000 × 2fs = 100 ps)
        temperature: Temperature in Kelvin
        report_interval: Steps between data reports
        output_prefix: File prefix for trajectory and log
    """
    # This example is unrestrained; add validated restraints before Context creation.
    if any("Barostat" in type(f).__name__ and f.getFrequency() > 0
           for f in simulation.system.getForces()):
        raise ValueError("NVT requires all barostats disabled or absent")

    # Set both the thermostat target and initial velocities.
    simulation.integrator.setTemperature(temperature*kelvin)
    simulation.context.setVelocitiesToTemperature(temperature*kelvin)

    # Add reporters
    simulation.reporters = []

    # Log file
    simulation.reporters.append(
        StateDataReporter(
            f"{output_prefix}_log.txt",
            report_interval,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            speed=True
        )
    )

    # DCD trajectory (compact binary format)
    simulation.reporters.append(
        DCDReporter(f"{output_prefix}_traj.dcd", report_interval)
    )

    duration = (n_steps * simulation.integrator.getStepSize()).value_in_unit(picoseconds)
    print(f"Running NVT equilibration: {n_steps} steps ({duration:.1f} ps)")
    simulation.step(n_steps)
    print("NVT equilibration complete")

    return simulation
```

### 4. NPT Equilibration and Production

Call this function first with a dedicated NPT equilibration prefix. Inspect density,
energy, structure and replicate stability before a separate production call; the
default duration is an example, not an equilibration or convergence criterion.

```python
def run_npt_production(simulation, n_steps=500000, temperature=300, pressure=1.0,
                        report_interval=5000, output_prefix="npt"):
    """
    NPT production run: constant N, P, T.

    Args:
        n_steps: Production steps (500000 × 2fs = 1 ns)
        temperature: Temperature in Kelvin
        pressure: Pressure in bar
        report_interval: Steps between reports
    """
    # Keep the Langevin thermostat and barostat at the same temperature.
    simulation.integrator.setTemperature(temperature*kelvin)
    system = simulation.system
    barostats = [f for f in system.getForces() if "Barostat" in type(f).__name__]
    if not system.usesPeriodicBoundaryConditions():
        raise ValueError("NPT requires a periodic system")
    if not barostats:
        system.addForce(MonteCarloBarostat(pressure*bar, temperature*kelvin, 25))
        simulation.context.reinitialize(preserveState=True)
    elif len(barostats) != 1 or not isinstance(barostats[0], MonteCarloBarostat):
        raise ValueError("This example supports one isotropic MonteCarloBarostat")
    else:
        barostats[0].setDefaultPressure(pressure*bar)
        barostats[0].setDefaultTemperature(temperature*kelvin)
        barostats[0].setFrequency(25)
    # Existing Context parameters must also be updated on repeated calls.
    simulation.context.setParameter(MonteCarloBarostat.Pressure(), pressure)
    simulation.context.setParameter(MonteCarloBarostat.Temperature(), temperature)

    # Update reporters
    simulation.reporters = []
    simulation.reporters.append(
        StateDataReporter(
            f"{output_prefix}_log.txt",
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True
        )
    )
    simulation.reporters.append(
        DCDReporter(f"{output_prefix}_traj.dcd", report_interval)
    )

    # Save checkpoints
    simulation.reporters.append(
        CheckpointReporter(f"{output_prefix}_checkpoint.chk", 50000)
    )

    duration = (n_steps * simulation.integrator.getStepSize()).value_in_unit(nanoseconds)
    print(f"Running NPT stage: {n_steps} steps ({duration:.3f} ns)")
    simulation.step(n_steps)
    simulation.saveCheckpoint(f"{output_prefix}_checkpoint.chk")
    simulation.saveState(f"{output_prefix}_state.xml")
    print("NPT stage complete")
    return simulation
```

## Trajectory Analysis with MDAnalysis

### 1. Load Trajectory

Use the solvated topology with exactly the DCD atom count and order (e.g. the
`minimized.pdb` above), not the original unsolvated input. MDAnalysis converts
lengths to Å and time to ps by default; OpenMM bare coordinates use nm. Preserve
frame box vectors and actual timestamps; do not infer time from frame number.

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, contacts
import numpy as np
import matplotlib.pyplot as plt

def load_trajectory(topology_file, trajectory_file):
    """
    Load an MD trajectory with MDAnalysis.

    Args:
        topology_file: PDB, PSF, or other topology file
        trajectory_file: DCD, XTC, TRR, or other trajectory
    """
    u = mda.Universe(topology_file, trajectory_file)
    print(f"Universe: {u.atoms.n_atoms} atoms, {u.trajectory.n_frames} frames")
    first, last = u.trajectory[0].time, u.trajectory[-1].time
    print(f"Time range: {first:g} to {last:g} ps")
    u.trajectory[0]
    return u
```

### 2. RMSD Analysis

```python
def compute_rmsd(u, selection="backbone", reference_frame=0):
    """
    Compute RMSD of selected atoms relative to reference frame.

    Args:
        u: MDAnalysis Universe
        selection: Atom selection string (MDAnalysis syntax)
        reference_frame: Frame index for reference structure

    Returns:
        numpy array with columns [frame index, time in ps, RMSD in Angstroms]
    """
    # RMSD performs its own fit; avoid rotating the stored trajectory here.
    # Make the selected molecule whole before analysis of periodic trajectories.
    R = rms.RMSD(u, select=selection, ref_frame=reference_frame)
    R.run()

    rmsd_data = R.results.rmsd  # columns: frame, time, RMSD
    return rmsd_data

def plot_rmsd(rmsd_data, title="RMSD over time", output_file="rmsd.png"):
    """Plot RMSD over simulation time."""
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(rmsd_data[:, 1] / 1000, rmsd_data[:, 2], 'b-', linewidth=0.5)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(title)
    ax.axhline(rmsd_data[:, 2].mean(), color='r', linestyle='--',
               label=f'Mean: {rmsd_data[:, 2].mean():.2f} Å')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    return fig
```

### 3. RMSF Analysis (Per-Residue Flexibility)

Make molecules whole and align a separate structural-analysis copy before RMSF;
RMSF does not align. Do not run periodic contacts on a rotated trajectory whose
box was not rotated. See [analysis reference](references/mdanalysis_analysis.md)
for alignment, PCA, DSSP, hydrogen bonds and population-derived free energy.

```python
def compute_rmsf(u, selection="protein and name CA", start_frame=0):
    """
    Compute per-residue RMSF (flexibility).

    Returns:
        residue_keys, rmsf_values; keys are (resindex, segid, resid, resname)
    """
    # Select atoms
    atoms = u.select_atoms(selection)

    if not len(atoms):
        raise ValueError("RMSF selection is empty")
    # A multi-atom selection returns mean atomic RMSF, not COM RMSF.
    R = rms.RMSF(atoms)
    R.run(start=start_frame)

    # Average by residue
    resids = []
    rmsf_per_res = []
    for res in u.select_atoms(selection).residues:
        res_atoms = res.atoms.intersection(atoms)
        if len(res_atoms) > 0:
            resids.append((res.ix, res.segid, res.resid, res.resname))
            # RMSF is indexed within the selected AtomGroup, not the Universe.
            rmsf_per_res.append(R.results.rmsf[atoms.resindices == res.ix].mean())

    return resids, np.array(rmsf_per_res)
```

### 4. Protein-Ligand Contacts

```python
def analyze_contacts(u, protein_sel="protein", ligand_sel="resname LIG",
                      radius=4.5, start_frame=0, periodic=True):
    """
    Track protein-ligand contacts over trajectory.

    Args:
        radius: Atom-pair cutoff in Angstroms; any pair defines a residue contact.
        periodic: Require box data and use minimum-image distances.
    Returns sets of (resindex, segid, resid, resname), one per analyzed frame.
    """
    from MDAnalysis.lib.distances import distance_array

    protein = u.select_atoms(protein_sel)
    ligand = u.select_atoms(ligand_sel)

    if not len(protein) or not len(ligand):
        raise ValueError("Protein and ligand selections must both be nonempty")
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be finite and positive")
    contact_frames = []
    for ts in u.trajectory[start_frame:]:
        if periodic and (ts.dimensions is None or not np.all(np.isfinite(ts.dimensions))
                         or np.any(ts.dimensions[:3] <= 0)):
            raise ValueError("Periodic contacts require valid frame box dimensions")
        # This includes hydrogen atoms unless the caller selects heavy atoms.
        distances = contacts.contact_matrix(
            distance_array(protein.positions, ligand.positions, box=ts.dimensions if periodic else None),
            radius,
        )
        contact_residues = set()
        for i in range(distances.shape[0]):
            if distances[i].any():
                res = protein[i].residue
                contact_residues.add((res.ix, res.segid, res.resid, res.resname))
        contact_frames.append(contact_residues)

    return contact_frames
```

## Force Fields and Preparation

Choose a validated protein/ligand/water/ion combination for the scientific system.
The working example retains AMBER14/TIP3P-FB; it is not a universal recommendation.
Current OpenMM also bundles `amber19-all.xml` (ff19SB, DNA OL21, RNA OL3, lipid21),
and `charmm36_2024.xml` with its own `charmm36_2024/water.xml`. Generic TIP3P and
CHARMM-modified TIP3P are not interchangeable. IDP ensembles are especially
sensitive to protein-water balance. Four-site waters need extra particles and a
matching Modeller geometry, unlike the three-site example above.

AMBER protein XML does not parameterize arbitrary ligands. Use a reviewed GAFF2
or OpenFF route with explicit stereochemistry, protonation, bond orders, atom
mapping and charge method. A ligand-only Interchange is not a protein-ligand system.
See [preparation reference](references/system_preparation.md) for conservative
PDBFixer and OpenFF examples, installation requirements and current release details.

## Validation and Provenance

- Start with minimization; finite energy does not establish a valid structure.
- Use staged NVT → NPT equilibration → production. This example has no restraints.
- A 2 fs step with constrained H-containing bonds is a starting choice. Larger
  steps/HMR require integrator, stability and observable validation for that system.
- Use per-frame boxes for periodic distances. Make molecules whole for structural
  observables using trusted bond topology; unwrapping cannot infer missing bonds.
- Assess burn-in, autocorrelation, effective samples, replicas and uncertainties.
  Short stable temperature/density traces do not establish conformational sampling.
- Save topology/atom order, System/Integrator XML, force-field files/versions,
  seeds, platform/precision (`platform_report()`), parameters and analysis
  selections. Re-run `python -m openmm.testInstallation` after a driver or
  OpenMM update before continuing a campaign on the same card. Checkpoints are
  platform/hardware/version specific; XML states are more portable but do not
  retain all internal random-generator state. Test restarting the actual setup.

## Additional Resources

- **OpenMM documentation**: https://openmm.org/documentation.html
- **OpenMM installation and GPU setup**: https://docs.openmm.org/latest/userguide/application/01_getting_started.html
- **OpenMM platform-specific properties (Precision, DeviceIndex, Threads)**: https://docs.openmm.org/latest/userguide/library/04_platform_specifics.html
- **MDAnalysis user guide**: https://docs.mdanalysis.org/
- **GROMACS** (alternative MD engine): https://manual.gromacs.org/
- **NAMD** (alternative): https://www.ks.uiuc.edu/Research/namd/
- **CHARMM-GUI** (web-based system builder): https://charmm-gui.org/
- **AmberTools** (free Amber tools): https://ambermd.org/AmberTools.php
- **OpenMM paper**: Eastman P et al. (2017) PLOS Computational Biology. PMID: 28278240
- **MDAnalysis paper**: Michaud-Agrawal N et al. (2011) J Computational Chemistry. PMID: 21500218
