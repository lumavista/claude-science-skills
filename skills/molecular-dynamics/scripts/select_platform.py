"""Platform selection for OpenMM: GPU by default, CPU as the fallback.

This is the OpenMM sibling of `skills/esm/scripts/_device.py` and follows the
same contract: the CPU path always produces a correct result, the GPU is only
an accelerator, and nothing here raises when OpenMM is missing or when a GPU
plugin is registered but unusable -- the caller gets the CPU answer and keeps
going. An explicit SCI_PLATFORM (or SCI_DEVICE) is honored verbatim and is NOT
verified, so a forced platform fails loudly at Context creation rather than
silently running somewhere the caller did not ask for.

The one thing this module adds over its torch sibling: verification. OpenMM
registers the CUDA platform whenever the plugin shared library loads, which
says nothing about the driver, the GPU, or the runtime kernel compiler. A
registered GPU plugin does not prove a working device. So the default path
builds a tiny system, actually computes forces and takes a few steps on the
candidate platform, and only then accepts it; on failure it records the reason
and tries the next platform down the list.

Limits worth knowing before trusting a pass:
- Verification catches what OpenMM raises (a missing driver, an nvrtc failure,
  an out-of-memory Context). A driver that hard-faults the process cannot be
  caught in-process by any Python code; if that happens, set SCI_PLATFORM=CPU.
- A tiny system passing says the platform runs, not that a large system fits in
  VRAM, and says nothing about throughput.
- A GPU plugin that fails to *load* is never registered, so verification never
  sees it. Its load error is copied into `notes` instead, so a CUDA build that
  quietly became OpenCL is still visible to the caller.

Windows: the pip CUDA runtime wheels (`nvidia-*`, notably the CUDA 13 layout
under `site-packages/nvidia/cu13/bin/x86_64`) are not on the DLL search path,
so the CUDA plugin fails with error 126 and OpenMM drops to OpenCL. Importing
this module puts those directories on the search path, which only helps if it
happens before the first `import openmm` in the process -- OpenMM loads its
plugins at import time. Import this module first.
"""

from __future__ import annotations

import math
import os
import sys
import sysconfig
from dataclasses import dataclass, field

__all__ = [
    "PLATFORM_PREFERENCE",
    "PlatformChoice",
    "select_platform",
    "verify_platform",
    "platform_properties",
    "get_platform",
    "registered_platforms",
    "plugin_load_failures",
    "prepare_cuda_runtime",
    "platform_report",
]

# Fastest first. CPU is a real, optimized, multithreaded platform; Reference is
# the slow correctness implementation and is last for that reason.
PLATFORM_PREFERENCE = ("CUDA", "OpenCL", "CPU", "Reference")

# Platforms whose Precision property is meaningful (see platform_properties).
_PRECISION_PLATFORMS = frozenset({"CUDA", "OpenCL", "HIP"})

# Canonical spellings used when OpenMM is not importable and nothing can be
# matched against the live registry. OpenMM platform names are case sensitive.
_KNOWN_NAMES = ("CUDA", "OpenCL", "HIP", "CPU", "Reference")

# Handles from os.add_dll_directory; a directory leaves the search path when
# its handle is closed, so they are kept for the life of the process.
_DLL_DIR_HANDLES: list = []


def _cuda_runtime_dirs(roots) -> list[str]:
    """Directories under `<root>/nvidia` that contain DLLs, shallowest first."""
    found = []
    for root in roots:
        base = os.path.join(root, "nvidia")
        if not os.path.isdir(base):
            continue
        for dirpath, _dirnames, filenames in os.walk(base):
            if any(name.lower().endswith(".dll") for name in filenames):
                found.append(os.path.normpath(dirpath))
    return sorted(set(found), key=lambda d: (d.count(os.sep), d))


def prepare_cuda_runtime(roots=None) -> tuple[str, ...]:
    """Put pip-installed CUDA runtime DLLs on the Windows DLL search path.

    No-op off Windows. Returns every runtime directory found, whether or not it
    was already on PATH. Both mechanisms are used
    because OpenMM loads its plugins with plain LoadLibrary, which reads PATH,
    while Python extension modules honor os.add_dll_directory. Never raises.
    """
    if os.name != "nt":
        return ()
    if roots is None:
        paths = sysconfig.get_paths()
        roots = {paths.get("purelib"), paths.get("platlib")} - {None}
    try:
        dirs = _cuda_runtime_dirs(roots)
    except OSError:
        return ()
    current = os.environ.get("PATH", "").split(os.pathsep)
    added = [d for d in dirs if d not in current]
    if added:
        os.environ["PATH"] = os.pathsep.join(added + current)
    for d in dirs:
        try:
            _DLL_DIR_HANDLES.append(os.add_dll_directory(d))
        except (AttributeError, OSError):
            pass
    return tuple(dirs)


# Recorded before prepare_cuda_runtime runs: if openmm is already imported, its
# plugins have already failed or succeeded and the DLL path cannot help.
_OPENMM_IMPORTED_FIRST = "openmm" in sys.modules
_CUDA_RUNTIME_DIRS = prepare_cuda_runtime()


@dataclass(frozen=True)
class PlatformChoice:
    """The selected platform, the properties to create its Context with, and why.

    `verified` is True only when a tiny system actually ran on this platform in
    this process. It is False for an env-var override (honored unverified, by
    contract), for `verify=False`, and when OpenMM could not be imported.
    """

    name: str
    properties: dict[str, str] = field(default_factory=dict)
    verified: bool = False
    notes: tuple[str, ...] = ()

    def platform(self):
        """The `openmm.Platform` object, or None if OpenMM/the platform is absent."""
        return get_platform(self.name)

    def simulation_kwargs(self) -> dict:
        """Kwargs for `Simulation(topology, system, integrator, **kwargs)`.

        Empty when OpenMM is absent, so the call site fails at its own import
        rather than here. With OpenMM present this resolves the name strictly:
        a platform that does not exist -- a typo in SCI_PLATFORM, say -- raises
        here instead of quietly handing the run to OpenMM's own autoselection.
        """
        try:
            import openmm as mm
        except ImportError:
            return {}
        platform = mm.Platform.getPlatformByName(self.name)
        return {"platform": platform, "platformProperties": dict(self.properties)}


def registered_platforms() -> tuple[str, ...]:
    """Names of the platforms OpenMM has registered, fastest-reported first.

    Empty when OpenMM is not installed or its plugin loading fails. Registration
    only means the plugin loaded; use `verify_platform` before relying on one.
    """
    try:
        import openmm as mm
    except ImportError:
        return ()
    try:
        platforms = [mm.Platform.getPlatform(i) for i in range(mm.Platform.getNumPlatforms())]
        platforms.sort(key=lambda p: p.getSpeed(), reverse=True)
        return tuple(p.getName() for p in platforms)
    except Exception:
        return ()


def plugin_load_failures() -> tuple[str, ...]:
    """OpenMM's own record of plugin libraries that failed to load.

    A failed plugin never registers, so this is the only place a broken CUDA
    install shows up. Empty when OpenMM is missing or reports nothing.
    """
    try:
        import openmm as mm

        return tuple(str(f) for f in mm.Platform.getPluginLoadFailures())
    except Exception:
        return ()


def _load_notes() -> list[str]:
    """Notes explaining plugins that never registered, for PlatformChoice.notes."""
    notes = [f"plugin failed to load: {f}" for f in plugin_load_failures()]
    if notes and os.name == "nt":
        if _OPENMM_IMPORTED_FIRST:
            notes.append("openmm was imported before select_platform, so the CUDA "
                         "runtime DLL path was not prepared; import select_platform first")
        elif not _CUDA_RUNTIME_DIRS:
            notes.append("no pip CUDA runtime DLLs found under site-packages/nvidia; "
                         "install the openmm[cuda12] or openmm[cuda13] extra")
    return notes


def _canonical_name(requested: str) -> str:
    """Fix only the capitalization of a requested name; never substitute another.

    Matched case-insensitively against the live registry first, then against the
    static table, so "cuda" from a torch-flavored SCI_DEVICE reaches OpenMM as
    "CUDA". An unrecognized name is returned unchanged so it fails loudly.
    """
    for candidate in registered_platforms() or _KNOWN_NAMES:
        if candidate.lower() == requested.lower():
            return candidate
    return requested


def get_platform(name: str):
    """The named `openmm.Platform`, or None if OpenMM or that platform is absent."""
    try:
        import openmm as mm

        return mm.Platform.getPlatformByName(name)
    except Exception:
        return None


def platform_properties(name: str, *, precision: str | None = None) -> dict[str, str]:
    """Context properties for this platform; empty dict for CPU and Reference.

    Sets Precision='mixed' on CUDA/OpenCL/HIP, overriding OpenMM's 'single'
    default. 'mixed' computes forces in single precision and accumulates and
    integrates in double, which is the right default for MD on a consumer
    NVIDIA card: GeForce parts carry only token FP64 hardware (reported at 1/64
    of the FP32 rate on Ada, for instance), so 'double' is
    drastically slower there, while plain 'single' degrades energy conservation
    in long integrations. Override with the `precision` argument or SCI_PRECISION
    when a specific study needs 'single' speed or 'double' accuracy.

    Keys the platform does not advertise in `getPropertyNames()` are dropped
    rather than passed, since an unknown property makes Context creation raise.
    Note OpenMM 8 renamed these keys: the property is 'Precision', with
    'CudaPrecision' kept only as a deprecated alias.
    """
    if name not in _PRECISION_PLATFORMS:
        return {}
    chosen = (precision or os.environ.get("SCI_PRECISION") or "mixed").strip()
    properties = {"Precision": chosen}
    platform = get_platform(name)
    if platform is None:
        return properties
    try:
        supported = set(platform.getPropertyNames())
    except Exception:
        return properties
    return {k: v for k, v in properties.items() if k in supported}


def _tiny_system(mm):
    """A small bonded + Lennard-Jones system and its positions, in OpenMM units.

    Deliberately minimal: enough to force kernel compilation and a real force
    evaluation on the candidate platform, small enough that it costs nothing and
    cannot fail for lack of memory. Charges are zero so no Ewald/PME setup is
    implied; this is a liveness check, not an accuracy benchmark.
    """
    n_particles = 8
    system = mm.System()
    nonbonded = mm.NonbondedForce()
    bonds = mm.HarmonicBondForce()
    positions = []
    for i in range(n_particles):
        system.addParticle(12.0)  # amu
        nonbonded.addParticle(0.0, 0.3, 0.5)  # charge e, sigma nm, epsilon kJ/mol
        positions.append(mm.Vec3(0.3 * i, 0.0, 0.0))  # nm
        if i:
            bonds.addBond(i - 1, i, 0.3, 1000.0)  # length nm, k kJ/mol/nm^2
    system.addForce(nonbonded)
    system.addForce(bonds)
    return system, positions


def verify_platform(name: str, properties: dict[str, str] | None = None) -> tuple[bool, str]:
    """Run a tiny system on this platform. Returns (ok, detail).

    `detail` carries the device description on success (DeviceName for CUDA and
    OpenCL, thread count for CPU) and the exception text on failure. Never
    raises: every failure mode here is a reason to fall back, not to stop.
    """
    try:
        import openmm as mm
        from openmm import unit
    except ImportError as err:
        return False, f"openmm not importable: {err}"

    context = None
    integrator = None
    try:
        platform = mm.Platform.getPlatformByName(name)
        system, positions = _tiny_system(mm)
        integrator = mm.VerletIntegrator(0.001)  # ps
        context = mm.Context(system, integrator, platform, dict(properties or {}))
        context.setPositions(positions)
        # getForces exercises the force kernels and the device-to-host copy.
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            return False, f"non-finite initial energy ({energy})"
        integrator.step(10)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not math.isfinite(energy):
            return False, f"non-finite energy after 10 steps ({energy})"
        return True, _device_detail(platform, context)
    except Exception as err:  # Driver, nvrtc, OOM, missing device, plugin stub.
        return False, f"{type(err).__name__}: {err}".strip()
    finally:
        # Release the GPU context before the next candidate or the real run.
        del context
        del integrator


def _device_detail(platform, context) -> str:
    """Whatever the live Context can say about the hardware it is running on."""
    parts = []
    try:
        available = set(platform.getPropertyNames())
    except Exception:
        available = set()
    for key in ("DeviceName", "DeviceIndex", "Precision", "Threads"):
        if key not in available:
            continue
        try:
            value = platform.getPropertyValue(context, key)
        except Exception:
            continue
        if value:
            parts.append(f"{key}={value}")
    return ", ".join(parts) or "no device properties reported"


def select_platform(prefer: str | None = None, *, precision: str | None = None,
                    verify: bool = True) -> PlatformChoice:
    """Pick the fastest platform that actually runs; fall back toward CPU.

    Priority: explicit argument > SCI_PLATFORM > SCI_DEVICE > CUDA > OpenCL >
    CPU > Reference. An explicit request is honored verbatim (capitalization
    aside) and is not verified, matching `_device.py`: a forced platform must
    fail loudly rather than quietly become something else.

    Without an explicit request and with `verify=True`, each candidate is run
    before it is accepted and a failing one is skipped with its reason recorded
    in `notes`. With no OpenMM at all the choice is CPU, unverified, noted --
    the caller's own `import openmm` is where that should surface.
    """
    requested = (prefer or os.environ.get("SCI_PLATFORM")
                 or os.environ.get("SCI_DEVICE") or "").strip()
    if requested:
        name = _canonical_name(requested)
        source = "argument" if prefer else "SCI_PLATFORM/SCI_DEVICE"
        return PlatformChoice(
            name=name,
            properties=platform_properties(name, precision=precision),
            verified=False,
            notes=(f"{name} requested explicitly ({source}); honored without verification",),
        )

    available = registered_platforms()
    notes: list[str] = _load_notes()
    if not available:
        notes.append("OpenMM not importable or no platform registered; "
                     "reporting CPU without verification")
        return PlatformChoice(name="CPU", properties={}, verified=False, notes=tuple(notes))

    candidates = [n for n in PLATFORM_PREFERENCE if n in available]
    candidates += [n for n in available if n not in PLATFORM_PREFERENCE]
    for name in candidates:
        properties = platform_properties(name, precision=precision)
        if not verify:
            notes.append(f"{name} selected without verification (verify=False)")
            return PlatformChoice(name, properties, False, tuple(notes))
        ok, detail = verify_platform(name, properties)
        if ok:
            notes.append(f"{name} verified: {detail}")
            return PlatformChoice(name, properties, True, tuple(notes))
        notes.append(f"{name} rejected: {detail}")

    # Every registered platform failed its check. Reference is the last resort
    # and is still worth returning: the caller may be on a machine where only
    # the GPU plugins are broken, and a loud failure there is more useful than
    # a silent None.
    notes.append("no registered platform passed verification; falling back to Reference")
    return PlatformChoice("Reference", {}, False, tuple(notes))


def platform_report(choice: PlatformChoice | None = None) -> str:
    """One-line summary for logs and provenance sidecars.

    Record this alongside the trajectory: platform and precision change
    trajectories in detail, and checkpoints are platform specific.
    """
    choice = choice or select_platform()
    precision = choice.properties.get("Precision")
    bits = [choice.name]
    if precision:
        bits.append(f"precision={precision}")
    bits.append("verified" if choice.verified else "unverified")
    return f"{', '.join(bits)}"


def _main() -> None:
    try:
        import openmm as mm

        print(f"OpenMM version: {mm.Platform.getOpenMMVersion()}")
    except ImportError as err:
        print(f"OpenMM not installed ({err}).")
        print("Install a CUDA-enabled build; see the skill's SKILL.md.")

    if _CUDA_RUNTIME_DIRS:
        print(f"CUDA runtime DLL directories on the search path: {', '.join(_CUDA_RUNTIME_DIRS)}")
    registered = registered_platforms()
    print(f"Registered platforms (by reported speed): {', '.join(registered) or 'none'}")
    print("  Registration means the plugin loaded, not that the device works.")
    print("  Plugin load failures, if any, are listed under the selection notes.")

    choice = select_platform()
    print()
    print(f"Selected: {platform_report(choice)}")
    if choice.properties:
        print(f"Properties: {choice.properties}")
    for note in choice.notes:
        print(f"  - {note}")


def _test_installation() -> int:
    """Run openmm.testInstallation after the DLL path is prepared (module import).

    Plain `python -m openmm.testInstallation` on Windows imports openmm before
    anything can fix the DLL search path, so it hides a working CUDA platform.
    """
    import runpy

    try:
        runpy.run_module("openmm.testInstallation", run_name="__main__")
    except ImportError as err:
        print(f"OpenMM not installed ({err}).")
        return 1
    return 0


if __name__ == "__main__":
    if "--test-installation" in sys.argv[1:]:
        sys.exit(_test_installation())
    _main()
