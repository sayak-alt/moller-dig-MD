# Final MD and shower-max digitisation package

MD and shower-max have separate classes. MD retains selectable Frank–Tamm or
thin-quartz LUT response. Shower-max retains the supplied particle/table
association, coordinate transform, fixed interpolation grid, missing-row
behavior, light-pipe factors and cuts. Both use the existing detector-level
Poisson fluctuations and electronics model.

## Build

Load your ROOT and Geant4 environments, then:

```bash
export REMOLL_DIR=/absolute/path/to/remoll
bash scripts/build.sh
```

REMOLL_DIR must contain include/ and lib64/libremoll.so, as required by this
project's CMakeLists.txt. Run from a normal shell; reroot is not required.

## Configuration

Edit db_gmn_conf_8gemmodules_GMN9.dat. Put these settings before endconfig:

```text
pe_method_mollerpmt 1
lut_ring_mollerpmt 0
lut_directory_mollerpmt /absolute/path/to/thin_quartz_tables
detectors_list mollerpmt showermax
lut_directory_showermax /absolute/path/to/this/package/showermax_tables
legacy_pion_names_showermax 1
rate_normalization_divisor 1
```

Use pe_method_mollerpmt 0 for Frank–Tamm. Thin-quartz tables are not bundled;
provide rN_PositionDepTable.csv and rN_EnergyDepTable.csv for each selected ring.
Shower-max tables are bundled unchanged. Keep mollerpmt in detectors_list to
run MD. Shower-max is enabled by the showermax token in detectors_list.

Optional maps:

```text
detector_map_mollerpmt /absolute/path/to/this/package/MOLLER_detector_electronics_map.csv
detector_map_showermax /absolute/path/to/this/package/MOLLER_showermax_electronics_map.csv
```

Omit mapping keys to use detector IDs alone. Example ROC/slot/channel assignments
are provisional wiring assignments, not a verified installed hardware map.

rate_normalization_divisor corresponds to goodFileCount in the standalone
analysis. Use the same divisor as that analysis: 1 for a single file; for an
ensemble normalized by its number of good files, set that number explicitly.
It is not the number of events processed. Partial event processing produces
partial weighted totals; the code does not compensate for missing events.

## Run

filelist.txt contains one input remoll ROOT path per line:

```text
/absolute/path/to/input1.root
/absolute/path/to/input2.root
```

```bash
bash scripts/run.sh db_gmn_conf_8gemmodules_GMN9.dat filelist.txt 1000
# Omit 1000 to process all entries.
```

The existing driver writes input1_dig.root beside input1.root using RECREATE;
rerunning replaces that output. MD is in pmt_digi; shower-max in showermax_digi.
The existing driver event-selection behavior is preserved.

## Output and plots

Both trees retain meanpe (expected PE), npe_poiss (sampled PE), ADC and mapping
branches. New branches are event_rate_hz (raw input rate), rate_GHz (normalized
rate), rate_x_meanpe and rate_x_npe (aligned with detid, units GHz times PE).
Rates never enter the Poisson draw or electronics calculations. Missing input
rate branches are reported as errors.

```bash
root -l 'scripts/plot_pe.C("/path/input1_dig.root","showermax_digi",73003)'
root -l 'scripts/plot_pe.C("/path/input1_dig.root","showermax_digi",0)'
root -l 'scripts/plot_pe.C("/path/input1_dig.root","pmt_digi",0,5)'
# Poisson-sampled PE:
root -l 'scripts/plot_pe.C("/path/input1_dig.root","showermax_digi",73003,0,true)'
```

The three panels show unweighted PE, PE weighted by rate, and the distribution
of rate times PE. Each entry is a detector total within one event. This is not
the standalone per-hit histogram: detector aggregation changes the PE
spectrum, while the summed rate times expected PE agrees for the same accepted
hits/events and normalization. Modules without accepted hits have no entry.
Ring filtering uses map metadata and therefore requires a valid map.

## Validation limits

The lookup equations were checked against direct reference calculations in
64 particle/energy cases. Class syntax and bookkeeping were checked with ROOT
stubs. ROOT/Geant4/remoll are unavailable in the development environment, so a
full real-dependency build and event comparison still need your installation.
The reference macro uses Float_t intermediate values and stored hit.r; this
implementation uses doubles and hypot(x,y), which may differ at cut boundaries.

## Shower-max only

Use exactly one active detectors_list line.
Set detectors_list mollerpmt showermax and the shower-max lookup directory. With no active
detectors_list line, MD and GEM remain disabled. The GEM handler is skipped
when no GEM detector is configured. The normal run command still applies.

## Updated detector-selection syntax (v3)

The supplied configuration now defaults to shower-max only:

```text
detectors_list showermax
```

For both PMT systems, use `detectors_list mollerpmt showermax`.
For all three implemented systems, use
`detectors_list mollergem mollerpmt showermax`.
The old enable_showermax key is deprecated and ignored, with a warning.
Remove it from existing configurations. Detector selection comes entirely
from detectors_list. Rebuild once after replacing src/mollerdig.cxx.
No detector class, PE calculation, or electronics code changed in v3.
