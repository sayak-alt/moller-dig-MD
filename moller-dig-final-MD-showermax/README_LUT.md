> Detector selection changed in v3: use `detectors_list showermax` (or add showermax alongside mollerpmt/mollergem). The old enable_showermax key is ignored. See README_FINAL.md and CHANGES_V3.md.

# LUT PMT response integration

Replace the matching files in your existing project and add include/MOLLERLookupResponse.h and src/MOLLERLookupResponse.cxx. Other uploaded sources are included unchanged. No CSV tables are bundled.

Configuration (before endconfig):

    pe_method_mollerpmt 0
    lut_directory_mollerpmt look_up_table
    lut_ring_mollerpmt 0

0 preserves the existing Frank-Tamm method (also the default if omitted). Set pe_method_mollerpmt to 1 to use lookup tables. lut_ring_mollerpmt 0 loads rings 1..6; set 1..6 for one ring. Paths are relative to the working directory of the executable. Tables must be rN_PositionDepTable.csv (13 columns in the original order) and rN_EnergyDepTable.csv (energy,scale), with headers. All requested tables are loaded and validated at setup; missing or malformed tables fail explicitly.

The new response reproduces the supplied standalone map, rotations, y offsets, electron/positron selection, laboratory pz>0 and kinetic energy k>1.1 MeV, angular clamp to 4 degrees, first matching position bin and energy window, and energy fallback scale 1. Unknown geometry IDs, unmatched positions, nonfinite or negative yields are rejected. Front-face filtering remains inactive, matching the supplied standalone macro.

Mean PE is accumulated per detector. The exact existing detector-level Poisson draw and waveform/ADC loop is shared by both methods: earliest detector hit timing, Gaussian timing jitter, exponential SPE charge, pedestal/noise, quantisation, clipping and zero-PE handling remain unchanged. In LUT mode earliest timing and hit output include accepted LUT hits only. Rate weighting is analysis-only and is not used to scale digitisation PE.

Existing PMT output branch names are retained. LUT meanpe is the detector sum of expected PE; npe_poiss is its Poisson sample. edep_mev uses sum deposits if available, otherwise positive deposits from accepted hits. leff_cm is NaN in LUT mode because no effective track length is inferred. The old vector-based DigitizeEvent interface remains unchanged; the driver selects DigitizeEventLUT when configured.

Build in your configured ROOT/Geant4/remoll environment using your normal CMake command; the new source and ROOT Physics component are included in CMakeLists.txt. For example:

    cmake -S . -B build -DREMOLL_DIR=/absolute/path/to/remoll
    cmake --build build -j4

Run your existing executable command with the updated database. Start with one ring and the same ROOT input used in the standalone validation. Compare sums of expected meanpe by event/detector before checking stochastic ADC results.

Validation here: C++17 syntax checks using lightweight ROOT API stubs; deterministic original/refactored Frank-Tamm output comparison; synthetic LUT tests for mean yield, rejected detector/position, and missing tables. These do NOT validate actual ROOT rotations, remoll field availability, linking, the entire driver, or physical response. ROOT/remoll and your real LUT CSVs were unavailable here, so a production build/run remains necessary. The geometry map is from the earlier uploaded standalone helper, not an independently verified latest geometry.

## Electronics mapping update

The bundled MOLLER_detector_electronics_map.csv contains 224 quartz IDs, using provisional ROC 1, slots 3..16 and channels 0..15. Set detector_map_mollerpmt to its absolute path or a path relative to the launch directory. An omitted key preserves operation without a map. The map does not modify PE, random sampling, or ADC signals.

Added aligned per-detector branches: rocid, slot, channel, segment, ring, segment_group, orientation, subdivision, mapping_valid, channel_enabled. segment_group is the detector-ID grouping; orientation is the actual FF/BF quartz designation from the checked online geometry. Added waveform branches: wf_rocid, wf_slot, wf_channel, wf_mapping_valid. Unknown IDs retain their signal and get mapping_valid=0, numeric map values -1 and string values unknown. The enabled column is metadata only in this version: channel_enabled=0 does NOT suppress digitisation.

The CSV parser requires the supplied header/order, accepts CRLF, and checks duplicate IDs, duplicate enabled addresses, field ranges and labels. No CODA bank encoding is added. Hardware addresses are provisional and slot 3..16 is a simulation convention, not a verified physical crate layout.

Mapping checks with lightweight ROOT stubs verified all 224 rows load, known/unknown mappings, waveform alignment, both response modes and unchanged deterministic Frank-Tamm ADC output. Full ROOT/remoll build and runtime checks must be performed in the actual environment.
