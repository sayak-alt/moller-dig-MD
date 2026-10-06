> Detector selection changed in v3: use `detectors_list showermax` (or add showermax alongside mollerpmt/mollergem). The old enable_showermax key is ignored. See README_FINAL.md and CHANGES_V3.md.

# Shower-max extension

This package extends the existing MD LUT/channel-map project. The main-detector configuration, output and algorithms remain intact. New sources: include/MOLLERShowerMaxLookupResponse.h and src/MOLLERShowerMaxLookupResponse.cxx. Updated: MOLLERPMTDigitizer.h/.cxx, mollerdig.cxx, CMakeLists.txt, database. Eight original uploaded CSV tables are bundled unchanged in showermax_tables. The provisional 28-channel map uses ROC 2, slot 3 channels 0..15 and slot 4 channels 0..11, avoiding the MD ROC 1 map.

Enable before endconfig:

    enable_showermax 1
    lut_directory_showermax /absolute/path/to/showermax_tables
    legacy_pion_names_showermax 1
    detector_map_showermax /absolute/path/to/MOLLER_showermax_electronics_map.csv

The map line is optional. enable_showermax defaults to 0. Shower-max enablement is independent of pe_method_mollerpmt (0/1 for the main detector). Both can run together. To run shower-max alone, omit mollerpmt from detectors_list while retaining enable_showermax 1. Electronics settings are now independent via *_showermax keys. Their default values and calculation match MD, but changing MD settings does not change shower-max. A separate digitiser/random generator keeps MD's random sequence unaffected. Both use Rseed; this can correlate their RNG sequences and is not intended as a realistic joint-noise correlation model.

Build:

    cmake -S . -B build -DREMOLL_DIR=/absolute/path/to/remoll
    cmake --build build -j4
    ./build/mollerdig db_gmn_conf_8gemmodules_GMN9.dat filelist.txt 1000

The current driver writes *_dig.root beside each input, as before. Output adds showermax_digi with the same PE, ADC, mapping and optional waveform branch names as pmt_digi. It has one entry per processed event. detid is the entrance-plane ID 73001..73028, identifying the predicted module signal, not an internal quartz-layer ID. Deposited energy and effective length are NaN for shower-max: entrance-plane total energy is not deposited energy. Provisional map ring is 7, segment is sector 1..28, FF/BF fields are NA.

Plot:

    showermax_digi->Draw("meanpe", "detid == 73003");
    showermax_digi->Draw("meanpe", "mapping_valid == 1 && rocid == 2 && slot == 3 && channel == 2");

Sector pattern follows the supplied document: 1 closed, 2 transition, 3 open, 4 transition, repeating. Open long-pass factor is .22, others .33.

Response reproduces the supplied nearest-module coordinate transform with radius 1100 mm, equally weighted x-linear/y-quadratic fits, linear energy interpolation and detector-ID filter factor. Selection: 73001..73028, 1020<radius<1180 mm, total energy e>10 MeV, lab pz>0 and one of eight supported PIDs. It does not sum responses from internal layers. Angular response is absent, matching the supplied helper.

Energy behaviour now matches the uploaded helper exactly: fixed grid 5,10,50,100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000 MeV. Energies outside its intervals (including exactly 9000) use fallback bounds 5 and 9000. Rows at 10000/11000 are loaded but not used. Missing endpoint rows evaluate to zero, as default TF2 coefficients do in the uploaded helper. Identical duplicate 1000 rows are merged without changing their values; conflicting duplicates fail validation. No new energy-range rejection is applied. This intentionally retains extrapolation and missing-row behaviour to reproduce the reference implementation.

IMPORTANT pion association: legacy_pion_names_showermax 1 preserves the supplied helper's PID +211 -> pi- filename and -211 -> pi+ filename. This is opposite standard PDG labels. Set 0 for PDG-consistent filenames only after confirming how these tables were generated. No silent charge reassignment is made.

Tests here used lightweight ROOT API stubs for class syntax and electronics, and the actual CSVs for direct numeric response tests. These checks cover known PE values, filter factors, particle range rejection, endpoint handling, summation, mapping and event clearing. Actual ROOT linking, remoll driver compilation, geometry/sign correctness and physical response need local build/run validation. Compare standalone and integrated expected PE using matching hit selection and per-detector event sums, not the standalone rate-weighted radial histogram.

## Separate digitiser classes
MD uses MOLLERPMTDigitizer exactly as in the pre-showermax channel-map package. Shower-max uses MOLLERShowerMaxDigitizer with its own Config, SetConfig, SetRandomSeed, SetEventNumber, BookBranches, Clear, DigitizeEvent, CSV mapping loader, RNG, waveform buffers and electronics functions. Its lookup model remains MOLLERShowerMaxLookupResponse. No shower-max response code is embedded in the MD class. The driver coordinates both classes in the same event loop and writes separate trees.

Independent shower-max electronics keys: store_waveforms_showermax, gatewidth_showermax, dt_showermax, nbits_showermax, vrange_showermax, R_showermax, ped_showermax, pedsigma_showermax, trigoffset_showermax, sigmatime_showermax, tau_showermax, Qpe_showermax. Values are included in the bundled database, matching existing MD settings. The two classes intentionally duplicate the short electronics calculation so they can be changed independently. No MD electronics behaviour was modified.
