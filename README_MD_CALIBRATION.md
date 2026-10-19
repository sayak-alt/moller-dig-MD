# MD channel pedestal and gain calibration

This package extends the uploaded framework; only MD calibration was changed.

## Run

Build once in your ROOT/Geant4/remoll environment:

    bash scripts/build.sh

The supplied .dat selects MD only and contains:

    detectors_list mollerpmt
    store_waveforms_mollerpmt 1
    detector_map_mollerpmt MOLLER_detector_electronics_map.csv
    calibration_file_mollerpmt MD_channel_calibration.csv

Use absolute paths if running outside this directory. Then:

    ./build/mollerdig db_gmn_conf_8gemmodules_GMN9.dat filelist_MD.txt 1000

Both PE methods use these calibrations. Existing output naming/overwrite rules
remain unchanged. Changing CSV constants requires another run, not a rebuild.
Omit calibration_file_mollerpmt to retain global pedestal/Qpe settings.

## CSV

224 unique channels, sorted by numeric segment, segment group, ring, subdivision,
and detector ID. Segment numbers follow the existing map (1 through 14); FF/BF
segment groups distinguish the associated assemblies. Geometry metadata is copied
unchanged from that map. Detector ID is the calibration key; electronics addresses
continue to come from the separate map. When a map is supplied, all IDs and
geometry metadata must match it. Missing/duplicate rows, bad numeric constants,
map mismatch, or a count other than 224 cause initialization errors.

All initial rows use channel-specific pedestal means uniformly drawn from 190–210 mV (200 mV +/-10 mV), sigma 5 ADC counts/sample, and gain 1e6. The conversion uses the supplied 1 V / 12-bit ADC: ADC counts = mV*4096/1000. Random seed: 2242026. These draws are fixed in the CSV, not redrawn each event.
These are placeholders reproducing nominal settings, not measured calibrations.
Edit rows for individual tubes, or all rows of a ring for a ring-level model.
If ADC range or bit depth changes, regenerate the voltage-derived pedestal constants.

With calibration enabled only the calibrated IDs are digitised: this excludes
non-MD deposits passed to the legacy Frank-Tamm path. Without a map, the 224
calibration IDs define membership; their physical identity cannot be cross-checked.
CSV enabled flags remain metadata; they do not disable digitisation.

## Pedestal width definition

pedestal_mean_adc is the mean baseline of ONE ADC SAMPLE, in ADC counts.
pedestal_sigma_adc is the standard deviation (RMS width) of ONE ADC SAMPLE,
in ADC counts, before rounding/clipping. Every generated sample receives an
independent Gaussian draw. There is no common-mode or time-correlated noise.

For N independent samples, the noise-only summed integral has mean N*mu and
sigma sqrt(N)*sigma_sample. For N=100, mu=300, sigma_sample=5: mean 30000,
sigma_integral about 50 counts. If the measured width is for an integrated
pedestal histogram, divide it by sqrt(N) before filling this sample-sigma column,
only if independent sample noise is a valid model. Gaussian FWHM=2.35482*sigma.
Measured widths after event-by-event baseline estimation have different statistics.

Zero-PE channels already present in the detector response now receive pedestal/noise
samples and normal integrated ADC values. The supplied configuration now enables all-224-channel readout; see below.

## Gain

Qpe=e*gain, with e=1.602176634e-19 C. This replaces the global Qpe_mollerpmt only
when a calibration file is active. PE means and Poisson draws are unchanged.
At gain 1e6 Qpe=1.602176634e-13 C, 0.136% above the old rounded 1.6e-13 C.
No single-PE multiplication fluctuations are added. To match the old rounded
charge exactly, set gain=998641.451913722.

## Output compatibility

New vectors aligned with detid: pedestal_mean_adc, pedestal_sigma_adc, pmt_gain,
qpe_c (coulombs). adc_int_pedsub changes from vector<int> to vector<double> so
fractional pedestal means are preserved. Update any typed reader accordingly.
Raw samples remain vector<unsigned short>; adc_int remains vector<int>.
Use pedestal_mean_adc[channel index], rather than a hardcoded 300, in replay.
Pulse shape, timing and gate construction remain unchanged.

## Modified files

include/MOLLERPMTDigitizer.h
src/MOLLERPMTDigitizer.cxx
src/mollerdig.cxx
 db_gmn_conf_8gemmodules_GMN9.dat
New file: MD_channel_calibration.csv

## Verification

C++17 syntax checks and deterministic ROOT-stub tests passed for all-224 loading,
per-channel overrides, gain scaling, waveform/raw/pedsub sum identities, non-MD
filtering, missing-file error, and buffer clearing. CSV uniqueness/sort checks
passed. No full ROOT/Geant4/remoll build or Gaussian-width measurement was possible
here; run those checks on ifarm. This is not a completed Podd replay implementation.

## Random pedestal variant

The CSV means are approximately 778–860 counts; 200 mV corresponds to 819.2 counts.
The signs refer to offsets above/below positive 200 mV, not negative absolute
baselines. Negative absolute voltages are not modeled by this unipolar ADC.
Sigma remains 5 counts (1.2207 mV), independent of the spread of channel means.
Global ped_mollerpmt 300 is only a fallback when the calibration file is omitted.

## Final zero-PE pedestal update

No peak finder, NSB/NSA settings, peak branches, or peak-window integration is
included. Those belong to replay. The full-gate sums remain digitiser diagnostics.

Stored detector responses with npe_poiss==0 now generate all samples of pedestal
and noise, with the channel-specific mean and sigma, and normal raw and
pedestal-subtracted integrals. With waveform storage off, sums are still computed.
All input events remain represented in pmt_digi, including empty events.
All channels are populated when readout_all_channels_mollerpmt is enabled. RNG sequencing changes because zero-PE
channels now consume Gaussian draws; subsequent draws may differ from old files.

Rebuild once, rerun digitisation, then reconstruct pulses and charge in replay.
Validation: syntax check and deterministic tests of zero-PE sample counts,
channel-calibrated pedestal sums, and positive-PE bookkeeping passed with ROOT
stubs. Full ROOT/Geant4/remoll and statistical noise checks require ifarm.

## Full 224-channel acquisition (latest)

The supplied configuration enables readout_all_channels_mollerpmt 1.
A complete 224-channel calibration CSV is mandatory in this mode. Every input
event now contains exactly 224 detector entries, ordered by ascending detector ID.
CSV rows remain sorted by segment. Both PE methods use this mode. Missing signal
yields zero expected/sampled PE; each channel still receives independent pedestal
noise and ADC quantisation. Every waveform uses the existing fixed gate.
At 100 samples/channel, wf_adc and each waveform companion vector contain
22400 elements per event when store_waveforms_mollerpmt 1 is set. Without waveform
storage, all-channel integration still runs but raw samples are not written.
Full-gate raw and pedestal-subtracted integrals include pedestal-only channels.
No fake hit times are added: hit/t0 branches retain their original truth semantics
and need not have 224 entries. In LUT mode leff_cm remains NaN. edep_mev for
channels with no accepted LUT response is zero in this implementation, even if
other nonaccepted particles deposited energy: it is not a complete truth ledger.
CSV enabled flags remain metadata and do not suppress readout.
Set readout_all_channels_mollerpmt 0 to return to sparse detector-response readout.
Pedestal and gain calibration, pulse timing, and no-peak-finding policy are retained.
All-channel loops consume more RNG draws; do not expect eventwise PE identity
against older output even with a fixed seed. Output size and runtime increase.
Unweighted PE plots now include all zero-signal channels/events; select meanpe>0
or npe_poiss>0 explicitly when a signal-only distribution is intended.

Validation: deterministic ROOT-stub tests of empty and nonempty events in both
PE modes, 224-channel counts, 22400 sample counts, signal retention, and waveform
sum identities passed. Full environment build/statistical verification remains
on ifarm.
