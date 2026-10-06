# Detector-list syntax update

Only src/mollerdig.cxx changes executable behavior:

1. In the detector setup loop, add:

```cpp
if(detectors_list[k] == "showermax") {
  doShowerMax = true;
}
```

2. Replace the old enable_showermax parser assignment with a deprecation
warning. It no longer changes doShowerMax. Remove that key from your .dat.

The existing if(doShowerMax) setup, tree creation, event processing and output
writing remain unchanged. The v2 guard that skips UnfoldEvent when no GEM is
configured remains in place.

The example .dat now uses detectors_list showermax for your first SM-only run.
Use exactly one active detectors_list line. Lookup paths, optional CSV maps,
rate normalization and electronics configuration remain unchanged.

Rebuild: bash scripts/build.sh
Run: ./build/mollerdig db_gmn_conf_8gemmodules_GMN9.dat filelist_showermax.txt 100
