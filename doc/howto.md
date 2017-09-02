# How to Use the Ntuplizer

## Run the Ntuplizer

Open the `macros/runNTGJ.C` in an editor.

First ignore the first six lines of Bash preamble wrapped within an `#if 0` and `#endif`. This allows the macro to be executed like a shell script (see below), while being legitimate C++ code.

The `NULL`-terminated `const char *package[]` enumerates the packages needed for the execution. They have to be legitimate packages on the CVMFS and Grid (not from a local `aliBuild`). If this macro is executed in a local `aliBuild` environment, this versions will be replaced to make sure the initial compilation of the ntuplizer succeeds.

Next, go to `plugin->SetAliROOTVersion()` and `plugin->SetAliPhysicsVersion()` and set them to valid versions.

Several `-1`-terminated run lists like `const int run_number_lhc15o[]` are defined. They are separately named to facilitate the storage of multiple periods in this macro. Then, after an appropriate `plugin->SetGridDataDir()` and `plugin->SetDataPattern()`, a `const int *run_number` points to the correct one.

Finally, there is a section that fixes the current practice of the output storage element being unpredictable instead of honoring `alien_CLOSE_SE`.

Then execute the Grid job as:

```bash
./macros/runNTGJ.C test
```

or:

```bash
./macros/runNTGJ.C full
```

The macro will automatically execute e.g. `root -l -b -q ./macros/runNTGJ.C("full")`.

## Run the Ntuple Merger

Do not attempt to merge ntuples on the Grid. They will inevitably fail as ALICE always call `hadd` to merge, which is in-memory. Instead, another macro is provided for disk-resident merging. For example, to merge `rx/lhc16h2a_bis-246994-1/001/AnalysisResults.root`, `rx/lhc16h2a_bis-246994-1/002/AnalysisResults.root`, etc., execute:

```bash
./macros/MergeNtuple.C "rx/lhc16h2a_bis-246994-1/*/*.root" merged_output.root
```

Note to quote the wildcard, in order to prevent its expansion by the shell before calling the macro. The wildcard is then expanded (macro internally) via shell, and is not limited to ROOT&rsquo;s `TChain` expansion (which does not allow merging across multiple directories).

`MergeNtuple.C` will prefer ROOT 6, if it is available in the local `aliBuild` environment. It is also possible to compile `MergeNtuple.C` via ROOT ACLiC.
