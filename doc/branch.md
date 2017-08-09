# Description of Ntuple Branches

## Data Types

The following C style notation are used below to describe the less obvious data types used in the ntuple.

* C strings: `const char *`
* C string array: `const char (*)[]`
* IEEE 754-2008 binary16 (&ldquo;half precision&rdquo;) stored in binary32: `__fp16`

## Run/Event Metadata

| Name                 | Type               | Description                               |
| -------------------- | ------------------ | ----------------------------------------- |
| `id_git`             | `const char *`     | Git blob name                             |
| `version_aliroot`    | `const char *`     | AliRoot version                           |
| `version_aliphysics` | `const char *`     | AliPhysics version                        |
| `version_jec`        | `const char *`     | Jet energy correction version             |
| `grid_data_dir`      | `const char *`     | AliEn data production directory           |
| `grid_data_pattern`  | `const char *`     | AliEn data pattern for the “find” command |
| `beam_particle`      | `int[2]`           | 1000&nbsp;_A_&nbsp;+&nbsp;_Z_ for each beam particle |
| `ntrigger_class`     | `unsigned long`    | Number of trigger classes                 |
| `trigger_class`      | `const char (*)[ntrigger_class]` | Trigger class names         |
| `run_number`         | `int`              | ALICE run number                          |
| `trigger_mask`       | `unsigned long[2]` | Trigger mask, each storing 50 bits        |
| `mixed_event`        | `char`             | True if ntuple underwent event mixing     |
| `multiplicity_v0`    | `float[64]`        | V0 multiplicity, per each of the 64 channels |
| `centrality_v0m`     | `float`            | V0M centrality, alias for `centrality[0]` |
| `centrality`         | `float[9]`         | Centrality for: {V0M, CL0, CL1, V0Mplus05, V0Mplus10, V0Mminus05, V0Mminus10, SPDClustersCorr, SPDTracklets} |
| `event_plane_psi_v0` | `float[3]`         | V0 event plane angle _&Psi;_ for _v_<sub>1</sub>, &hellip; _v_<sub>3</sub> |
| `event_plane_q_v0`   | `double[3][2]`     | V0 event plane _Q_ vector for _v_<sub>1</sub>, &hellip; _v_<sub>3</sub> |
| `has_misalignment_matrix` | `bool` | True if the EMCAL misalignment matrix was loaded |
| `eg_ntrial`          | `int`              | Number of trial until the current event is generated |
| `eg_perp_hat`        | `float`            | pT hat of the current event               |
| `eg_cross_section`   | `float`            | Cross section for the current event (mb)  |
| `primary_vertex`     | `double[3]`        | Primary vertex                            |

## Muon Tracks

| Name                                | Type                  | Description                                                         |
| ----------------------------------- | --------------------- | ------------------------------------------------------------------- |
| `nmuon_track`                       | `unsigned long`       | Number of muon tracks                                               |
| `muon_track_pt`                     | `__fp16[nmuon_track]` | Muon track transverse momentum _p_<sub>T</sub>                      |
| `muon_track_eta`                    | `__fp16[nmuon_track]` | Muon track pseudorapidity _&eta;_                                   |
| `muon_track_phi`                    | `__fp16[nmuon_track]` | Muon track azimuth _&phi;_                                          |
| `muon_track_p_dca`                  | `__fp16[nmuon_track]` | Muon track p&nbsp;&times;&nbsp;DCA (see below)                      |
| `muon_track_sigma_p_dca`            | `__fp16[nmuon_track]` | Uncorrected muon track &sigma;(p&nbsp;&times;&nbsp;DCA) (see below) |
| `muon_track_sigma_slope_p`          | `__fp16[nmuon_track]` | (see below)                                                         |
| `muon_track_distance_sigma_slope_p` | `__fp16[nmuon_track]` | (see below)                                                         |

To make an arbitrary, resolution corrected &ldquo;n&sigma;&rdquo; cut on p&nbsp;&times;&nbsp;DCA, the boolean `pass_cut` is calculated as:

```C++
const float nrp = nsigma_p_dca × muon_track_delta_sagitta_p[i]
const float p_resolution_effect = muon_track_sigma_p_dca[i] / (1 - nrp / (1 + nrp))
const float slope_resolution_effect = muon_track_distance_sigma_slope_p[i]
const float sigma_p_dca_resolution_effect = sqrt(std::pow(p_resolution_effect, 2) + std::pow(slope_resolution_effect, 2))
const bool pass_cut = muon_track_p_dca[i] < sigma_p_dca_resolution_effect
```
