# Description of Ntuple Branches

## Type description

The following C style naming are used below for data types

* C strings: const char \*
* C string array: const char (\*)[]
* IEEE 754-2008 binary16 (“half precision”) stored in binary32: \_\_fp16
* IEEE 754-2008 binary32: float

## Run/Event Metadata

| Name                  | Type             | Description                               |
| --------------------- | ---------------- | ----------------------------------------- |
| id\_git               | const char *     | Git blob name                             |
| version\_aliroot      | const char *     | AliRoot version                           |
| version\_aliphysics   | const char *     | AliPhysics version                        |
| version\_jec          | const char *     | Jet energy correction version             |
| grid\_data\_dir       | const char *     | AliEn data production directory           |
| grid\_data\_pattern   | const char *     | AliEn data pattern for the “find” command |
| beam\_particle        | int[2]           | A×1000+Z for each beam particle           |
| ntrigger\_class       | unsigned long    | Number of trigger classes                 |
| trigger\_class        | const char (*)[ntrigger\_class] | Trigger class names        |
| run\_number           | int              | ALICE run number                          |
| trigger\_mask         | unsigned long[2] | Trigger mask, each storing 50 bits        |
| mixed\_event          | char             | True if ntuple underwent event mixing     |
| multiplicity\_v0      | float[64]        | V0 multiplicity, divided by channel       |
| centrality\_v0m       | float            | V0M centrality, alias for centrality[0]   |
| centrality            | float[9]         | Centrality for: {V0M, CL0, CL1, V0Mplus05, V0Mplus10, V0Mminus05, V0Mminus10, SPDClustersCorr, SPDTracklets} |
| event\_plane\_psi\_v0 | float[3]         | V0 event plane angle Ψ for v1, … v3       |
| event\_plane\_q\_v0   | double[3][2]     | V0 event plane Q vector for v1, … v3      |
| has\_misalignment\_matrix | bool | True if the EMCAL misalignment matrix was loaded  |
| eg\_ntrial            | int              | Number of trial until the current event is generated |
| eg\_perp\_hat         | float            | pT hat of the current event               |
| eg\_cross\_section    | float            | Cross section for the current event (mb)  |
| primary\_vertex       | double[3]        | Primary vertex                            |
