# Test structure fixture

`GID1A_AFDB.pdb` is the AlphaFold DB predicted model for GID1A (UniProt **Q9MAA7**),
downloaded from https://alphafold.ebi.ac.uk/entry/Q9MAA7.

It is bundled only so the `test` profile can exercise the interactive variant effect
inspection tool end-to-end without structure prediction.

**Licence:** AlphaFold DB structures are distributed under **CC-BY-4.0**
(© EMBL-EBI & DeepMind). If reused, attribute AlphaFold DB / UniProt Q9MAA7 accordingly.

> TODO (release): move this fixture to the `deepmutscan` branch of `nf-core/test-datasets`
> (e.g. `deepmutscan/testdata/GID1A_AFDB.pdb`) with the same attribution, then reference it in
> `conf/test.config` via `params.pipelines_testdata_base_path` and remove this committed copy.
