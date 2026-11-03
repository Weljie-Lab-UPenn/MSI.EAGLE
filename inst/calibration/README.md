# Offline mass-calibration ion library

Version 2026.10.27.1 contains 19,085 positive-mode and 10,123 negative-mode ion entries, including lipid classes, small metabolites and Leu-Enk lock references. `compounds.tsv` links compounds/formulas with the class-appropriate adducts available in each polarity. This is broad biochemical coverage, not an exhaustive chemical database or a list of confirmed identities in a sample.

## GUI workflow

Load a saved dataset from the working-directory dropdown, choose polarity and select adducts appropriate to the acquisition solvent. Common adducts are enabled by default; formate and neutral-loss options can be enabled explicitly. Click Find anchors. The default table shows eligible matched candidates and selects a recommendation balanced across eight mass intervals, with no more than two recommendations per interval and one per reference family. Click rows to change the selection. Show weaker and ambiguous matches exposes other observed matches and their rejection reasons. Compare fits, select a model row to update its curve, then Apply and Save.

The anchor table shows mean_intensity (mean across loaded pixels, in the dataset's existing intensity units). Dominance is the ratio to the second strongest observed peak in the search window. No competing peak explicitly displays an unbounded ratio; finite ratios remain numerically sortable.

Candidate selection does not require an evidence-entry step. Selecting a rejected candidate does not bypass its peak-quality or ambiguity checks. Changing polarity, adducts, references or search mode requires Find anchors again. Changing the displayed candidate view resets row selection to the recommendations in that view.

## Library construction and provenance

`scripts/build_calibration_library.py` accepts broad positive/negative formula-bearing source tables and their database-verification tables. It writes the two ion libraries, a common compound index and a version/hash manifest. The initial seed comes from the broad computed manuscript libraries, rather than its small anchor tables. All source files were read only. No manuscript folder is needed at runtime.

Every generated ion mass is calculated from its ion elemental formula, explicit charge and electron mass. Generated masses were checked against the source tables within 0.000002 Da. The instrument-specific Leu-Enk references are retained explicitly as overrides. `source_scope` distinguishes compounds with database records from computational species supported by database examples of the class rule. `source_ids` retain the corresponding LIPID MAPS, PubChem, KEGG or ChEBI identifiers as present in the source verification snapshots. These are provenance, not claims of newly performed live database verification or confirmation of every enumerated species.

Rebuild offline with:

    python3 scripts/build_calibration_library.py POSITIVE_SOURCE.tsv NEGATIVE_SOURCE.tsv POSITIVE_CHECKS.tsv NEGATIVE_CHECKS.tsv inst/calibration

The source tables require name, class, adduct, formula, mz and prior. Verification tables require name, cls, check, db1 and db2. Source hashes and final asset hashes are recorded in manifest.json. Updating sources is an explicit maintenance operation; Find anchors does not contact a database.

## Candidate discovery

Enabled ions are grouped by elemental ion formula and charge, preserving alternative names and adduct pathways. They are one calibration target even when their chemical identities cannot be distinguished. Different ion compositions within 5 ppm are counted automatically against the full enabled library; historical alternative counts in custom files are ignored. All alternatives are retained in the report. Candidate matching uses ±20 ppm, or ±60 ppm below m/z 300 in extended mode. Peak dominance, abundance, available isotope agreement and M+2 interference are checked. Only the nearest eligible assignment to a shared observed feature can be used. Recommendations prioritize eligible, abundant candidates and coverage; they do not constitute independent mass identifications.

Custom TSVs require name, adduct, formula (base molecular or permanent-ion formula) and theo. Ion formulas are computed for supported singly charged adducts when absent. Optional columns include ion_formula, charge, reference_family, prior, source_ids and library_version. Same-composition targets with conflicting theoretical masses are rejected. Unsupported charge states and inconsistent polarity are rejected.

## Calibration and persistence

One shared processed feature axis is corrected per dataset. Intensities and pixel metadata are preserved. Outside the anchor range, the endpoint correction is held constant. LOO is internal model-selection validation after anchor selection, not an independent accuracy estimate. Existing identifications require mass review. Save uses the same file-save dialog as the other tabs. The chosen name creates a new bundle folder containing the imzML/ibd dataset, Cardinal feature/pixel sidecars, a calibration-details text report and all TSV reports. No RDS files are generated. Original feature masses support undo.


Run `Rscript tests/calibration-checks.R` and `Rscript tests/calibration-library-checks.R` from the package root. Large-data performance, complete moved-bundle portability and live interactive use still require integration verification for production release.

## Cross-tab use

To analyze a calibration result in Statistics, save the calibrated imzML bundle, then select its imzML file in Statistics and load it. Statistics uses saved files only; there is no in-memory calibration shortcut. Open result in Peak Picking remains available in Mass Calibration.

Histology Overlay now offers Saved imzML bundle alongside current Overview data and the legacy RDS upload. Select a bundle under the working directory and click Load bundle; the loaded location is shown. Switching datasets clears mapped results, fit caches and ion selections. This does not replace the Overview dataset. Registration parameter exports remember the direct bundle path and mass/calibration signature; unavailable or changed bundles require reselection.

Calibrated imzML saves across Peak Picking, UMAP, segmentation, phenotyping and Histology persist mean/count/frequency for the actual saved pixels. Phenotyping retains the full feature table. Incompatible calibrated feature lists cannot be appended as though their identities agree. Dataset pickers exclude historical calibration-detail RDS files. New calibration exports remain imzML-only.

Regression checks: `Rscript tests/calibration-save-checks.R`, `Rscript tests/calibration-compatibility-checks.R`, and `Rscript tests/calibration-loader-checks.R`. These exercise module server logic and file-backed data; they do not replace live visual inspection of every tab.

## Feature-column policy

Calibration adds only `mz_original` and `correction_ppm` to the loaded/saved feature table. Existing source columns and standard mean/count/frequency summaries are preserved. Observed masses, outside-range flags, model/library details and report feature indices are written to the extra calibration reports instead. Earlier calibration-specific columns are removed in memory when older calibrated datasets are opened through supported dataset loaders; existing files are not rewritten. Source-provided feature IDs are preserved, while synthetic feature IDs from earlier calibration exports are removed. Calibration guards and append checks use the retained original masses and corrections, rather than requiring a calibration-ID column.
