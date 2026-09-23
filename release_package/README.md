# WINERED DIB spectrum data package

This package contains sanitized copies of processed WINERED spectra selected
from the DIBproject database and release inventory.

## Download

The package can be downloaded from the following tagged GitHub Release:

- [WINERED DIB spectrum data package](https://github.com/SatoshiHamano/DIBanalyzer/releases/tag/winered-dib-spectra-candidate-2026-09-12)

## Contents

- `spectra/`: FITS files grouped by object name and `combineID`.
- `LICENSE.txt`: CC BY 4.0 license notice for the data package.
- `MANIFEST_sanitized.csv`: one row per FITS file, including object metadata,
  `combineID`, echelle order, source filename, package-relative release path,
  wavelength range, and publication tags.
- `OBJECTS.csv`: one row per object, with the database object ID, preferred
  name, aliases, RA and Dec in sexagesimal and decimal-degree forms, spectral
  type, E(B-V), included `combineID` values, and publication tags.

## Scope

- Number of objects: 50
- Number of FITS spectra: 1000
- Instrument mode: WIDE
- Processing level: processed, telluric-corrected, wavelength-corrected,
  heliocentric spectra as referenced by the DIBproject database inventory.

### Relationship to published studies

This package is based on a processed-spectrum collection that was previously
assembled for sharing with a research collaborator. That collection combined
the main 2022 study sample with additional DIB targets:

- The package contains spectra for all **31 reddened target objects** analyzed
  in Hamano et al. (2022). The unreddened reference star Rigel is not included.
- **19 objects** are supplementary DIB targets included in that collaborator
  dataset in addition to the 2022 sample.
- Of those 19 additional objects, **five** are targets included in Hamano et al.
  (2015) or Hamano et al. (2016): Cyg OB2 No. 11, HD 12953, HD 21389,
  zeta Ori A, and HD 50064.
- The remaining **14 objects** are unpublished DIB targets.

For each object, `OBJECTS.csv` records the publication tags and DOI mappings.
`MANIFEST_sanitized.csv` repeats these tags per spectrum.

Relevant papers currently identified are:

- Hamano et al. (2022), *Survey of near-infrared diffuse interstellar bands in
  Y and J bands. I. Newly identified bands*,
  DOI: https://doi.org/10.3847/1538-4365/ac7567
- Hamano et al. (2016), *Near Infrared Diffuse Interstellar Bands Toward the
  Cygnus OB2 Association*,
  DOI: https://doi.org/10.3847/0004-637X/821/1/42
- Hamano et al. (2015), *Near-infrared diffuse interstellar bands in
  0.91-1.32 micrometers*,
  DOI: https://doi.org/10.1088/0004-637X/800/2/137

## License and citation

This data package is licensed under the **Creative Commons Attribution 4.0
International License (CC BY 4.0)**. The data may be used, modified, and
redistributed under the terms of that license. See `LICENSE.txt` and
https://creativecommons.org/licenses/by/4.0/.

When using the spectra in scientific work, cite the relevant publication or
publications listed above and identified for each object in `OBJECTS.csv`. A
separate citation to the GitHub repository or GitHub Release is not required.
Prior permission, collaboration, and co-authorship are not conditions of use.

## Known validation limitations

These files are historical processed products selected through the current
DIBproject database. The following technical validation work remains:

- reconstruct and verify the observing-run, instrument-mode, and pipeline-version
  mapping for wavelength-correction parameter sets;
- validate wavelength WCS consistency in all released FITS files;
- document the legacy IRAF telluric-correction path and residual limitations.

One historical product inspected during software auditing contained conflicting
`CDELT1` and `CD1_1` values, and the applied run-specific wavelength-calibration
file is not recorded in historical output metadata. Until the package-wide
checks are complete, this data package should not be treated as validated
for precision radial-velocity or line-center measurements. This caveat does not
by itself establish that all or most spectra are affected.

## FITS header sanitization

The original FITS files were not modified. Public-sharing copies were created
with the following internal/private header keywords removed when present:

- `FITSFILE`
- `OBSERVER`
- `WODBPI`
- `WODBOBS`
- `WODBPROP`

The sanitizer also adds a `RELEASE` keyword and HISTORY entries noting that the
copy was sanitized.

## Directory layout

```text
spectra/
  OBJECT_NAME/
    COMBINE_ID/
      *_mORDER_*.fits
```

Use `MANIFEST_sanitized.csv` as the authoritative index for object IDs,
aliases, `combineID`, echelle order, wavelength range, and provenance tags.
Paths in this public manifest are relative to the package; local source paths
are deliberately omitted.

Use `OBJECTS.csv` for object identification. Its coordinates come from the
DIBproject database columns `object.ra` and `object.decli`. The coordinate
frame and equinox were not recorded in the release inventory and should be
verified before the final public release.
