This repository contains small utilities and configuration files used in our manuscript:

> Kawaguchi Y. W. et al. Improved genome assembly of whale shark, the world’s biggest fish: revealing “chromocline” in intragenomic heterogeneity (in prep.)

If you use these scripts, please cite the paper above and the original tools referenced below.

# Contents

## extrinsic.M.RM.E.P.cfg
Customized AUGUSTUS extrinsic configuration that integrates three hint sources:
* Repeat (e.g., RepeatMasker/TE hints)
* Protein (alignments to proteins)
* cDNA (transcript/cDNA alignments)

### Reference
> Keller, O., Kollmar, M., Stanke, M., & Waack, S. (2011). A novel hybrid gene prediction method employing protein multiple sequence alignments. Bioinformatics, 27(6), 757–763. https://doi.org/10.1093/bioinformatics/btr010

## pafCoordsDotPlotly_coloredeachaliginment_pdf.R
Modified from: pafCoordsDotPlotly.R in [tpoorten/dotPlotly] (MIT License)

Nature of change: lightweight derivative, keeping CLI behavior as close to upstream as possible.

### Added features
* PDF export support
* Per-alignment coloring (each alignment drawn in a distinct color)
* Full-chromosome axes (include unaligned regions in the plotted range)

Usage: Same as the original DotPlotly script. See upstream docs: 
> https://github.com/tpoorten/dotPlotly


## Attribution & License

- Repository license: MIT (see `LICENSE`).
- AUGUSTUS config (`extrinsic.M.RM.E.P.cfg`): Artistic-2.0 (see `LICENSES/Artistic-2.0.txt` and `THIRD_PARTY_NOTICES.md`).
- DotPlotly derivative (`pafCoordsDotPlotly_coloredeachaliginment_pdf.R`): MIT.
