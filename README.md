This repository contains small utilities and configuration files used in our manuscript:

> Kawaguchi Y. W. et al. Improved genome assembly of whale shark, the world's biggest fish: revealing intragenomic heterogeneity in molecular evolution (in prep.)

If you use these scripts, please cite the paper above and the original tools referenced below.

# Contents

## extrinsic.M.RM.E.P.cfg
A customized AUGUSTUS configuration file based on the sample extrinsic.cfg.
It integrates three hint sources: repeat, protein, and cDNA.

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

## HiC_plot.py
A small Python script to generate Hi‑C contact-map visualizations aligned to scaffold/orderings from assembly files. The script reads a three-column Hi‑C matrix, an original assembly file and a curated scaffold/order file, then produces contact-map images showing the raw assembly order and the curated scaffold order for visual comparison.

### Usage
Create the Hi-C matrix TSV first (example using Juicer Tools):
```
java -jar juicer_tools.jar dump observed KR <hic_file> assembly assembly BP 2500000 > <hic_matrix.tsv>
```
Then run the plotting script:
```
python HiC_plot.py -a <assembly_original> -c <assembly_curated> -m <hic_matrix.tsv> -o <output.png>
```
Note: the example uses a resolution of 2,500,000 bp. To change the resolution, adjust the BP value in the Juicer Tools command above and the BIN value at L1 of` HiC_plot.py`

#### Options:
* -a, --assembly0 : Original assembly file (required).
* -c, --assembly-curated : Curated assembly (scaffold order) file (required).
* -m, --hic-matrix : Hi‑C matrix TSV file (required).
* -o, --output : Output filename with extension (png or pdf). The base name is used and the script writes files named <base>.<ext>.
* -s, --scale-factor : Float multiplier applied to matrix coordinates before binning (default: 2).
* --max-scaffolds : Limit number of curated scaffolds used (set 0 for all).

### Requirements
* Python:: 3.+ 
* Libraries:: pandas, numpy, matplotlib
* JuicerTools:: v1.9.9
  * required. Outputs of newer versions are not compatible; use v1.9.9 when generating the initial .hic and assembly files with juicer_tools.jar.

## Attribution & License

- Repository license: MIT (see `LICENSE`).
- AUGUSTUS config (`extrinsic.M.RM.E.P.cfg`): Artistic-2.0 (see `LICENSE.Artistic-2.0` and `THIRD_PARTY_NOTICES.md`).
- DotPlotly derivative (`pafCoordsDotPlotly_coloredeachaliginment_pdf.R`): MIT.
- HiC contact map (`HiC_plot.py`): MIT.
