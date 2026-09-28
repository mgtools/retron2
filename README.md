## Tools for investigating array-associated reverse transcriptases

This repository contains a candidate retron ncRNA scanner, a browser-based
alignment and structure viewer, and a small CT-to-dot-bracket converter.

## Requirements

- Python 3 for the command-line tools. They use only the Python standard library.
- A modern browser for the viewer. It is a static HTML file and does not need a
	build step or server; open `viewer/index.html` directly.

## Retron ncRNA scanner

[`scanner/retron_ncrna_finder.py`](scanner/retron_ncrna_finder.py) searches one
or more DNA sequences in a FASTA file for candidate retron ncRNA (msr-msd)
regions. It looks for an inverted repeat with configurable arm and loop lengths,
allows a configurable number of arm mismatches, and evaluates the guanines at
the inner arm junctions. It searches both strands by default and removes
substantially overlapping hits, keeping the higher-scoring hit.

Run a scan and write FASTA, TSV, and GFF results:

```bash
python3 scanner/retron_ncrna_finder.py \
	--input genome.fasta \
	--output candidates.fasta \
	--tsv candidates.tsv \
	--gff candidates.gff \
	--verbose
```

The FASTA output contains each candidate sequence and a header with its source
sequence, strand, 1-based inclusive coordinates, length, stem/loop properties,
junction bases, and score. The TSV provides the same information in columns;
the GFF output contains genomic intervals. The score is a heuristic ranking,
not a probability or a substitute for biological validation.

Useful options (run with `--help` for the complete list):

- `--min-arm`, `--max-loop`, and `--max-mismatch` control the inverted-repeat
	search. Defaults are 8 nt, 400 nt, and 1 mismatch.
- `--min-length` and `--max-length` constrain candidate length (defaults: 50
	and 500 nt).
- `--region-start` and `--region-end` limit the search to a 1-based inclusive
	interval in each input sequence.
- `--require-g` keeps only candidates with G at both inner arm-junction
	positions. Without it, candidates with missing junction Gs may still be
	reported.
- `--min-score` sets a score threshold and `--top N` limits the final output to
	the highest-scoring candidates.
- `--forward-only` searches only the input strand; otherwise both strands are
	searched.

### Example genomes

The [`scanner/`](scanner/) directory includes two complete phage genome FASTA
files: `MarsHill.fna` (Staphylococcus phage MarsHill, accession MW248466.1) and
`SA1.fna` (Staphylococcus phage vB_StaM_SA1, accession MW218148.1). The example
search windows are taken from `scanner/commands.sh`. Run the following from the
repository root to write separate results for each window:

```bash
python3 scanner/retron_ncrna_finder.py -i scanner/MarsHill.fna \
	-o MarsHill-ncrna.fna --gff MarsHill-ncrna.gff --tsv MarsHill-ncrna.tsv \
	--region-start 213847 --region-end 215077 --min-arm 8 --require-g

python3 scanner/retron_ncrna_finder.py -i scanner/SA1.fna \
	-o SA1-9447-ncrna.fna --gff SA1-9447-ncrna.gff --tsv SA1-9447-ncrna.tsv \
	--region-start 9447 --region-end 10635 --min-arm 8 --require-g

python3 scanner/retron_ncrna_finder.py -i scanner/SA1.fna \
	-o SA1-99715-ncrna.fna --gff SA1-99715-ncrna.gff --tsv SA1-99715-ncrna.tsv \
	--region-start 99715 --region-end 100172 --min-arm 8 --require-g
```

These commands scan both strands, as the scanner does by default. Remove the
`--region-start` and `--region-end` options to search each full genome instead.

## Alignment and structure viewer

[`viewer/index.html`](viewer/index.html) displays a multiple-sequence alignment
alongside the secondary structure of a reference sequence. It includes an
alignment conservation view, a structure arc plot, and a 2D structure view.
Hovering over alignment positions links the highlighted sequence columns and
structure positions; the structure view can be zoomed and panned. A theme
toggle and alignment text-size controls are also available.

Open the page in a browser, then load:

1. An aligned FASTA file containing sequences of equal alignment length. Use
	 `-` or `.` for gaps.
2. A structure file with an optional `>` header, then an ungapped reference
	 sequence, then its dot-bracket string. The sequence and structure string
	 should have equal lengths and correspond to one of the aligned sequences.

The structure file is optional in the interface, but is needed to show the
structure panels. The viewer currently recognizes standard `()` base pairs;
other bracket types used to encode pseudoknots are not interpreted as pairs.
Select **Load sample data** to try the example files in [`viewer/`](viewer/).
The hosted viewer is also available at
[omics.informatics.indiana.edu/myRT/ART/SA1](https://omics.informatics.indiana.edu/myRT/ART/SA1/).

## CT to dot-bracket converter

[`misc/CT2dot.py`](misc/CT2dot.py) converts a CT secondary-structure file into a
three-line dot-bracket file accepted by the viewer:

```bash
python3 misc/CT2dot.py input.ct output.dot
```

Omit `output.dot` to print the result to standard output. Nested pairs are
written as `()`; crossing pairs are assigned additional bracket types
(`[]`, `{}`, `<>`) to represent pseudoknots, up to four levels. The viewer
currently renders only the `()` pairs.

## Development note

Claude.ai was used during development. Prompts used for these tools were
intended to be shared here as examples.

