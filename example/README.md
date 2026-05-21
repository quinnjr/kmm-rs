# kmm-rs example

Self-contained, fully-deterministic fixture for the KMM plugin. Two synthetic
"classes" with biased nucleotide composition demonstrate that the
*k*-Markov-model scorer correctly classifies reads back to the distribution
they came from.

| Class | Nucleotide bias | Training length | Test reads |
|---|---|---:|---:|
| `GC_rich` | 40 % G, 40 % C, 10 % A, 10 % T | 50 000 bp | 10 |
| `AT_rich` | 40 % A, 40 % T, 10 % G, 10 % C | 50 000 bp | 10 |

Order: 3 (4-mer transition probabilities).

## Files

```
example/
├── parameters/kmm.txt        PluMA parameter file (paths relative to example/)
├── models/
│   ├── GC_rich.model         Trained k-Markov model — basename = class label
│   └── AT_rich.model
├── training/
│   ├── GC_rich.fa            Committed training "genomes" (auditable)
│   └── AT_rich.fa
├── sequences.fasta           20 unlabeled reads to classify
├── expected_output.tsv       Reference output for regression testing
├── config.txt                Minimal PluMA pipeline (Prefix + one Plugin line)
├── run_pluma_example.sh      One-shot wrapper: builds .so, runs pluma, diffs output
└── README.md                 (this file)
```

## Running through the pluma binary

```bash
./example/run_pluma_example.sh                    # uses `pluma` from $PATH
./example/run_pluma_example.sh /path/to/pluma     # explicit binary
PLUMA=/path/to/pluma ./example/run_pluma_example.sh
```

The script builds the release `.so`, stages it under a temporary
`PLUMA_PLUGIN_PATH` as `KMM/libKMMPlugin.so`, runs `pluma example/config.txt`
from the repo root, and diffs `example/out/results.tsv` against the committed
`expected_output.tsv`. Returns non-zero on mismatch.

The plugin's `PluMAPlugin::input()` impl resolves the `models` and `input`
paths inside `parameters/kmm.txt` against the pipeline prefix (the parent of
`parameters/`), so the fixture's relative paths work without modification.

Requires a `pluma` built with `--with-rust` and the `HAVE_RUST`-define fix
from [FIUBioRG/PluMA#13](https://github.com/FIUBioRG/PluMA/pull/13).

## Regenerating the synthetic data

The fixture is fully deterministic — same seed (xorshift64\*, `0xC0FFEE` /
`0xBADC0DE`) yields the same bytes every run. To regenerate the entire
fixture (training genomes, models, test reads, parameter file, expected
output):

```bash
cargo run --release --example generate_fixture
```

The generator is a Rust example binary (not a unit test); it uses the kmm-rs
library API to train fresh `.model` files from the seeded training FASTAs and
asserts that all 20 test reads classify back to their true class. Use it
when you change the algorithm or want to rebuild the fixture from scratch —
**the example itself is exercised end-to-end through the pluma binary**, not
through this generator.

## Expected output

Tab-separated `(read_index, best_class, log_score)`, e.g.:

```
1	GC_rich	-232.3413
2	GC_rich	-228.17935
...
11	AT_rich	-249.1132
12	AT_rich	-241.8867
...
```

The first 10 reads should all classify as `GC_rich`, the next 10 as `AT_rich`.
The generator fails if any read is misclassified.

## What KMM does

For each model `M_c`, training over a genome `G_c` estimates the conditional
probability `P(b | b_{-1} b_{-2} ... b_{-k})` of nucleotide `b` given the
previous `k` bases. Scoring a read `r` against `M_c` yields
`log P(r | M_c) = Σ_i log P(r_i | r_{i-1} ... r_{i-k})`. The classifier picks
`argmax_c log P(r | M_c)`. See
[movingpictures83/KMM](https://github.com/movingpictures83/KMM) for the
original C++ implementation this Rust port replaces.
