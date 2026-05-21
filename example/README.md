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
└── README.md                 (this file)
```

## Regenerating

The fixture is fully deterministic — the same seed (xorshift64*, `0xC0FFEE` /
`0xBADC0DE`) yields the same bytes every run. To regenerate:

```bash
cargo run --release --example generate_fixture
```

The generator also asserts that all 20 reads classify correctly against their
true class; if you break the algorithm, it fails fast.

## Running the example

### Via the CLI binary (`kmm`)

```bash
# From the example/ directory so the relative paths resolve:
cd example
cargo run --release --bin kmm parameters/kmm.txt /tmp/kmm.out
diff /tmp/kmm.out expected_output.tsv
```

### Via PluMA (real Rust plugin)

Once [`FIUBioRG/PluMA#13`](https://github.com/FIUBioRG/PluMA/pull/13) lands so
`pluma` is built with `-DHAVE_RUST`, drop the `.so` and `Cargo.toml` under
`PluMA/plugins/KMM/`:

```bash
cd /path/to/PluMA/plugins
mkdir -p KMM
ln -sfn /path/to/kmm-rs/target/release/libkmm_rs.so KMM/libKMMPlugin.so
ln -sfn /path/to/kmm-rs/Cargo.toml KMM/Cargo.toml
```

…and run a one-Plugin pipeline:

```
# config.txt
Prefix /path/to/kmm-rs/example/
Plugin KMM  inputfile parameters/kmm.txt  outputfile /tmp/kmm.out
```

```bash
pluma config.txt
diff /tmp/kmm.out /path/to/kmm-rs/example/expected_output.tsv
```

The plugin's `PluMAPlugin::input()` impl resolves the `models` and `input`
paths inside `parameters/kmm.txt` against the pipeline prefix (the parent of
`parameters/`), so no absolute paths or `$cwd` tricks are needed.

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
