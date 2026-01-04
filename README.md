# kmm-rs

Rust implementation of k-Markov Models (KMM) for DNA sequence classification (PluMA plugin).

Based on the original C++ implementation by [movingpictures83/KMM](https://github.com/movingpictures83/KMM) and work by Vanessa Aguiar-Pulido at FIU Bioinformatics Research Group.

## Description

KMM uses k-th order Markov models to classify DNA sequences. It can:
1. Build probabilistic models from reference genomes
2. Score query sequences against multiple models
3. Classify sequences to the best-matching model

## Building

```bash
cargo build --release
```

## Testing

```bash
cargo test
```

## Usage

### Build a model from a genome
```bash
cargo run --release -- build <order> <genome.fasta> <output.model>
```

Example:
```bash
cargo run --release -- build 3 ecoli.fasta ecoli.model
```

### Score sequences against models
```bash
cargo run --release -- score <order> <models_dir> <input.fasta> <output.txt>
```

### Run as PluMA plugin
```bash
cargo run --release -- plugin <config.txt> <output.txt>
```

Config file format:
```
order 3
models /path/to/models/
input /path/to/sequences.fasta
```

## Algorithm

The k-Markov model stores:
- Conditional probabilities P(b_{k+1} | b_1...b_k) for all (k+1)-mers
- Initial probabilities P(b_1...b_k) for all k-mers

For sequence classification:
1. Calculate the log-likelihood of each sequence under each model
2. Assign the sequence to the model with the highest score

## Input Format

FASTA format for both genomes (model building) and query sequences (scoring).

```fasta
>sequence1
ACGTACGTACGT
>sequence2
TTAACCGGTTAA
```

## Output Format

Tab-separated file with columns: ReadID, BestModel, Score

## License

MIT License
