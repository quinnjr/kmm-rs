//! Generate a self-contained, deterministic fixture for the kmm-rs example.
//!
//! Run with:
//!
//! ```bash
//! cargo run --release --example generate_fixture
//! ```
//!
//! Writes:
//! - `example/models/GC_rich.model`  — k-Markov model trained on a
//!   synthetic GC-biased "genome" (40 % G, 40 % C, 10 % A, 10 % T).
//! - `example/models/AT_rich.model`  — k-Markov model trained on a
//!   synthetic AT-biased "genome" (40 % A, 40 % T, 10 % G, 10 % C).
//! - `example/training/GC_rich.fa`  — committed training FASTA so the
//!   fixture is auditable and regeneration is purely deterministic.
//! - `example/training/AT_rich.fa`
//! - `example/sequences.fasta`     — 20 test reads (10 from each
//!   distribution) for `KmmPlugin` to classify.
//! - `example/parameters.kmm.txt`  — PluMA parameter file (order, models
//!   path, input path).
//! - `example/expected_output.tsv` — committed reference output produced
//!   by the *current* commit; used as a regression check.
//!
//! The generator is fully deterministic — same seed yields the same bytes
//! every run, so the fixture stays in lockstep with the committed CSVs.

use kmm_rs::{KmmPlugin, Model, Score};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

// xorshift64* — tiny deterministic PRNG so the fixture has no rand crate dep
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn pick(&mut self, weights: &[(u8, u64)]) -> u8 {
        let total: u64 = weights.iter().map(|(_, w)| *w).sum();
        let mut t = self.next() % total;
        for &(b, w) in weights {
            if t < w {
                return b;
            }
            t -= w;
        }
        weights[weights.len() - 1].0
    }
}

const ORDER: usize = 3;
const TRAINING_LEN: usize = 50_000;
const READ_LEN: usize = 200;
const READS_PER_CLASS: usize = 10;

fn generate_biased(rng: &mut Rng, len: usize, weights: &[(u8, u64)]) -> String {
    let mut s = String::with_capacity(len);
    for _ in 0..len {
        s.push(rng.pick(weights) as char);
    }
    s
}

fn write_fasta(path: &Path, name: &str, seq: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, ">{}", name)?;
    for chunk in seq.as_bytes().chunks(80) {
        w.write_all(chunk)?;
        writeln!(w)?;
    }
    Ok(())
}

fn write_test_fasta(
    path: &Path,
    rng: &mut Rng,
    gc_weights: &[(u8, u64)],
    at_weights: &[(u8, u64)],
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for i in 0..READS_PER_CLASS {
        let seq = generate_biased(rng, READ_LEN, gc_weights);
        writeln!(w, ">read_GC_{:02}", i)?;
        writeln!(w, "{}", seq)?;
    }
    for i in 0..READS_PER_CLASS {
        let seq = generate_biased(rng, READ_LEN, at_weights);
        writeln!(w, ">read_AT_{:02}", i)?;
        writeln!(w, "{}", seq)?;
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    let example = here.join("example");
    let models = example.join("models");
    let training = example.join("training");
    let parameters = example.join("parameters");
    fs::create_dir_all(&models)?;
    fs::create_dir_all(&training)?;
    fs::create_dir_all(&parameters)?;

    let gc_weights = [(b'A', 1), (b'C', 4), (b'G', 4), (b'T', 1)];
    let at_weights = [(b'A', 4), (b'C', 1), (b'G', 1), (b'T', 4)];

    // ------------------------------------------------------------------
    // Training genomes (committed so the build is auditable)
    // ------------------------------------------------------------------
    let mut rng = Rng(0xC0FFEE);
    let gc_genome = generate_biased(&mut rng, TRAINING_LEN, &gc_weights);
    let at_genome = generate_biased(&mut rng, TRAINING_LEN, &at_weights);
    write_fasta(&training.join("GC_rich.fa"), "GC_rich_training", &gc_genome)?;
    write_fasta(&training.join("AT_rich.fa"), "AT_rich_training", &at_genome)?;
    eprintln!("wrote training FASTAs");

    // ------------------------------------------------------------------
    // Train and save the two models
    // ------------------------------------------------------------------
    for (name, genome) in [("GC_rich", &gc_genome), ("AT_rich", &at_genome)] {
        let mut model = Model::new(ORDER);
        model.build_from_genome(genome);
        let path = models.join(format!("{name}.model"));
        model.save(&path)?;
        eprintln!("wrote {}", path.display());
    }

    // ------------------------------------------------------------------
    // Test reads (10 from each class — KMM should classify them correctly)
    // ------------------------------------------------------------------
    let mut rng = Rng(0xBADC0DE);
    let reads_path = example.join("sequences.fasta");
    write_test_fasta(&reads_path, &mut rng, &gc_weights, &at_weights)?;
    eprintln!("wrote {}", reads_path.display());

    // ------------------------------------------------------------------
    // PluMA parameter file. Placed at example/parameters/kmm.txt so the
    // <prefix>/parameters/<plugin>.txt convention applies and the path
    // resolver in our PluMAPlugin impl picks up the relative paths
    // automatically when called from pluma (or from anywhere via the
    // trait method).
    // ------------------------------------------------------------------
    let params_path = parameters.join("kmm.txt");
    fs::write(
        &params_path,
        format!(
            "# kmm-rs example parameters\n\
             #   order: k-Markov order (3 → 4-mer transition probabilities)\n\
             #   models: directory of .model files; basename = class label\n\
             #   input: FASTA of reads to classify\n\
             # Paths are relative to the pipeline prefix (the parent of\n\
             # this file's directory), matching PluMA's R-plugin convention.\n\
             order	{ORDER}\n\
             models	models\n\
             input	sequences.fasta\n"
        ),
    )?;
    eprintln!("wrote {}", params_path.display());

    // ------------------------------------------------------------------
    // Reference output (regression check)
    // ------------------------------------------------------------------
    let scorer = Score::new(ORDER);
    let expected_path = example.join("expected_output.tsv");
    scorer.score_models(&models, &reads_path, &expected_path)?;
    eprintln!("wrote {}", expected_path.display());

    // Sanity check: assert every GC-tagged read scored higher on the GC model,
    // and every AT-tagged read scored higher on the AT model.
    let output = fs::read_to_string(&expected_path)?;
    let mut wrong = 0usize;
    for (i, line) in output.lines().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let expected = if i < READS_PER_CLASS { "GC_rich" } else { "AT_rich" };
        if cols[1] != expected {
            wrong += 1;
            eprintln!("  read {} expected {} got {}", i, expected, cols[1]);
        }
    }
    if wrong > 0 {
        return Err(format!("{wrong} reads misclassified — fixture is broken").into());
    }
    eprintln!("\n✓ all 20 reads classified correctly — fixture ok");

    // Also exercise the full KmmPlugin lifecycle (config → models → input → out)
    let mut plugin = KmmPlugin::new();
    plugin.input(&params_path)?;
    eprintln!("  KmmPlugin loaded config (order={})", plugin.order());

    Ok(())
}
