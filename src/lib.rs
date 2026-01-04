//! k-Markov Models (KMM) Plugin - Rust implementation for PluMA
//!
//! DNA sequence classification using k-th order Markov models.
//! Based on work by Vanessa Aguiar-Pulido, FIU Bioinformatics Research Group.
//!
//! Original C++ implementation: movingpictures83/KMM

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Map base to index (A=0, C=1, G=2, T=3)
#[inline(always)]
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Generate all k-mers of a given order
pub struct Kmers {
    order: usize,
    kmer_list: HashMap<String, usize>,
}

impl Kmers {
    pub fn new(order: usize) -> Self {
        let mut kmers = Kmers {
            order,
            kmer_list: HashMap::with_capacity(4_usize.pow((order + 1) as u32)),
        };
        kmers.generate_kmer_list();
        kmers
    }

    fn generate_kmer_list(&mut self) {
        let k = self.order + 1;
        let total = 4_usize.pow(k as u32);
        let bases = [b'A', b'C', b'G', b'T'];

        // Generate k-mers directly without recursion
        let mut kmer_bytes = vec![b'A'; k];

        for i in 0..total {
            // Convert index to k-mer
            let mut idx = i;
            for j in (0..k).rev() {
                kmer_bytes[j] = bases[idx % 4];
                idx /= 4;
            }

            // Safety: we only use ASCII bytes
            let kmer = unsafe { String::from_utf8_unchecked(kmer_bytes.clone()) };
            self.kmer_list.insert(kmer, i);
        }
    }

    pub fn get_kmer_list(&self) -> &HashMap<String, usize> {
        &self.kmer_list
    }

    pub fn order(&self) -> usize {
        self.order
    }
}

/// Markov model for DNA sequences
pub struct Model {
    order: usize,
    log_probabilities: Vec<f32>,
    size: usize,
}

impl Model {
    pub fn new(order: usize) -> Self {
        let size = 4_usize.pow((order + 1) as u32) + 4_usize.pow(order as u32);
        Model {
            order,
            log_probabilities: vec![f32::NEG_INFINITY; size],
            size,
        }
    }

    /// Load model from file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::open(path.as_ref())
            .map_err(|e| format!("Failed to open model file: {}", e))?;
        let reader = BufReader::new(file);

        let mut log_probs = Vec::new();

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
            let prob: f32 = line
                .trim()
                .parse()
                .map_err(|e| format!("Failed to parse probability: {}", e))?;
            log_probs.push(prob);
        }

        // Determine order from size: size = 4^(k+1) + 4^k
        let mut order = 0_usize;
        let mut found = false;
        for k in 0..=10 {
            let expected_size = 4_usize.pow((k + 1) as u32) + 4_usize.pow(k as u32);
            if expected_size == log_probs.len() {
                order = k;
                found = true;
                break;
            }
            if expected_size > log_probs.len() {
                return Err(format!("Invalid model size: {}", log_probs.len()));
            }
        }

        if !found {
            return Err("Could not determine model order".to_string());
        }

        Ok(Model {
            order,
            log_probabilities: log_probs,
            size: 4_usize.pow((order + 1) as u32) + 4_usize.pow(order as u32),
        })
    }

    /// Build model from genome sequence - optimized version
    pub fn build_from_genome(&mut self, genome: &str) {
        let k = self.order;
        let k1 = k + 1;
        let size_k1 = 4_usize.pow(k1 as u32);
        let size_k = 4_usize.pow(k as u32);

        let mut freqs = vec![0u64; self.size];
        let genome_bytes = genome.as_bytes();
        let len = genome_bytes.len();

        if len < k1 {
            return;
        }

        // Count k+1 mers using sliding window with direct index calculation
        let mut idx: Option<usize> = None;
        let mask_k1 = size_k1 - 1; // For wrapping index

        for i in 0..=len - k1 {
            // Check for N in current window
            let mut has_n = false;
            let mut new_idx = 0_usize;

            if let Some(prev) = idx {
                // Sliding window: shift left by 4 bits and add new base
                let new_base = base_to_idx(genome_bytes[i + k]);
                if let Some(b) = new_base {
                    new_idx = ((prev << 2) & mask_k1) | b;
                } else {
                    has_n = true;
                }
            } else {
                // Calculate full index for first window or after N
                for j in 0..k1 {
                    if let Some(b) = base_to_idx(genome_bytes[i + j]) {
                        new_idx = (new_idx << 2) | b;
                    } else {
                        has_n = true;
                        break;
                    }
                }
            }

            if has_n {
                idx = None;
            } else {
                freqs[new_idx] += 1;
                idx = Some(new_idx);
            }
        }

        // Count k-mers for initial probabilities
        if k > 0 {
            let mask_k = size_k - 1;
            idx = None;

            for i in 0..=len - k {
                let mut has_n = false;
                let mut new_idx = 0_usize;

                if let Some(prev) = idx {
                    let new_base = base_to_idx(genome_bytes[i + k - 1]);
                    if let Some(b) = new_base {
                        new_idx = ((prev << 2) & mask_k) | b;
                    } else {
                        has_n = true;
                    }
                } else {
                    for j in 0..k {
                        if let Some(b) = base_to_idx(genome_bytes[i + j]) {
                            new_idx = (new_idx << 2) | b;
                        } else {
                            has_n = true;
                            break;
                        }
                    }
                }

                if has_n {
                    idx = None;
                } else {
                    freqs[size_k1 + new_idx] += 1;
                    idx = Some(new_idx);
                }
            }
        }

        // Calculate probabilities
        // For (k+1)-mers: P(last char | first k chars)
        for i in 0..size_k {
            let base_idx = i * 4;
            let total: u64 = freqs[base_idx] + freqs[base_idx + 1] + freqs[base_idx + 2] + freqs[base_idx + 3];
            if total > 0 {
                let total_f = total as f32;
                for j in 0..4 {
                    let idx = base_idx + j;
                    if freqs[idx] > 0 {
                        self.log_probabilities[idx] = (freqs[idx] as f32 / total_f).ln();
                    }
                }
            }
        }

        // For k-mers: initial probabilities
        let total_k: u64 = freqs[size_k1..].iter().sum();
        if total_k > 0 {
            let total_f = total_k as f32;
            for i in 0..size_k {
                let idx = size_k1 + i;
                if freqs[idx] > 0 {
                    self.log_probabilities[idx] = (freqs[idx] as f32 / total_f).ln();
                }
            }
        }
    }

    /// Save model to file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let file = File::create(path.as_ref())
            .map_err(|e| format!("Failed to create model file: {}", e))?;
        let mut writer = BufWriter::new(file);

        for prob in &self.log_probabilities {
            writeln!(writer, "{}", prob)
                .map_err(|e| format!("Failed to write probability: {}", e))?;
        }

        Ok(())
    }

    pub fn order(&self) -> usize {
        self.order
    }

    pub fn log_probabilities(&self) -> &[f32] {
        &self.log_probabilities
    }

    pub fn size(&self) -> usize {
        self.size
    }
}

/// DNA sequence collection
pub struct Sequences {
    sequences: Vec<String>,
}

impl Sequences {
    pub fn new() -> Self {
        Sequences {
            sequences: Vec::new(),
        }
    }

    /// Load sequences from FASTA file
    pub fn load<P: AsRef<Path>>(&mut self, path: P) -> Result<(), String> {
        let file = File::open(path.as_ref())
            .map_err(|e| format!("Failed to open sequence file: {}", e))?;
        let reader = BufReader::new(file);

        let mut current_seq = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
            let line = line.trim();

            if line.starts_with('>') {
                if !current_seq.is_empty() {
                    self.sequences.push(std::mem::take(&mut current_seq));
                }
            } else {
                // Filter valid nucleotides - process bytes directly
                for &b in line.as_bytes() {
                    let c = b.to_ascii_uppercase();
                    if matches!(c, b'A' | b'C' | b'G' | b'T' | b'N') {
                        current_seq.push(c as char);
                    }
                }
            }
        }

        if !current_seq.is_empty() {
            self.sequences.push(current_seq);
        }

        Ok(())
    }

    pub fn get_sequences(&self) -> &[String] {
        &self.sequences
    }

    /// Load entire genome as single string
    pub fn load_genome<P: AsRef<Path>>(path: P) -> Result<String, String> {
        let file = File::open(path.as_ref())
            .map_err(|e| format!("Failed to open genome file: {}", e))?;
        let reader = BufReader::new(file);

        let mut genome = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
            let line = line.trim();

            if !line.starts_with('>') {
                for &b in line.as_bytes() {
                    let c = b.to_ascii_uppercase();
                    if matches!(c, b'A' | b'C' | b'G' | b'T' | b'N') {
                        genome.push(c as char);
                    }
                }
            }
        }

        Ok(genome)
    }
}

impl Default for Sequences {
    fn default() -> Self {
        Self::new()
    }
}

/// Score sequences against models
pub struct Score;

impl Score {
    pub fn new(_order: usize) -> Self {
        Score
    }

    /// Score a single read against a model - optimized version
    #[inline]
    pub fn score_read(&self, read: &str, model: &Model) -> f32 {
        let k = model.order();
        let k1 = k + 1;
        let size_k1 = 4_usize.pow(k1 as u32);

        let read_bytes = read.as_bytes();
        let len = read_bytes.len();

        if len < k1 {
            return f32::NEG_INFINITY;
        }

        let log_probs = model.log_probabilities();
        let mask = size_k1 - 1;

        let mut score: f32 = 0.0;
        let mut valid_kmers = 0;
        let mut idx: Option<usize> = None;

        for i in 0..=len - k1 {
            let mut has_n = false;
            let new_idx;

            if let Some(prev) = idx {
                // Sliding window update
                if let Some(b) = base_to_idx(read_bytes[i + k]) {
                    new_idx = ((prev << 2) & mask) | b;
                } else {
                    idx = None;
                    continue;
                }
            } else {
                // Full calculation
                let mut calc_idx = 0_usize;
                for j in 0..k1 {
                    if let Some(b) = base_to_idx(read_bytes[i + j]) {
                        calc_idx = (calc_idx << 2) | b;
                    } else {
                        has_n = true;
                        break;
                    }
                }
                if has_n {
                    continue;
                }
                new_idx = calc_idx;
            }

            let prob = log_probs[new_idx];
            if prob.is_finite() {
                score += prob;
                valid_kmers += 1;
            }

            idx = Some(new_idx);
        }

        if valid_kmers == 0 {
            f32::NEG_INFINITY
        } else {
            score
        }
    }

    /// Score all reads against all models in a directory
    pub fn score_models<P1: AsRef<Path>, P2: AsRef<Path>, P3: AsRef<Path>>(
        &self,
        models_dir: P1,
        input_file: P2,
        output_file: P3,
    ) -> Result<(), String> {
        // Load models
        let mut models: Vec<(String, Model)> = Vec::new();

        let entries = fs::read_dir(models_dir.as_ref())
            .map_err(|e| format!("Failed to read models directory: {}", e))?;

        for entry in entries {
            let entry = entry.map_err(|e| format!("Failed to read directory entry: {}", e))?;
            let path = entry.path();

            if path.is_file() {
                if let Some(name) = path.file_stem() {
                    let name = name.to_string_lossy().to_string();
                    match Model::from_file(&path) {
                        Ok(model) => {
                            println!("[KMM] Loaded model: {}", name);
                            models.push((name, model));
                        }
                        Err(e) => {
                            eprintln!("[KMM] Warning: Could not load {}: {}", path.display(), e);
                        }
                    }
                }
            }
        }

        if models.is_empty() {
            return Err("No valid models found".to_string());
        }

        // Load sequences
        let mut sequences = Sequences::new();
        sequences.load(input_file.as_ref())?;
        let reads = sequences.get_sequences();

        println!("[KMM] Scoring {} reads against {} models", reads.len(), models.len());

        // Score each read against all models
        let out_file = File::create(output_file.as_ref())
            .map_err(|e| format!("Failed to create output file: {}", e))?;
        let mut writer = BufWriter::new(out_file);

        for (i, read) in reads.iter().enumerate() {
            let mut best_score = f32::NEG_INFINITY;
            let mut best_model = String::new();

            for (name, model) in &models {
                let score = self.score_read(read, model);
                if score > best_score {
                    best_score = score;
                    best_model = name.clone();
                }
            }

            writeln!(writer, "{}\t{}\t{}", i + 1, best_model, best_score)
                .map_err(|e| format!("Failed to write result: {}", e))?;
        }

        Ok(())
    }
}

/// KMM Plugin - main interface
pub struct KmmPlugin {
    order: usize,
    models_path: String,
    input_file: String,
}

impl Default for KmmPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl KmmPlugin {
    pub fn new() -> Self {
        KmmPlugin {
            order: 3,
            models_path: String::new(),
            input_file: String::new(),
        }
    }

    /// Read configuration from file
    pub fn input<P: AsRef<Path>>(&mut self, config_file: P) -> Result<(), String> {
        let file = File::open(config_file.as_ref())
            .map_err(|e| format!("Failed to open config file: {}", e))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read config: {}", e))?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 2 {
                match parts[0].to_lowercase().as_str() {
                    "order" => {
                        self.order = parts[1]
                            .parse()
                            .map_err(|_| "Invalid order value".to_string())?;
                    }
                    "models" | "pathtmodels" => {
                        self.models_path = parts[1].to_string();
                    }
                    "input" | "inputname" => {
                        self.input_file = parts[1].to_string();
                    }
                    _ => {}
                }
            }
        }

        Ok(())
    }

    pub fn run(&self) {
        // Scoring happens in output phase for compatibility with PluMA pattern
    }

    pub fn output<P: AsRef<Path>>(&self, output_file: P) -> Result<(), String> {
        let scorer = Score::new(self.order);
        scorer.score_models(&self.models_path, &self.input_file, output_file)
    }

    pub fn order(&self) -> usize {
        self.order
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmers_generation() {
        let kmers = Kmers::new(1); // 2-mers
        let list = kmers.get_kmer_list();

        assert_eq!(list.len(), 16); // 4^2 = 16
        assert!(list.contains_key("AA"));
        assert!(list.contains_key("AT"));
        assert!(list.contains_key("TT"));
    }

    #[test]
    fn test_kmers_order_0() {
        let kmers = Kmers::new(0); // 1-mers
        let list = kmers.get_kmer_list();

        assert_eq!(list.len(), 4);
        assert!(list.contains_key("A"));
        assert!(list.contains_key("C"));
        assert!(list.contains_key("G"));
        assert!(list.contains_key("T"));
    }

    #[test]
    fn test_model_creation() {
        let model = Model::new(2);
        // size = 4^3 + 4^2 = 64 + 16 = 80
        assert_eq!(model.size(), 80);
    }

    #[test]
    fn test_sequences_loading() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "TTAA").unwrap();

        let mut seqs = Sequences::new();
        seqs.load(file.path()).unwrap();

        assert_eq!(seqs.get_sequences().len(), 2);
        assert_eq!(seqs.get_sequences()[0], "ACGT");
        assert_eq!(seqs.get_sequences()[1], "TTAA");
    }

    #[test]
    fn test_score_simple() {
        let mut model = Model::new(1);
        // Set uniform probabilities for testing
        let size_k1 = 16; // 4^2
        for i in 0..size_k1 {
            model.log_probabilities[i] = -1.0; // ln(e^-1)
        }

        let scorer = Score::new(1);
        let score = scorer.score_read("ACGT", &model);

        assert!(score.is_finite());
        assert!(score < 0.0); // Log probs are negative
    }

    #[test]
    fn test_model_build() {
        let mut model = Model::new(2);
        model.build_from_genome("ACGTACGTACGT");

        // Should have counted some k-mers
        let probs = model.log_probabilities();
        let finite_count = probs.iter().filter(|p| p.is_finite()).count();
        assert!(finite_count > 0);
    }
}
