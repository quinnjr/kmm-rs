//! KMM CLI - k-Markov Models for DNA sequence classification

use kmm_rs::{KmmPlugin, Model, Score, Sequences};
use std::env;
use std::process;

fn print_usage(prog: &str) {
    eprintln!("Usage: {} <command> [options]", prog);
    eprintln!();
    eprintln!("Commands:");
    eprintln!("  score <order> <models_dir> <input.fasta> <output.txt>");
    eprintln!("      Score sequences against models");
    eprintln!();
    eprintln!("  build <order> <genome.fasta> <output.model>");
    eprintln!("      Build a model from a genome");
    eprintln!();
    eprintln!("  plugin <config.txt> <output.txt>");
    eprintln!("      Run as PluMA plugin");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage(&args[0]);
        process::exit(1);
    }

    match args[1].as_str() {
        "score" => {
            if args.len() != 6 {
                eprintln!("Usage: {} score <order> <models_dir> <input.fasta> <output.txt>", args[0]);
                process::exit(1);
            }

            let order: usize = args[2].parse().unwrap_or_else(|_| {
                eprintln!("Invalid order value");
                process::exit(1);
            });

            let scorer = Score::new(order);
            if let Err(e) = scorer.score_models(&args[3], &args[4], &args[5]) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }

            println!("[KMM] Results written to {}", args[5]);
        }

        "build" => {
            if args.len() != 5 {
                eprintln!("Usage: {} build <order> <genome.fasta> <output.model>", args[0]);
                process::exit(1);
            }

            let order: usize = args[2].parse().unwrap_or_else(|_| {
                eprintln!("Invalid order value");
                process::exit(1);
            });

            println!("[KMM] Building order-{} model from {}", order, args[3]);

            let genome = match Sequences::load_genome(&args[3]) {
                Ok(g) => g,
                Err(e) => {
                    eprintln!("Error loading genome: {}", e);
                    process::exit(1);
                }
            };

            println!("[KMM] Loaded genome: {} bp", genome.len());

            let mut model = Model::new(order);
            model.build_from_genome(&genome);

            if let Err(e) = model.save(&args[4]) {
                eprintln!("Error saving model: {}", e);
                process::exit(1);
            }

            println!("[KMM] Model saved to {}", args[4]);
        }

        "plugin" => {
            if args.len() != 4 {
                eprintln!("Usage: {} plugin <config.txt> <output.txt>", args[0]);
                process::exit(1);
            }

            let mut plugin = KmmPlugin::new();

            if let Err(e) = plugin.input(&args[2]) {
                eprintln!("Error reading config: {}", e);
                process::exit(1);
            }

            plugin.run();

            if let Err(e) = plugin.output(&args[3]) {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        }

        _ => {
            print_usage(&args[0]);
            process::exit(1);
        }
    }
}
