use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use kmm_rs::{Kmers, Model, Score, Sequences};
use std::hint::black_box;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_genome(size: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    (0..size)
        .map(|i| bases[i % 4])
        .collect()
}

fn create_test_fasta(num_seqs: usize, seq_len: usize) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    for i in 0..num_seqs {
        writeln!(file, ">seq{}", i).unwrap();
        let seq = create_test_genome(seq_len);
        writeln!(file, "{}", seq).unwrap();
    }
    file
}

fn bench_kmers_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmers_generation");

    for order in [1, 2, 3, 4].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(order), order, |b, &order| {
            b.iter(|| {
                let kmers = Kmers::new(black_box(order));
                black_box(kmers.get_kmer_list().len())
            })
        });
    }

    group.finish();
}

fn bench_model_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("model_build");

    for genome_size in [1000, 10000, 100000].iter() {
        let genome = create_test_genome(*genome_size);

        group.bench_with_input(
            BenchmarkId::new("order3", genome_size),
            genome_size,
            |b, _| {
                b.iter(|| {
                    let mut model = Model::new(3);
                    model.build_from_genome(black_box(&genome));
                    black_box(model.size())
                })
            },
        );
    }

    group.finish();
}

fn bench_score_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("score_read");

    // Build a test model
    let genome = create_test_genome(10000);
    let mut model = Model::new(3);
    model.build_from_genome(&genome);

    let scorer = Score::new(3);

    for read_len in [100, 500, 1000, 5000].iter() {
        let read = create_test_genome(*read_len);

        group.bench_with_input(BenchmarkId::from_parameter(read_len), read_len, |b, _| {
            b.iter(|| {
                let score = scorer.score_read(black_box(&read), black_box(&model));
                black_box(score)
            })
        });
    }

    group.finish();
}

fn bench_sequences_load(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequences_load");

    for num_seqs in [10, 100, 1000].iter() {
        let file = create_test_fasta(*num_seqs, 500);

        group.bench_with_input(BenchmarkId::from_parameter(num_seqs), num_seqs, |b, _| {
            b.iter(|| {
                let mut seqs = Sequences::new();
                seqs.load(black_box(file.path())).unwrap();
                black_box(seqs.get_sequences().len())
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_kmers_generation,
    bench_model_build,
    bench_score_read,
    bench_sequences_load
);
criterion_main!(benches);
