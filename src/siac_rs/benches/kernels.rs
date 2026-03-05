use criterion::{criterion_group, criterion_main, Criterion};

fn noop_benchmark(_c: &mut Criterion) {}

criterion_group!(benches, noop_benchmark);
criterion_main!(benches);
