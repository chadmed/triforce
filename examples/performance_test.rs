use triforce;

fn usage<T, A>(_: T) -> A {
    panic!(
        "Usage: perf_test <duration> <blocksize>\n\
         Simulate processing <duration> seconds of microphone input \
         (default to 60 seconds), \
         using blocks of <blocksize> samples (defaults to 1024)."
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seconds: f32 = if args.len() < 2 {
        60.0
    } else {
        args[1].parse().unwrap_or_else(usage)
    };
    let blocksize: usize = if args.len() < 3 {
        1024
    } else {
        args[2].parse().unwrap_or_else(usage)
    };
    let sample_rate = 48000.0;
    let mut inst = triforce::Triforce::with_sample_rate(48000.0);
    let i1: Vec<f32> = (0..blocksize).map(|x| (x as f32) / 1024.0).collect();
    let i2: Vec<f32> = (0..blocksize)
        .map(|x| ((x as f32) + 10.0) / 1024.0)
        .collect();
    let i3: Vec<f32> = (0..blocksize)
        .map(|x| ((x as f32) - 10.0) / 1024.0)
        .collect();
    let mut out: Vec<f32> = Vec::new();
    out.resize(blocksize, 0.0);
    for _ in 0..(sample_rate * seconds / blocksize as f32) as i32 {
        inst.process_slice(&i1, &i2, &i3, &mut out, 100.0, blocksize);
    }
}
