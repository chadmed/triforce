use triforce;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seconds: f32 = if args.len() < 2 {
        60.0
    } else {
        match args[1].parse() {
            Ok(x) => x,
            Err(_) => {
                panic!(
                    "Usage: perf_test <duration>\n\
                        Simulate processing <duration> seconds of microphone input \
                        (default to 60 seconds)."
                )
            }
        }
    };
    let sample_rate = 48000.0;
    let mut inst = triforce::Triforce::with_sample_rate(48000.0);
    let i1: Vec<f32> = (0..1024).map(|x| (x as f32) / 1024.0).collect();
    let i2: Vec<f32> = (0..1024).map(|x| ((x as f32) + 10.0) / 1024.0).collect();
    let i3: Vec<f32> = (0..1024).map(|x| ((x as f32) - 10.0) / 1024.0).collect();
    let mut out: Vec<f32> = vec![0.0; 1024];
    for _ in 0..(sample_rate * seconds / 1024.0) as i32 {
        inst.process_slice(&i1, &i2, &i3, &mut out, 100.0, 1024);
    }
}
