//! Using ULTRA as a library: find repeats in an in-memory string with
//! `Ultra::find_repeats_in_string`. Port of upstream
//! `examples/library_example.cpp`; run with
//! `cargo run --release --example library_example`.

use ultra_repeats::cli::{string_to_args, Settings};
use ultra_repeats::ultra::Ultra;

fn print_repeats(ultra: &mut Ultra, s: &str) {
    println!("Finding repeats in string ``{}''", s);
    let repeats = ultra
        .find_repeats_in_string(s)
        .expect("string fits the model");

    println!("Start Length Pattern");
    for r in &repeats {
        println!(
            "{} {} {} {}",
            r.window_start,
            r.repeat_length,
            r.repeat_period,
            String::from_utf8_lossy(&r.get_consensus())
        );
    }
    // The caller owns the returned repeats; they are dropped here.
    println!("---------");
}

fn main() {
    println!("Starting test of library.");

    // These are the settings we will use with ULTRA
    let mut settings = Settings::new();
    settings.run_without_reader = true; // This must be set to true when using as lib
    settings.window_size = 1000000; // window size should be equal to the largest sequence you will analyze
    settings.overlap = 0; // overlap should be turned off

    // We can also pass in an argument string to settings. As in the C++
    // example, the first token is taken as the program name, so "-p" is
    // swallowed, "25" becomes the input file name and only "-i 3 -d 3" apply.
    // The first word stands in for the program name.
    let argv = string_to_args("ultra -p 25 -i 3 -d 3");
    settings.prepare_settings(); // This should be called before running parse_input()
    if !settings.parse_input(&argv) {
        std::process::exit(0);
    }

    // Finally, we call assign_settings()
    settings.assign_settings();

    // We create a reusable ULTRA object
    let mut ultra = Ultra::new(&settings);

    print_repeats(&mut ultra, "aaaaaaaaaaaaaaaaaaaaaaaaaa");
    print_repeats(
        &mut ultra,
        "aggtaaggtaaggtaaggtaaggtaagcggtataacatacagatctgactactactactactactactactactac",
    );

    println!("All positions reported are 0 indexed.");
}
