//! The `ultra` command (`main.cpp`).

use std::fs::File;
use std::io::{self, BufWriter, Write};

use ultra::cli::{
    large_tune_settings, medium_tune_settings, small_tune_settings, tune_settings_for_path,
    Settings,
};
use ultra::mask::output_masked_fasta;
use ultra::ultra::Ultra;
use ultra::writers::{fmt_f, fmt_g};

/// C++ `main_wrapper`. Returns the process exit status.
fn main_wrapper(argv: &[String]) -> i32 {
    // Prepare settings
    let mut settings = Settings::new();
    settings.prepare_settings();
    if !settings.parse_input(argv) {
        exit(0);
    }

    settings.assign_settings();
    settings.print_memory_usage();
    if settings.show_memory {
        exit(0);
    }

    // Tuning stuff.
    let mut coverage: Vec<u64> = Vec::new();
    let mut shuffled_coverage: Vec<u64> = Vec::new();
    let mut seq_length: u64 = 0;
    if settings.tune {
        // Get the parameters we are going to run with
        let param_strings: Vec<String>;
        if settings.tune_param_path.is_empty() {
            if settings.tune_medium {
                param_strings = medium_tune_settings();
            } else if settings.tune_large {
                param_strings = large_tune_settings();
            } else {
                param_strings = small_tune_settings();
            }
        } else {
            param_strings = tune_settings_for_path(&settings.tune_param_path);
        }

        // Run ULTRA on each parameter set
        println!("Performing parameter search...");
        println!("Coverage, False Discovery Rate, Parameters");

        let mut best_coverage_index: i32 = -1;
        let mut best_coverage: f64 = 0.0;
        for arg_string in &param_strings {
            let mut search_setting = Settings::new();
            search_setting.prepare_settings();
            search_setting.set_multi_option();
            if !search_setting.parse_multi_input(argv, arg_string) {
                exit(0);
            }
            search_setting.assign_settings();

            if !search_setting.tune_indels {
                search_setting.max_insert = 0;
                search_setting.max_delete = 0;
            }
            search_setting.suppress_out = true;
            search_setting.hide_settings = true;
            search_setting.produce_mask = true;
            search_setting.no_split = true;
            search_setting.max_split = 0;
            search_setting.ultra_out = false;
            search_setting.json_out = false;
            search_setting.bed_out = false;

            // We first get the normal coverage
            let mut ultra = Ultra::new(&search_setting);
            ultra.analyze_file();
            ultra.output_repeats(true);
            seq_length = reader_total_seq_length(&ultra);
            let cvg_tmp = ultra.coverage();
            let real_coverage = cvg_tmp as f64 / seq_length as f64;
            coverage.push(cvg_tmp);
            drop(ultra);

            // Next we get our shuffled coverage
            let mut ultra = Ultra::new(&search_setting);
            ultra.shuffle_sequence = true;
            ultra.analyze_file();
            ultra.output_repeats(true);
            let cvg_tmp = ultra.coverage();
            let false_coverage = cvg_tmp as f64 / seq_length as f64;
            shuffled_coverage.push(cvg_tmp);
            drop(ultra);

            let mut fdr = 1.0;
            if real_coverage > 0.0 {
                fdr = false_coverage / real_coverage;
            }

            if fdr <= settings.tune_fdr && real_coverage > best_coverage {
                best_coverage = real_coverage;
                best_coverage_index = coverage.len() as i32 - 1;
            }

            println!(
                "({}/{}): {}, {}, {}",
                coverage.len(),
                param_strings.len(),
                fmt_f(real_coverage, 5),
                fmt_f(fdr, 5),
                arg_string
            );
        }

        println!("----------------------------");
        if best_coverage_index >= 0 {
            let best = best_coverage_index as usize;
            let real_coverage = coverage[best] as f64 / seq_length as f64;
            let false_coverage = shuffled_coverage[best] as f64 / seq_length as f64;
            let fdr = false_coverage / real_coverage;

            println!(
                "Best coverage within FDR limit ({}): {}, {}, {}",
                fmt_f(settings.tune_fdr, 3),
                fmt_f(real_coverage, 5),
                fmt_f(fdr, 5),
                param_strings[best]
            );

            settings = Settings::new();
            settings.prepare_settings();
            settings.set_multi_option();
            if !settings.parse_multi_input(argv, &param_strings[best]) {
                exit(0);
            }
            settings.assign_settings();
        } else {
            println!(
                "No parameters found within FDR limit ({}).",
                fmt_f(settings.tune_fdr, 3)
            );
            exit(0);
        }

        if settings.tune_only {
            exit(0);
        }
    }

    // Perform the actual run
    let mut ultra = Ultra::new(&settings);
    ultra.analyze_file();
    ultra.output_repeats(true);
    let true_coverage: u64 = ultra.coverage();
    // Check to see if we are making a masked file
    if settings.produce_mask {
        write_mask(&settings, &mut ultra);
    }
    if !settings.disable_summary {
        seq_length = reader_total_seq_length(&ultra);
        let coverage_ratio = true_coverage as f64 / seq_length as f64;
        println!("----------------------------");
        println!(
            "Annotation Coverage: {} ({} / {})",
            fmt_f(coverage_ratio, 4),
            true_coverage,
            seq_length
        );
    }

    drop(ultra);
    if settings.estimate_fdr {
        settings.suppress_out = true;
        settings.hide_settings = true;
        settings.produce_mask = true;
        settings.no_split = true;
        settings.max_split = 0;
        settings.ultra_out = false;
        settings.json_out = false;
        settings.bed_out = false;
        settings.out_file = String::new();
        let mut ultra = Ultra::new(&settings);
        ultra.shuffle_sequence = true;
        ultra.analyze_file();
        ultra.output_repeats(true);
        let shuff_coverage = ultra.coverage();
        let fdr: f32 = shuff_coverage as f32 / true_coverage as f32;
        println!("Estimated false discovery rate: {}", fmt_g(fdr as f64));
        // The C++ never deletes this instance.
        drop(ultra);
    }

    0
}

fn reader_total_seq_length(ultra: &Ultra) -> u64 {
    ultra
        .reader
        .as_ref()
        .expect("reader")
        .fasta_reader
        .total_seq_length
}

/// The masked FASTA output of `main`. If the mask file cannot be opened this
/// prints an error and exits with status 255; upstream 1.2.2 passes the null
/// `FILE *` on and crashes.
fn write_mask(settings: &Settings, ultra: &mut Ultra) {
    let f = match File::create(&settings.mask_file) {
        Ok(f) => f,
        Err(_) => {
            let _ = io::stdout().flush();
            eprintln!("Unable to open mask file {}", settings.mask_file);
            std::process::exit(255);
        }
    };
    let mut w = BufWriter::new(f);
    let result = output_masked_fasta(
        &settings.in_file,
        &mut w,
        &mut ultra.masks_for_seq,
        settings.mask_with_n,
    )
    .and_then(|_| w.flush());
    if let Err(e) = result {
        let _ = io::stdout().flush();
        eprintln!("Unable to write mask file {}: {}", settings.mask_file, e);
        std::process::exit(255);
    }
}

/// C `exit`: flushes stdout first.
fn exit(code: i32) -> ! {
    let _ = io::stdout().flush();
    std::process::exit(code)
}

fn main() {
    let argv: Vec<String> = std::env::args_os()
        .map(|a| a.to_string_lossy().into_owned())
        .collect();
    let r = main_wrapper(&argv);
    let _ = io::stdout().flush();
    std::process::exit(r);
}
