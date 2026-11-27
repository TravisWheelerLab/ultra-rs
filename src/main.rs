//! The `ultra` command (`main.cpp`).

use std::fs::File;
use std::io::{self, BufWriter, Write};

use ultra_repeats::cli::{
    large_tune_settings, medium_tune_settings, small_tune_settings, tune_settings_for_path,
    Settings,
};
use ultra_repeats::mask::output_masked_fasta;
use ultra_repeats::ultra::Ultra;
use ultra_repeats::writers::{fmt_f, fmt_g, on_write_error};

/// `println!` that ends the program quietly on a closed pipe (and reports
/// other write errors) instead of panicking; see `writers::on_write_error`.
macro_rules! say {
    ($($arg:tt)*) => {
        ultra_repeats::writers::print_stdout(&format!("{}\n", format_args!($($arg)*)))
    };
}

/// Runs `ultra` over its input and writes the remaining repeats. A write
/// error ends the program as [`on_write_error`] says, before `ultra` is
/// dropped, so its buffered output is not written after the error.
fn run(ultra: &mut Ultra) {
    if let Err(e) = ultra.analyze_file().and_then(|()| ultra.output_repeats(true)) {
        on_write_error(&e);
    }
}

/// Fresh settings, with `upstream_compat` taken from the environment
/// variable `ULTRA_UPSTREAM_COMPAT` (see the field's documentation).
fn new_settings() -> Settings {
    let mut settings = Settings::new();
    settings.upstream_compat = std::env::var_os("ULTRA_UPSTREAM_COMPAT")
        .is_some_and(|v| !v.is_empty() && v != "0");
    settings
}

/// C++ `main_wrapper`. Returns the process exit status.
fn main_wrapper(argv: &[String]) -> i32 {
    // Prepare settings
    let mut settings = new_settings();
    settings.prepare_settings();
    if !settings.parse_input(argv) {
        exit(ultra_repeats::cli::usage_error_code());
    }

    settings.assign_settings();
    if let Err(e) = settings.print_memory_usage() {
        on_write_error(&e);
    }
    if settings.show_memory {
        exit(0);
    }
    check_memory(&mut settings);

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
        say!("Performing parameter search...");
        say!("Coverage, False Discovery Rate, Parameters");

        let mut best_coverage_index: i32 = -1;
        let mut best_coverage: f64 = 0.0;
        for arg_string in &param_strings {
            let mut search_setting = new_settings();
            search_setting.prepare_settings();
            search_setting.set_multi_option();
            if !search_setting.parse_multi_input(argv, arg_string) {
                exit(ultra_repeats::cli::usage_error_code());
            }
            search_setting.assign_settings();
            // The search runs write nothing, so they must not open (and so
            // empty) the -o file, as the C++ did.
            search_setting.out_file.clear();

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
            search_setting.bpaf_out = false;

            // We first get the normal coverage
            let mut ultra = Ultra::new(&search_setting);
            run(&mut ultra);
            seq_length = reader_total_seq_length(&ultra);
            let cvg_tmp = ultra.coverage();
            let real_coverage = cvg_tmp as f64 / seq_length as f64;
            coverage.push(cvg_tmp);
            drop(ultra);

            // Next we get our shuffled coverage
            let mut ultra = Ultra::new(&search_setting);
            ultra.shuffle_sequence = true;
            run(&mut ultra);
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

            say!(
                "({}/{}): {}, {}, {}",
                coverage.len(),
                param_strings.len(),
                fmt_f(real_coverage, 5),
                fmt_f(fdr, 5),
                arg_string
            );
        }

        say!("----------------------------");
        if best_coverage_index >= 0 {
            let best = best_coverage_index as usize;
            let real_coverage = coverage[best] as f64 / seq_length as f64;
            let false_coverage = shuffled_coverage[best] as f64 / seq_length as f64;
            let fdr = false_coverage / real_coverage;

            say!(
                "Best coverage within FDR limit ({}): {}, {}, {}",
                fmt_f(settings.tune_fdr, 3),
                fmt_f(real_coverage, 5),
                fmt_f(fdr, 5),
                param_strings[best]
            );

            settings = new_settings();
            settings.prepare_settings();
            settings.set_multi_option();
            if !settings.parse_multi_input(argv, &param_strings[best]) {
                exit(ultra_repeats::cli::usage_error_code());
            }
            settings.assign_settings();
        } else {
            say!(
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
    run(&mut ultra);
    let true_coverage: u64 = ultra.coverage();
    // Check to see if we are making a masked file
    if settings.produce_mask {
        write_mask(&settings, &mut ultra);
    }
    if !settings.disable_summary {
        seq_length = reader_total_seq_length(&ultra);
        // an empty input: 0, where the C++ prints -nan
        let coverage_ratio = if seq_length == 0 && !settings.upstream_compat { 0.0 } else { true_coverage as f64 / seq_length as f64 };
        say!("----------------------------");
        say!(
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
        settings.bpaf_out = false;
        settings.out_file = String::new();
        let mut ultra = Ultra::new(&settings);
        ultra.shuffle_sequence = true;
        run(&mut ultra);
        let shuff_coverage = ultra.coverage();
        let fdr: f32 = shuff_coverage as f32 / true_coverage as f32;
        say!("Estimated false discovery rate: {}", fmt_g(fdr as f64));
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

/// Physical memory in bytes, where the system says (Linux `/proc/meminfo`).
fn physical_memory() -> Option<u64> {
    let text = std::fs::read_to_string("/proc/meminfo").ok()?;
    let kb = text.lines().find_map(|l| l.strip_prefix("MemTotal:"))?.trim().strip_suffix("kB")?.trim().parse::<u64>().ok()?;
    kb.checked_mul(1024)
}

/// Stop with a message when the settings need more memory than the machine has (e.g.
/// a very large `--period` or `--overlap`); the allocation would abort the process.
fn check_memory(settings: &mut Settings) {
    let need = settings.memory_usage_bytes();
    let have = physical_memory();
    // the estimate can run to about twice what is used, so only refuse well past it
    if need.is_none() || matches!((need, have), (Some(n), Some(h)) if n / 2 > h) {
        const GB: f64 = 1024.0 * 1024.0 * 1024.0;
        let _ = io::stdout().flush();
        match (need, have) {
            (Some(n), Some(h)) => eprintln!(
                "These settings need about {:.1} GB of memory and this machine has {:.1} GB; lower --period, --win_size, --overlap, --windows or --threads",
                n as f64 / GB,
                h as f64 / GB
            ),
            _ => eprintln!("These settings need more memory than can be addressed; lower --period, --win_size, --overlap or --windows"),
        }
        std::process::exit(ultra_repeats::cli::usage_error_code());
    }
}

/// C `exit`: flushes stdout first.
fn exit(code: i32) -> ! {
    let _ = io::stdout().flush();
    std::process::exit(code)
}

fn main() {
    // Arguments that are not UTF-8 (possible in Linux file names) used to be
    // converted with replacement characters, so the input could not be
    // found, or output went to a file with another name.
    let mut argv: Vec<String> = Vec::new();
    for a in std::env::args_os() {
        match a.into_string() {
            Ok(s) => argv.push(s),
            Err(a) => {
                eprintln!("Argument '{}' is not valid UTF-8; ULTRA needs UTF-8 arguments and file names", a.to_string_lossy());
                std::process::exit(255);
            }
        }
    }
    let r = main_wrapper(&argv);
    let _ = io::stdout().flush();
    std::process::exit(r);
}
