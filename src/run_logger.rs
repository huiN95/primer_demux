use crate::cli::Cli;
use metrics;
use metrics_exporter_prometheus::{PrometheusBuilder, PrometheusHandle};
// use once_cell::sync::Lazy;
// use std::sync::Mutex;
// use std::thread;
// use std::time::Duration;
use std::{
    error::Error,
    fs::{create_dir_all, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};
use tracing_subscriber;
// static LOG_GUARD: Lazy<Mutex<Option<tracing_appender::non_blocking::WorkerGuard>>> =
//     Lazy::new(|| Mutex::new(None));
// use std::{fs, path::Path};
// use tracing_subscriber::fmt::time::ChronoLocal;
use tracing_subscriber::{fmt, layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

pub fn init_tracing_log(cli: &Cli) -> tracing_appender::non_blocking::WorkerGuard {
    let log_path = Path::new(&cli.log_folder);

    // Create recursively if it does not exist
    if !log_path.exists() {
        create_dir_all(log_path).unwrap();
    }
    let run_log_path = log_path.join(format!("{}_primer_demux_run.log", cli.log_folder));

    // let run_log_filepath = format!("{}_primer_demux_run.log", log_path);
    let log_file = File::create(&run_log_path).unwrap_or_else(|e| {
        panic!(
            "create run log file error: {} ({})",
            run_log_path.display(),
            e
        )
    });
    // print!("run log file: {}", &run_log_filepath);
    let (non_blocking, _guard) = tracing_appender::non_blocking(log_file);
    let filter = if let Some(s) = &cli.log {
        EnvFilter::try_new(s).unwrap()
    } else {
        // Default when no parameters or environment variables are passed
        EnvFilter::new("info")
    };
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .with_timer(tracing_subscriber::fmt::time::ChronoLocal::rfc_3339())
                .with_ansi(false)
                .with_writer(non_blocking)
                .with_target(true)
                .with_thread_ids(true),
        )
        .try_init()
        .ok();
    // tracing_subscriber::fmt()
    //     // .with_timer(UtcTime::rfc_3339())
    //     .with_timer(ChronoLocal::rfc_3339())
    //     .with_ansi(false) // remove color info
    //     .with_writer(non_blocking)
    //     .try_init()
    //     .expect("tracing_subscriber already initialized!");
    _guard
    // *LOG_GUARD.lock().unwrap() = Some(_guard);
}

pub struct MetricsGuard {
    handle: PrometheusHandle,
    path: PathBuf,
}

impl Drop for MetricsGuard {
    fn drop(&mut self) {
        // collect metrics
        let report = self.handle.render();

        // try to write; if failed, only print log, do not panic
        if let Err(e) = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.path)
            .and_then(|mut f| f.write_all(report.as_bytes()))
        {
            eprintln!("metrics dump failed: {e}");
        }
    }
}

/// Called at program start, returns a MetricsGuard.
pub fn init_metrics<P: AsRef<Path>>(dir: P) -> Result<MetricsGuard, Box<dyn Error>> {
    // 1. Build recorder
    let recorder = PrometheusBuilder::new()
        .set_buckets(
            vec![
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 13.0, 21.0, 34.0, 55.0, 100.0,
                1000.0,
            ]
            .as_slice(),
        )
        .unwrap()
        .build_recorder();
    let handle = recorder.handle();
    // 2. Register as global. If already registered, return Err instead of panic.
    metrics::set_global_recorder(Box::new(recorder)).unwrap();

    // 3. Generate log path
    // let path = dir.as_ref().join("primer_metrics.prom"); // can change name/add timestamp as needed
    // let mut f = std::fs::File::create(&path).unwrap();
    let in_path = dir.as_ref();
    let out_path = if in_path.is_dir() {
        // Passed a directory
        in_path.join("_primer_metrics.log")
    } else {
        // Passed a file; change filename
        let stem = in_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("metrics");

        let mut buf = PathBuf::from(in_path);
        buf.set_file_name(format!("{stem}_primer_metrics.log"));
        buf
    };
    Ok(MetricsGuard {
        handle,
        path: out_path,
    })
}
