#![allow(unused)]
use clap::{App, Arg};
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, BufWriter};
use std::thread;
use std::time::Duration;

// TODO: Make it multithreaded via threadpools?
// TODO: Optimise modifying every vec element
// NOTE: This could be done via:
//       Get length of each bed region,
//       assign var to 0
//       the first x amount of elements will have the name according to index one
//       add length to var, then start from that +1 until the next region length,
//       assign all of these to the next index of the bed names

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone)]
struct BedRegion {
    chromosome: String,
    start: i64,
    end: i64,
    name: String,
}

struct OutputType {
    file: BufWriter<File>,
    stdout: bool,
}

#[derive(Debug)]
struct DepthInfo {
    chromosome: String,
    basenumber: i64,
    reads: i64,
    name: String,
}

impl DepthInfo {
    fn add_name(&mut self, name: &str) {
        self.name.push_str(name);
    }
}

fn main() {
    let matches = App::new("Adding names to samtools depth output from bedfiles")
        .version("0.1")
        .author("Me")
        .about(
            "Adds names defined in a bedfile to the output of \"samtools depth\" command output, as the resulting depth file does not contain the name of what genic region the base belongs to.
It also assumes that both the bed file and the depth file are sorted by start position / base position respectively",
        )
        .arg(
            Arg::with_name("depth")
                .short('d')
                .long("depth")
                .value_name("depth")
                .help("Depth file from \"samtools depth\" command to add names to")
                .takes_value(true)
                .multiple(true)
                .required(true),
        )
        .arg(
            Arg::with_name("bedfile")
                .short('b')
                .long("bed")
                .value_name("bedfile")
                .help("bed file containing the start and end values, as well as the names of the regions")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short('o')
                .long("output")
                .value_name("output")
                .help("Either define the output file name or the output will be written to stdout")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("threads")
                .short('t')
                .long("threads")
                .value_name("threads")
                .help("How many threads to use for the program")
                .takes_value(true),
        )
        .get_matches();

    let bedfile = matches.value_of("bedfile").unwrap();
    let depthfiles: Vec<String> = matches
        .get_many::<String>("depth")
        .expect("Womw")
        .map(|s| s.to_string())
        .collect();

    // Read bed regions once, so if there are more than one depth file to look at, no need to read the bed file again
    let bed_regions: Vec<BedRegion> = read_bed(bedfile);

    // Multithread configuration
    let n_threads: usize = matches
        .get_one("threads")
        .unwrap_or(&"3".to_string())
        .parse::<usize>()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads.to_owned() as usize)
        .build_global()
        .unwrap();

    if depthfiles.len() > 1 {
        // If more than one file given, automatically output to different files
        // NOTE: This is the multithreaded version using rayon
        depthfiles.par_iter().for_each(|i| {
            let mut depths: Vec<DepthInfo> = read_depths(i);
            add_name_to_depth(depths.as_mut(), &bed_regions);
            write_depthn(&depths, i, false);
        });
    } else {
        // Only one depth file to look at and write/print, stdout or outputfile if given
        let mut depths = read_depths(&depthfiles[0]);
        add_name_to_depth(&mut depths, &bed_regions);
        if matches.is_present("output") {
            let o = matches.value_of("output").unwrap();
            write_depthn(&depths, o, false);
        } else {
            write_depthn(&depths, "", true);
        }
    }
}

/// Writes the depth info to either stdout or a file
fn write_depthn(depths: &Vec<DepthInfo>, filename: &str, stdout: bool) {
    // Writes the depth file along with the name column
    let mut writer = BufWriter::new(File::create(&format!("{}.depthn", filename)).unwrap());
    let mut stdout_writer = BufWriter::new(std::io::stdout());
    for d in depths {
        if stdout {
            // println!(
            //     "{}\t{}\t{}\t{}",
            //     d.chromosome, d.basenumber, d.reads, d.name
            // );
            writeln!(
                stdout_writer,
                "{}\t{}\t{}\t{}",
                d.chromosome, d.basenumber, d.reads, d.name
            )
            .unwrap();
        } else {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                d.chromosome, d.basenumber, d.reads, d.name
            )
            .unwrap();
        }
    }
    println!("Wrote {f} to file {f}.depthn", f = filename);
}

/// Adds the name of the region to the depth file, based on the bed file
/// input: vector of depth info, vector of bed regions
/// output: mutated original vector of depth info with the name of the region added
fn add_name_to_depth(depths: &mut Vec<DepthInfo>, bed_regions: &Vec<BedRegion>) {
    let bar = ProgressBar::new_spinner();
    bar.enable_steady_tick(Duration::from_millis(10));
    bar.set_message("Modifying Depth files...");
    let mut idx = 0;
    for depth in depths {
        if depth.basenumber <= bed_regions[idx].end {
            depth.name = bed_regions[idx].name.clone();
        } else {
            if idx < bed_regions.len() - 1 {
                idx += 1;
            }
            depth.name = bed_regions[idx].name.clone();
        }
    }
    bar.finish();
}

/// Read a .bed file that contains choromome, start, end and name of region. Should only be read once
fn read_bed(bedfile: &str) -> Vec<BedRegion> {
    let content = File::open(bedfile).expect("Unable to open file");
    let reader = BufReader::new(content);
    let mut bed_regions: Vec<BedRegion> = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split_line = line.split('\t').collect::<Vec<&str>>();
        match split_line.len() {
            4 => {
                let chromosome = split_line[0].to_string();
                let start = split_line[1]
                    .parse::<i64>()
                    .unwrap_or_else(|_| panic!("Unable to parse start value {}", line));
                let end = split_line[2]
                    .parse::<i64>()
                    .unwrap_or_else(|_| panic!("Unable to parse end value {}", line));
                let name = split_line[3].to_string();
                bed_regions.push(BedRegion {
                    chromosome,
                    start,
                    end,
                    name,
                });
            }
            _ => {
                println!(
                    "Error in bed file; incorrectly structured at line: {}. ",
                    line
                );
                std::process::exit(1);
            }
        }
    }
    bed_regions
}

/// Read a .depth file from the output of the samtools depth command
fn read_depths(filename: &str) -> Vec<DepthInfo> {
    let content = File::open(filename).expect("Unable to open file");
    let mut depths: Vec<DepthInfo> = Vec::new();
    let reader = BufReader::new(content);
    let bar = ProgressBar::new_spinner();
    bar.enable_steady_tick(Duration::from_millis(10));
    bar.set_message("Reading Depth files...");
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split_line = line.split('\t');
        let chromosome = split_line.next().unwrap();
        let position = split_line.next().unwrap().parse::<i64>().unwrap();
        let depth = split_line.next().unwrap().parse::<i64>().unwrap();
        depths.push(DepthInfo {
            chromosome: chromosome.to_string(),
            basenumber: position,
            reads: depth,
            name: "".to_string(),
        });
    }
    bar.finish();
    depths
}
