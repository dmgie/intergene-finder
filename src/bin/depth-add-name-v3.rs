#![allow(unused)]
use clap::{App, Arg};
use std::boxed;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, BufWriter};

// TODO: Make a nice progress loading bar to show how far through the vec it is
//       or have ... that increase every second to show it is still running
// TODO: Optimise modifying every vec element
// NOTE: This could be done via:
//       Get length of each bed region,
//       assign var to 0
//       the first x amount of elements will have the name according to index one
//       add length to var, then start from that +1 until the next region length,
//       assign all of these to the next index of the bed names

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
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
        .get_matches();

    let bedfile = matches.value_of("bedfile").unwrap();
    let depthfiles = matches.values_of("depth").unwrap().collect::<Vec<&str>>();

    // Read bed regions once, so if there are more than one depth file to look at, no need to read the bed file again
    let bed_regions: Vec<BedRegion> = read_bed(bedfile);

    if matches.value_of("depth").unwrap().len() > 1 {
        // If more than one file given, automatically output to different files
        for i in depthfiles {
            println!("Reading and manipulating depth file {}", i);
            let mut depths: Vec<DepthInfo> = read_depths(i);
            add_name_to_depth(depths.as_mut(), &bed_regions);
            let mut writer = BufWriter::new(File::create(&format!("{}.depthn", i)).unwrap());
            for d in depths {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    d.chromosome, d.basenumber, d.reads, d.name
                )
                .unwrap();
            }
            println!("Wrote {} to file {}.depthn", i, i);
        }
    } else {
        // Only one depth file to look at and write/print, stdout or outputfile if given
        let mut depths = read_depths(depthfiles[0]);
        add_name_to_depth(&mut depths, &bed_regions);
        if matches.is_present("output") {
            let o = matches.value_of("output").unwrap();
            let mut writer = BufWriter::new(File::create(&format!("{}.depthn", o)).unwrap());
            for i in depths {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    i.chromosome, i.basenumber, i.reads, i.name
                )
                .unwrap();
            }
        } else {
            println!("No output file defined, writing to stdout");
            for i in depths {
                println!(
                    "{}\t{}\t{}\t{}",
                    i.chromosome, i.basenumber, i.reads, i.name
                );
            }
        }
    }
}

fn write(depths: Vec<DepthInfo>, filename: &str, stdout: bool) {
    // Writes the depth file along with the name column
    let mut output = String::new();
    for depth in &depths {
        let _ = writeln!(
            output,
            "{}\t{}\t{}\t{}\n",
            depth.chromosome, depth.basenumber, depth.reads, depth.name
        );
    }
    if stdout {
        println!("{}", output);
    } else {
        let mut file = File::create(filename).unwrap();
        file.write_all(output.as_bytes()).unwrap();
    }
}

fn add_name_to_depth(depths: &mut Vec<DepthInfo>, bed_regions: &Vec<BedRegion>) {
    // Adds the name of the region to the depth file, based on the bed file

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
    // for i in 0..depths.len() {
    //     if depths[i].basenumber < bed_regions[idx].end {
    //         depths[i].name = bed_regions[idx].name.clone();
    //     } else {
    //         if idx < bed_regions.len() - 1 {
    //             idx += 1;
    //         }
    //         depths[i].name = bed_regions[idx].name.clone();
    //     }
    // }
}

/// Read a .bed file that contains choromome, start, end and name of region
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
    depths
}
