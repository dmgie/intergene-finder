#![allow(unused)]
use clap::{App, Arg};
use std::fs::File;
use std::io::{self, prelude::*, BufReader, BufWriter};

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
struct BedRegion {
    // chromosome: String,
    // start: u32,
    end: u32,
    name: String,
}
#[derive(Debug)]
struct DepthInfo {
    chromosome: String,
    basenumber: u32,
    reads: u32,
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
            add_name_to_depth(&mut depths, &bed_regions);
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
    for i in 0..depths.len() {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            depths[i].chromosome, depths[i].basenumber, depths[i].reads, depths[i].name,
        ));
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
    for i in 0..depths.len() {
        if depths[i].basenumber < bed_regions[idx].end {
            depths[i].name = bed_regions[idx].name.clone();
        } else {
            if idx < bed_regions.len() - 1 {
                idx += 1;
            }
            depths[i].name = bed_regions[idx].name.clone();
        }
    }
}

fn read_bed(bedfile: &str) -> Vec<BedRegion> {
    let content = File::open(bedfile).expect("Unable to open file");
    let reader = BufReader::new(content);
    let mut bed_regions: Vec<BedRegion> = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split_line = line.split("\t");
        let chromosome = split_line.next().unwrap().to_string();
        let start = split_line.next().unwrap().parse::<u32>().unwrap();
        let end = split_line.next().unwrap().parse::<u32>().unwrap();
        let name = split_line.next().unwrap().to_string();
        bed_regions.push(BedRegion {
            // chromosome,
            // start,
            end,
            name,
        });
    }
    bed_regions
}

fn read_depths(filename: &str) -> Vec<DepthInfo> {
    let content = File::open(filename).expect("Unable to open file");
    let mut depths: Vec<DepthInfo> = Vec::new();
    let reader = BufReader::new(content);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split_line = line.split("\t");
        let chromosome = split_line.next().unwrap();
        let position = split_line.next().unwrap().parse::<u32>().unwrap();
        let depth = split_line.next().unwrap().parse::<u32>().unwrap();
        depths.push(DepthInfo {
            chromosome: chromosome.to_string(),
            basenumber: position,
            reads: depth,
            name: "".to_string(),
        });
    }
    depths
}
