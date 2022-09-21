use clap::{App, Arg};
use std::env;
use std::fs;

#[derive(Debug)]
#[warn(dead_code)]
struct bed_region {
    chromosome: String,
    start: u32,
    end: u32,
    name: String,
}
#[derive(Debug)]
struct depth_info {
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
            "Adds names defined in a bedfile to the output of \"samtools depth\" command output, as the resulting depth file does not contain the name of what genic region the base belongs to",
        )
        .arg(
            Arg::with_name("depth")
                .short('d')
                .long("depth")
                .value_name("depth")
                .help("Depth file from \"samtools depth\" command to add names to")
                .takes_value(true)
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
    let depthfile = matches.value_of("depth").unwrap();
    let mut tRNA_regions: Vec<bed_region> = read_bed(bedfile);
    let mut depths: Vec<depth_info> = read_depths(depthfile);
    // Check each tRNA depth basenumber against each tRNA region start and end
    // If the basenumber is within the start and end, add the name of the tRNA to the tRNA_depths struct
    for i in 0..depths.len() {
        for j in 0..tRNA_regions.len() {
            if depths[i].basenumber >= tRNA_regions[j].start
                && depths[i].basenumber <= tRNA_regions[j].end
            {
                depths[i].name = tRNA_regions[j].name.clone();
            }
        }
    }
    depths.sort_by(|a, b| a.basenumber.cmp(&b.basenumber));

    if matches.is_present("output") {
        let outputfile = matches.value_of("output").unwrap();
        let mut file = fs::File::create(outputfile).expect("Unable to create file");
        write_with_name(depths, outputfile)
    } else {
        // let mut stdout = String::from("");
        for i in 0..depths.len() {
            print!(
                "{}\t{}\t{}\t{}\n",
                depths[i].chromosome, depths[i].basenumber, depths[i].reads, depths[i].name,
            );
        }
    }
}

fn write_with_name(depths: Vec<depth_info>, filename: &str) {
    let mut output = String::new();
    for i in 0..depths.len() {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            depths[i].chromosome, depths[i].basenumber, depths[i].reads, depths[i].name,
        ));
    }
    fs::write(filename, output).expect("Unable to write file");
}

fn read_bed(bedfile: &str) -> Vec<bed_region> {
    let content = fs::read_to_string(bedfile).expect("Something went wrong reading the file");
    let mut bed_regions: Vec<bed_region> = Vec::new();
    for line in content.lines() {
        let mut split_line = line.split("\t");
        let chromosome = split_line.next().unwrap().to_string();
        let start = split_line.next().unwrap().parse::<u32>().unwrap();
        let end = split_line.next().unwrap().parse::<u32>().unwrap();
        let name = split_line.next().unwrap().to_string();
        bed_regions.push(bed_region {
            chromosome,
            start,
            end,
            name,
        });
    }
    bed_regions
}

fn read_depths(filename: &str) -> Vec<depth_info> {
    let content = fs::read_to_string(filename).expect("Something went wrong reading the file");
    let mut depths: Vec<depth_info> = Vec::new();
    for line in content.lines() {
        let mut split_line = line.split("\t");
        let chromosome = split_line.next().unwrap();
        let position = split_line.next().unwrap().parse::<u32>().unwrap();
        let depth = split_line.next().unwrap().parse::<u32>().unwrap();
        depths.push(depth_info {
            chromosome: chromosome.to_string(),
            basenumber: position,
            reads: depth,
            name: "".to_string(),
        });
    }
    depths
}
