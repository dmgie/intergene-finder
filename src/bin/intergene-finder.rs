// Read file
use clap::{App, Arg};
use std::fs;
// use std::fs::File;
// use std::io::{self, prelude::*, BufReader, BufWriter};

#[derive(Debug)]
enum GFFErrors {
    FileNotFound,
    // FileParseError,
    // InvalidExtension,
    // InvalidGFFEntry,
    GFFColumnMismatch,
}

#[derive(Debug)]
enum FastaErrors {
    FileNotFound,
}

#[allow(dead_code)]
#[derive(Debug, PartialEq, Eq, Ord, PartialOrd, Clone)]
struct GffEntry {
    seqid: String,
    source: String,
    r#type: String,
    start: i32,
    end: i32,
    score: String,
    strand: String,
    phase: String,
    attributes: String,
    seq: String,
}

impl GffEntry {
    fn add_seq(&mut self, seq: String) {
        self.seq = seq;
    }
}

struct GFF {
    header: String,
    entries: Vec<GffEntry>,
}

#[allow(dead_code)]
struct Seq {
    header: String,
    seq: String,
}

// TODO: If given a fasta, also extract sequences from the fasta, if not don't, just create a .gff file including the intergenic regions
// TODO: allow for adding a padding/buffer option in which there either has to be a min length for IGR's, or the distance between the IGR and the previous/next
// gene has to at least be that minimum distance apart i.e GENE ---40nt--> IGR ---40nt---> Gene

fn main() {
    let matches = App::new("Finds intergenic regions, creates FASTA and GFF file")
        .version("0.1")
        .author("Daniel Giesel")
        .about(
            "Extracts intergenic region from a gff file, creating a new gff file and (optionally) a fasta file with the sequences",
        )
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("input")
                .help("Input GFF file to find the intergenic regions from")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("fasta")
                .short('f')
                .long("fasta")
                .value_name("fasta")
                .help("The genome sequence as a fasta file from which to extract sequence from")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("types")
                .short('t')
                .long("types")
                .value_name("types")
                .help("The entry types to extract sequences from, separated by commas")
                .takes_value(true)
                .multiple(true)
                .value_delimiter(',')
        )
        .arg(
            Arg::with_name("min_distance")
                .short('d')
                .long("min_distance")
                .value_name("min_distance")
                .help("The minimum distance between two genes to be considered an intergenic region")
                .default_value("3")
                .default_missing_value("3")
                .takes_value(true)
        )
        .get_matches();

    let refgff = parse_gff(matches.value_of("input").expect("Expect input file"))
        .expect("Error parsing GFF file");

    // TODO: Change it somehow to look at if the fasta file is given, if not, don't extract sequences
    let reffasta =
        parse_fasta(matches.value_of("fasta").unwrap()).expect("Error parsing FASTA file");

    let refseq = &reffasta[0].seq; // Get the first sequence from the fasta file i.e the only
    let refseqlen = refseq.len() as i32;

    let gff_entries = refgff.entries;

    // Get intergenic regions (start, end) in a vector
    let min_distance = matches
        .value_of("min_distance")
        .unwrap()
        .parse::<i32>()
        .unwrap();
    let intergenic_regions: Vec<(i32, i32)> =
        get_intergenic_regions(&gff_entries, refseqlen, min_distance);
    let intergenic_entries =
        create_intergenic_entries(intergenic_regions, gff_entries[0].seqid.clone());

    // Merge intergenic and gff entries and sort them
    let mut merged_entries: Vec<GffEntry> = gff_entries.to_vec();
    merged_entries.extend(intergenic_entries);
    merged_entries.sort_by(|a, b| a.start.cmp(&b.start));

    // Extract sequences from fasta file and add them to the entries
    let merged_with_seq = add_seq_to_entries(&mut merged_entries, refseq);

    // Write new gff file, fasta files
    write_gff_from_vec(&refgff.header, &merged_with_seq, "reference+intergenic.gff");

    // Write fasta files depending on the types given
    // go through the entries and collect all the type fields
    let mut valid_types: Vec<String> = Vec::new();
    for entry in &merged_with_seq {
        if !valid_types.contains(&entry.r#type) {
            valid_types.push(entry.r#type.clone());
        }
    }
    // Only write fasta files for valid types
    for entry_type in matches.values_of("types").unwrap() {
        if valid_types.contains(&entry_type.to_string()) {
            write_fasta_to_file(
                entry_type,
                &merged_with_seq,
                &format!("{}.fasta", entry_type),
            );
        } else {
            println!(
                "\x1b[91mERROR: Invalid entry type:\x1b[0m \x1b[94m{t}\x1b[0m. Not creating a fasta file for {t}. Please check if it was spelled correctly.",
                t = entry_type
            );
        }
    }
}

/// For each intergenic region create a GFFEntry that has the same format as the other entries
/// mainly start, end, type, and attributes i.e (ID,Name,locus_tag) defined
fn create_intergenic_entries(intergenic_regions: Vec<(i32, i32)>, seqid: String) -> Vec<GffEntry> {
    let mut igr_counter = 0;
    let intergenic_entries: Vec<GffEntry> = intergenic_regions
        .iter()
        .map(|(start, end)| {
            igr_counter += 1;
            // For getting the seq for all intergenic regions
            let seq = String::from("");
            // let seq = refseq
            //     .get(*start as usize..*end as usize)
            //     .unwrap()
            //     .to_string();
            GffEntry {
                // same seqid as the rest
                seqid: seqid.clone(),
                source: "intergene-finder".to_string(),
                r#type: "intergenic".to_string(),
                start: *start,
                end: *end,
                score: ".".to_string(),
                strand: ".".to_string(),
                phase: ".".to_string(),
                attributes: format!(
                    "ID=IGR-{};Name=INTERGENIC_{};locus_tag=INTERGENIC_{}",
                    igr_counter, igr_counter, igr_counter
                )
                .to_string(),
                seq,
            }
        })
        .collect();
    intergenic_entries
}

// Given a vector of GFF entries, add the sequence to each entry
fn add_seq_to_entries(entries: &mut Vec<GffEntry>, refseq: &String) -> Vec<GffEntry> {
    let mut entries_with_seq: Vec<GffEntry> = Vec::new();
    for entry in entries {
        let seq = refseq
            .get(entry.start as usize - 1..entry.end as usize)
            .unwrap()
            .to_string();
        entry.add_seq(seq);
        entries_with_seq.push(entry.clone());
    }
    entries_with_seq
}

#[allow(dead_code)]
fn parse_fasta(file: &str) -> Result<Vec<Seq>, FastaErrors> {
    let file = match fs::read_to_string(file) {
        Ok(file) => file,
        Err(_) => return Err(FastaErrors::FileNotFound),
    };
    let mut entries: Vec<Seq> = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in file.lines() {
        let line = line.to_string();
        if line.starts_with(">") {
            if !header.is_empty() {
                entries.push(Seq { header, seq });
            }
            header = line;
            seq = String::new();
        } else {
            seq.push_str(&line);
        }
    }
    entries.push(Seq { header, seq });
    Ok(entries)
}

fn parse_gff(file: &str) -> Result<GFF, GFFErrors> {
    // Take a filename and parse it into a GFF struct

    let file = match fs::read_to_string(file) {
        Ok(file) => file,
        Err(_) => return Err(GFFErrors::FileNotFound),
    };

    // Header info
    let mut header = String::new();
    for line in file.lines() {
        let line = line.to_string();
        if line.starts_with("#") {
            header.push_str(&format!("{}\n", &line));
        } else {
            break;
        }
    }

    // Entry info
    let mut entries: Vec<GffEntry> = Vec::new();
    for line in file.lines() {
        let line = line.to_string();
        // Skip comments & header lines
        if line.starts_with("#") {
            continue;
        }

        let parts: Vec<&str> = line.split("\t").collect();

        // Check if the line has the correct number of columns, create entry if so
        // error if not
        let gffrecord = match parts.len() {
            9 => GffEntry {
                seqid: parts[0].to_string(),
                source: parts[1].to_string(),
                r#type: parts[2].to_string(),
                start: parts[3].parse::<i32>().unwrap(),
                end: parts[4].parse::<i32>().unwrap(),
                score: parts[5].to_string(),
                strand: parts[6].to_string(),
                phase: parts[7].to_string(),
                attributes: parts[8].to_string(),
                seq: String::new(),
            },
            _ => return Err(GFFErrors::GFFColumnMismatch),
        };
        entries.push(gffrecord);
    }
    Ok(GFF { header, entries })
}

/// Given a (vector of) GFF entry struct(s) return a the start and end of the intergenic regions
/// as a vector of tuples
fn get_intergenic_regions(gff: &Vec<GffEntry>, end: i32, buffer: i32) -> Vec<(i32, i32)> {
    // We obtain all the intergenic regions by going through a vector of GFFEntries
    let mut regions: Vec<(i32, i32)> = Vec::new();

    // should it start with 0 or 1? 0 means it includes (0,1) as intergenic range, if 1 is the first pos
    // if 1 is the first
    let mut last_end = match gff.first().unwrap().start {
        1 => 1,
        _ => 0,
    };

    for entry in gff {
        // skip if type is any of:
        if entry.r#type == "region"
            || entry.r#type == "sequence_feature"
            || entry.r#type == "start_codon"
            || entry.r#type == "stop_codon"
        {
            continue;
        }
        // If it happens right after stop codon, automatically update last end to be + 1
        if entry.start == last_end + 1 {
            last_end = last_end + 1;
        }
        // if next entry is less than (i.e within) Xnt of previous end/entry, do not mark it as an intergenic
        // if entry.start < last_end + 40 {
        //     continue;
        // }

        if entry.start > last_end + buffer {
            regions.push((last_end + 1, entry.start - 1)); // NOTE: Remove if include end of previous seq
                                                           // regions.push((last_end + 1, entry.start - 1));
        }
        last_end = entry.end;
    }
    // Fill in the last if the end of the gff hasnt reached the end refernece sequence
    if last_end < end {
        regions.push((last_end, end)); //NOTE: Remove if include end of previous seq
                                       // regions.push((last_end + 1, end));
    }
    regions
}

/// More generic fasta writer
fn write_fasta_to_file(entry_type: &str, gff_entries: &Vec<GffEntry>, filename: &str) {
    let mut to_write = String::new();
    // For entries matching type, write their sequences to a file
    for entry in gff_entries {
        if entry.r#type == entry_type {
            to_write.push_str(&format!(
                ">{} length: {}\n",
                entry.attributes,
                entry.seq.len()
            ));
            let mut seq = entry.seq.clone();
            while seq.len() > 80 {
                to_write.push_str(&format!("{}\n", &seq[..80]));
                seq = seq[80..].to_string();
            }
            to_write.push_str(&format!("{}\n", &seq));
        }
        // to_write.push_str(&format!(">{}\n{}\n", entry.attributes, entry.seq));
    }
    fs::write(filename, to_write).expect("Unable to write file");
}

// Create a new GFF file that includes the (intergenic) regions that we added
fn write_gff_from_vec(header: &String, gff_entries: &Vec<GffEntry>, fname: &str) {
    // Recreate a gff file from the header, entries to include and write it to a file
    let mut gff_file = String::new();
    gff_file.push_str(&format!("{}\n", header));
    for entry in gff_entries {
        gff_file.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            entry.seqid,
            entry.source,
            entry.r#type,
            entry.start,
            entry.end,
            entry.score,
            entry.strand,
            entry.phase,
            entry.attributes,
        ));
    }
    fs::write(fname, gff_file).unwrap();
}
