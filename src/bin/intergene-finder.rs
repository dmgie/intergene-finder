// Read file
use clap::{App, Arg};
use std::fs;
use std::io::Error;
// use std::fs::File;
// use std::io::{self, prelude::*, BufReader, BufWriter};
//
// TODO: Strand-specific intergene finding. This is because the start/stop are switchedin location dependingon the strand

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
    start: i64,
    end: i64,
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
            Arg::with_name("strandedness")
                .short('s')
                .long("strandedness")
                .value_name("strandedness")
                .help("Flag whether to keep strandedness in mind or not when calculating intergenic regions")
                .required(false)
                .takes_value(false)
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
        .arg(
            Arg::with_name("output")
                .short('o')
                .long("output")
                .value_name("output")
                .help("The output folder, if not given, will be in the same folder as the input file")
                .takes_value(true)
                .required(false)
                .default_value(".")
        )
        .get_matches();

    let refgff = parse_gff(
        matches
            .get_one::<String>("input")
            .expect("Expect input file"),
    )
    .expect("Error parsing GFF file");
    // TODO: Change it somehow to look at if the fasta file is given, if not, don't extract sequences
    let reffasta =
        parse_fasta(matches.get_one::<String>("fasta").unwrap()).expect("Error parsing FASTA file");
    let refseq = &reffasta[0].seq; // Get the first sequence from the fasta file i.e the only
    let refseqlen = refseq.len() as i64;

    let gff_entries = refgff.entries;

    // Get intergenic regions (start, end) in a vector
    let min_distance = matches
        .get_one::<String>("min_distance")
        .unwrap()
        .parse::<i64>()
        .unwrap();

    // Add intergenic entries for each strand separately to the GFF entries
    let mut merged_entries: Vec<GffEntry> = gff_entries.to_vec();
    if matches.is_present("strandedness") {
        for strand in ["+", "-"] {
            let intergenic_region: Vec<(i64, i64)> =
                get_intergenic_regions(&gff_entries, refseqlen, min_distance, strand);
            let intergenic_entries =
                create_intergenic_entries(intergenic_region, gff_entries[0].seqid.clone(), strand);
            merged_entries.extend(intergenic_entries);
        }
    } else {
        let intergenic_region: Vec<(i64, i64)> =
            get_intergenic_regions(&gff_entries, refseqlen, min_distance, ".");
        let intergenic_entries =
            create_intergenic_entries(intergenic_region, gff_entries[0].seqid.clone(), ".");
        merged_entries.extend(intergenic_entries);
    }

    // Merge intergenic and gff entries and sort them
    merged_entries.sort_by(|a, b| a.end.cmp(&b.end));

    // Extract sequences from fasta file and add them to the entries
    let merged_with_seq = add_seq_to_entries(&mut merged_entries, refseq);

    // Write new gff file, fasta files
    write_gff_from_vec(&refgff.header, &merged_with_seq, "reference+intergenic.gff")
        .expect("Unable to create the GFF file");

    // Write fasta files depending on the types given
    // go through the entries and collect all the type fields
    let mut valid_types: Vec<String> = Vec::new();
    for entry in &merged_with_seq {
        if !valid_types.contains(&entry.r#type) {
            valid_types.push(entry.r#type.clone());
        }
    }

    // Only write fasta files for valid types
    for entry_type in matches.get_many::<String>("types").expect("No types given") {
        if valid_types.contains(&entry_type.to_string()) {
            let output = matches
                .get_one::<String>("output")
                .expect("No output path given");
            write_fasta_to_file(
                entry_type,
                &merged_with_seq,
                &format!("{out}/{ttype}.fasta", out = output, ttype = entry_type),
            )
            .expect("Unable to write fasta file(s)");
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
fn create_intergenic_entries(
    intergenic_regions: Vec<(i64, i64)>,
    seqid: String,
    strand: &str,
) -> Vec<GffEntry> {
    let mut igr_counter = 0;
    let intergenic_entries: Vec<GffEntry> = intergenic_regions
        .iter()
        .map(|(start, end)| {
            igr_counter += 1;
            let seq = String::from("");
            GffEntry {
                // same seqid as the rest
                seqid: seqid.clone(),
                source: "intergene-finder".to_string(),
                r#type: "intergenic".to_string(),
                start: *start,
                end: *end,
                score: ".".to_string(),
                strand: strand.to_string(),
                phase: ".".to_string(),
                attributes: format!(
                    "ID=IGR_{a}({strand});Name=INTERGENIC_{a}({strand});locus_tag=INTERGENIC_{a}({strand})",
                    a = igr_counter,
                    strand = strand
                )
                    .to_string(),
                seq,
            }
        })
        .collect();
    intergenic_entries
}

/// Given a (vector of) GFF entry struct(s) return a the start and end of the intergenic regions
/// as a vector of tuples
fn get_intergenic_regions(
    gff: &Vec<GffEntry>,
    end: i64,
    buffer: i64,
    strand: &str,
) -> Vec<(i64, i64)> {
    // We obtain all the intergenic regions by going through a vector of GFFEntries
    let mut regions: Vec<(i64, i64)> = Vec::new();

    // should it start with 0 or 1? 0 means it includes (0,1) as intergenic range, if 1 is the first pos
    // if 1 is the first
    let mut last_end = match gff.first().unwrap().start {
        1 => 1,
        _ => 0,
    };

    // Filter if stranded and given, otherwise just go through all the entries i.e non-strand-specific IGRs
    let filtered = match strand {
        "+" => gff
            .iter()
            .filter(|x| x.strand == "+" || x.strand == ".")
            .collect::<Vec<&GffEntry>>(),
        "-" => gff
            .iter()
            .filter(|x| x.strand == "-" || x.strand == ".")
            .collect::<Vec<&GffEntry>>(),
        _ => gff.iter().collect::<Vec<&GffEntry>>(),
    };

    for entry in filtered {
        // skip if type is any of:
        if entry.r#type == "region" || entry.r#type == "sequence_feature" {
            continue;
        }
        // If it happens right after stop codon, automatically update last end to be + 1
        if entry.start == last_end + 1 {
            last_end = last_end + 1;
        }
        // if next entry is less than (i.e within) Xnt of previous end/entry, do not mark it as an intergenic
        if entry.start < last_end + buffer {
            continue;
        }

        // if entry.start > last_end + buffer {
        if entry.start > last_end {
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
                start: parts[3].parse::<i64>().unwrap(),
                end: parts[4].parse::<i64>().unwrap(),
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
/// More generic fasta writer
fn write_fasta_to_file(
    entry_type: &str,
    gff_entries: &Vec<GffEntry>,
    filename: &str,
) -> Result<(), Error> {
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
    fs::write(filename, to_write)?;
    Ok(())
}

// Create a new GFF file that includes the (intergenic) regions that we added
fn write_gff_from_vec(
    header: &String,
    gff_entries: &Vec<GffEntry>,
    fname: &str,
) -> Result<(), Error> {
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
    fs::write(fname, gff_file)?;
    Ok(())
}
