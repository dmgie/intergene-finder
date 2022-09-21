// Read file
use std::env;
use std::fs;

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

#[allow(dead_code)]
fn parse_fasta(file: &str) -> Vec<Seq> {
    let mut seqs: Vec<Seq> = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in file.lines() {
        let line = line.to_string();
        if line.starts_with(">") {
            if !header.is_empty() {
                seqs.push(Seq { header, seq });
            }
            header = line;
            seq = String::new();
        } else {
            seq.push_str(&line);
        }
    }
    seqs.push(Seq { header, seq });
    seqs
}

fn parse_gff(file: &str) -> GFF {
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
        if line.starts_with("#") {
            continue;
        }
        let mut parts = line.split("\t");
        let seqid = parts.next().unwrap().to_string();
        let source = parts.next().unwrap().to_string();
        let r#type = parts.next().unwrap().to_string();
        // Ignore region parts
        // if r#type == "region" {
        //     continue;
        // }
        let start = parts.next().unwrap().to_string().parse::<i32>().unwrap();
        let end = parts.next().unwrap().to_string().parse::<i32>().unwrap();
        let score = parts.next().unwrap().to_string();
        let strand = parts.next().unwrap().to_string();
        let phase = parts.next().unwrap().to_string();
        let attributes = parts.next().unwrap().to_string();
        entries.push(GffEntry {
            seqid,
            source,
            r#type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
            seq: String::new(),
        });
    }

    GFF { header, entries }
}

fn get_intergenic_regions(gff: &Vec<GffEntry>, end: i32) -> Vec<(i32, i32)> {
    // We obtain all the intergenic regions by going through a vector of GFFEntries
    let mut regions: Vec<(i32, i32)> = Vec::new();
    // let mut last_end = 0; // should it be 0 or 1? 0 means it includes (0,1) as intergenic range
    let mut last_end = 1;
    for entry in gff {
        // skip if type is region, since that just is info about the entire genome
        if entry.r#type == "region" || entry.r#type == "sequence_feature" {
            continue;
        }
        if entry.start > last_end {
            regions.push((last_end, entry.start)); // NOTE: Remove if include end of previous seq
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

fn main() {
    let args: Vec<String> = env::args().collect();

    let refgff = &args[1];
    let reffasta = &parse_fasta(&fs::read_to_string(&args[2]).unwrap())[0];
    let refseq = &reffasta.seq;
    let refseqlen = refseq.len() as i32;

    let gff = parse_gff(&fs::read_to_string(refgff).unwrap());
    let mut gff_entries = gff.entries;

    // Get intergenic regions (start, end) in a vector
    let intergenic_regions = get_intergenic_regions(&gff_entries, refseqlen);
    let mut igr_counter = 0;
    // For each intergenic region create a GFFEntry that has the same format as the other entries
    // mainly start, end, type, and attributes i.e (ID,Name,locus_tag) defined
    let intergenic_entries: Vec<GffEntry> = intergenic_regions
        .iter()
        .map(|(start, end)| {
            igr_counter += 1;
            // For getting the seq for all intergenic regions
            // let seq = String::from("");
            let seq = refseq
                .get(*start as usize..*end as usize)
                .unwrap()
                .to_string();
            GffEntry {
                // same seqid as the rest
                seqid: reffasta.header.clone().split(" ").next().unwrap()[1..].to_string(),
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

    // Add sequence for GFF entry struc using the given reference seqence
    for entry in &mut gff_entries {
        let seq = &refseq[(entry.start) as usize..entry.end as usize];
        entry.add_seq(seq.to_string());
    }

    // Merge intergenic and gff entries and sort them
    let mut merged_entries: Vec<GffEntry> = gff_entries.to_vec();
    merged_entries.extend(intergenic_entries);
    merged_entries.sort_by(|a, b| a.start.cmp(&b.start));

    // Creat a GFF file from the merged entries
    write_gff_from_vec(&gff.header, &merged_entries, "reference+intergenic.gff");
    write_intergenic_to_file(&merged_entries, "intergenic.fasta")
}

// write intergenic regions as fasta
fn write_intergenic_to_file(gff_entries: &Vec<GffEntry>, filename: &str) {
    // Write the header and the sequence underneath with a maximum characters of 80
    let mut to_write = String::new();
    for entry in gff_entries {
        if entry.r#type == "intergenic" {
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

fn write_gff_from_vec(header: &String, gff_entries: &Vec<GffEntry>, fname: &str) {
    // Recreate a gff file from the header, merged entries and write it to a file
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
