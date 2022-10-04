use csv::Reader;
use rayon::prelude::*;
use std::{env, error::Error, fs::File, io};

type BRecord = (String, u64, u64, String);
type DRecord = (String, u64, u32);

// enum Recs {
//     BRecord,
//     DRecord,
// }

fn main() {
    let args: Vec<String> = env::args().collect();
    let d_file = &args[1];
    let b_file = &args[2];
    add_name(d_file.to_string(), b_file.to_string()).unwrap();
}

fn add_name(d_file: String, b_file: String) -> Result<(), Box<dyn Error>> {
    let d_file = File::open(d_file)?;
    let b_file = File::open(b_file)?;
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());
    let mut d_rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(d_file);
    let mut b_rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(b_file);

    // parallelise the code below using rayon

    for entry in b_rdr.deserialize() {
        let bed_entry: BRecord = entry?;
        let end = bed_entry.2;
        for result in d_rdr.deserialize() {
            let record: DRecord = result?;
            let basepos = record.1;
            if basepos > end {
                break;
            }
            wtr.serialize((&record.0, &record.1, &record.2, &bed_entry.3))?;
        }
    }
    wtr.flush()?;
    Ok(())
}
