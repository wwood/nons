extern crate bio;
use bio::io::fasta;

use std::collections::HashMap;
use std::str;

// read a fasta file into a HashMap of sequence name to sequence.
#[allow(dead_code)]
pub fn read_fasta(reference_fasta: &str) -> HashMap<String,String>{
    let reader = fasta::Reader::from_file(reference_fasta).unwrap();
    let mut contig_to_seq = HashMap::new();
    for record in reader.records() {
        let rec = record.unwrap();
        //println!("{} {}", rec.id().unwrap(), str::from_utf8(rec.seq()).unwrap());
        // Using String::from here is a cludge, but I cannot work out how to copy a str
        contig_to_seq.insert(
            String::from(rec.id()),
            String::from(str::from_utf8(rec.seq()).unwrap())
        );
    }
    return contig_to_seq;
}


// read a fasta file into a HashMap of sequence name to sequence.
#[allow(dead_code)]
pub fn read_fasta2(reference_fasta: &str) -> HashMap<String,Vec<u8>>{
    let reader = fasta::Reader::from_file(reference_fasta).unwrap();
    let mut contig_to_seq = HashMap::new();
    for record in reader.records() {
        let rec = record.unwrap();
        //println!("{} {}", rec.id().unwrap(), str::from_utf8(rec.seq()).unwrap());
        // Using String::from here is a cludge, but I cannot work out how to copy a str
        let my_seq = rec.seq().to_owned();
        contig_to_seq.insert(
            String::from(rec.id()),
            my_seq
        );
    }
    return contig_to_seq;
}
