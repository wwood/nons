extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::prelude::*;

extern crate bio;

use std::env;
use std::str;

mod util;

fn main() {
    let args: Vec<String> = env::args().collect();

    match args.len() {
        3 => huh(&args[1], &args[2], &mut std::io::stdout()),
        _ => print_help()
    }
}

fn print_help(){
    println!("Usage:");
    println!("");
    println!(" consensusm <reference_fasta> <aligned_reads.bam>");
    println!("");
}


fn huh<T: std::io::Write>(reference_fasta_path: &str, indexed_bam_file_path: &str, mut stream: &mut T){
    // println!("fasta: {}", reference_fasta_path);
    // println!("bam: {}", indexed_bam_file_path);
    let contig_to_seq = util::read_fasta(reference_fasta_path);
    let bam = bam::Reader::from_path(indexed_bam_file_path).unwrap();

    // pileup over all covered sites
    let mut last_tid: u32 = u32::max_value();
    let mut tid: u32;
    let mut ref_iterator = "".chars();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let mut counts: [u32; 256] = [0; 256];
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        tid = pileup.tid();
        if tid != last_tid {
            if last_tid != u32::max_value() {
                writeln!(&mut stream).unwrap();
            }

            let contig_name = str::from_utf8(
                bam.header().target_names()[pileup.tid() as usize]
            ).unwrap();
            writeln!(&mut stream, ">{}", contig_name).unwrap();
            last_tid = tid;
            ref_iterator = contig_to_seq[contig_name].chars();
        }
        let ref_base = ref_iterator.next().unwrap();
        if ref_base == 'N' {
            for alignment in pileup.alignments() {
                match alignment.qpos() {
                    Some(qpos) => {
                        let base: u8 = alignment.record().seq()[qpos];
                        counts[base as usize] += 1;
                    }
                    None => {} // ignore indels for the moment.
                }
            }
            {
                let mut max_count = 0;
                let mut max_base: char = '?';
                for (i, count) in counts.iter().enumerate() {
                    if *count > max_count {
                        max_count = *count;
                        max_base = char::from(i as u8);
                    }
                }
                write!(&mut stream, "{}",max_base).unwrap();
            }
        } else {
            write!(&mut stream, "{}",ref_base).unwrap();
        }
    }
    writeln!(&mut stream).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_no_n(){
        let mut stream = Cursor::new(Vec::new());
        huh("test-data/random500.fna",
            "test-data/random500.random500.reads.bam",
            &mut stream);
        assert_eq!(
            ">random_sequence_length_500
ACAACTGTCCCACAGAGGAATTTCGTCCTCCGCGTGTTTATGTCGATCCGGAGTACAAGAAGATTCTCACGTTGAAGATTGAATAGAAGTAGCAGTTTAGCGATATTGACAGTATTGAATACTTAAGCGCAACCCTCAAGATCTCTCCAGGAAACTATCGTAGAAAGGTCCATGTAGCGTAGCCCTAGTGCGTAACGCTCTGGCTATAACGCCTCGGCATGGTCGTGCGAGTAATACAAGCCAGACACAGAAAGAAACAAAGTAACTGGGATACAATTGCTCAAATCGCACCATCCCTTGCCTACCCGCTG
",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_n(){
        let mut stream = Cursor::new(Vec::new());
        huh("test-data/random500.one_N.fna",
            "test-data/random500.random500.reads.bam",
            &mut stream);
        assert_eq!(
            ">random_sequence_length_500
ACAACTGTCCCACAGAGGAATTTCGTCCTCCGCGTGTTTATGTCGATCCGGAGTACAAGAAGATTCTCACGTTGAAGATTGAATAGAAGTAGCAGTTTAGCGATATTGACAGTATTGAATACTTAAGCGCAACCCTCAAGATCTCTCCAGGAAACTATCGTAGAAAGGTCCATGTAGCGTAGCCCTAGTGCGTAACGCTCTGGCTATAACGCCTCGGCATGGTCGTGCGAGTAATACAAGCCAGACACAGAAAGAAACAAAGTAACTGGGATACAATTGCTCAAATCGCACCATCCCTTGCCTACCCGCTG
",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
