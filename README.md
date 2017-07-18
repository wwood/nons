`nons` is a simple and hopefully not too slow program to fill Ns in contigs, by taking the consensus of reads aligned to these contigs.

To run it, take your sorted BAM file (perhaps created with BamM) and run nons:
```
nons contigs_with_Ns.fna reads_aligned_to_contigs_with_Ns.bam >filled_contigs.fna
```
`filled_contigs.fna` should now contain no `N`s, unless there was no Ns aligned at those positions.

## Install
After installing [rust](http://rust-lang.org/), install nons through `cargo`
```
cargo install nons
```

## Help
Probably the easiest would be to raise an [issue on GitHub](https://github.com/wwood/nons/issues) or contact the author [directly](http://ecogenomic.org/personnel/dr-ben-woodcroft).

## License
nons is written by [Ben Woodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics (UQ)](http://ecogenomic.org/) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html). A copy of the LICENSE is included in the repository.
