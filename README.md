# nimblefin
A (terribly named) fast (and barebone!) BAM parser and pileup rust crate

This is a rust-based BAM parser with pileup capabilities and aimed to be as fast as possible while using as little memory as possible.
Full pileup approximately takes ~90 seconds on a 1.7GB BAM (102,650,976 reads)
It is/was a learning exercise and it is absoluely missing some features and may not handle all BAM files.

# Example
See src/main.rs for example usage

# Project aims (Reasons for existing?)
- To be as fast as possible
- To use as little memory as possible
- To make use of tokio/async
- To be a learning exercise
