# bubbles
A terribly named, "fast" BAM parser and pileup rust crate making use of SIMD

Originally, this crate started out life as a just a learning curve to better understand the BAM specification and learn SIMD.

It has since evolved, supporting the following:

- BAI for reading just portions of the BAM
- Requesting just specific types of data from the BAM, for example, just the read names and sequences, or tags. This can significantly
  increase parsing speed if you're only interested in just portions of the data.
- Reading BAMs from a tokio BufReader, enabling parsing of BAMs directly from memory, a file, etc
- A unique, conditional-based pileup system, as reads are provided within the fetch_reads, function, if you don't return the read and
  instead return None, the read won't be included the final pileup scores.

# Example

```rust
use bubbles::*;

#[tokio::main]
async fn main()
{
  let mut read_count = 0;

  // Automatically searches for a sample.bam.bai, but this can optionall be supplied with set_bai if the path/name differs
	let pileup = bam::Builder::from_path("sample.bam")
	.await
  .expect("Unable to open BAM")
  // Fetch all of chromosome 3
  // Only available when a BAI is present, omitting add_fetch_regions reads the entire BAM
  // NamedTid also supports transcript_ids if the BAM is aligned to the transcriptome
	.add_fetch_region(bam::FetchRegion::NamedTid("chr3"))
  .expect("Unable to add search region")
  // Fetch a section of chromosome 2
	.add_fetch_region(bam::FetchRegion::NamedTidRegion(
	  "chr2", 119998389, 201064349,
	)
  .expect("Unable to add search region")
  // what features to request, omitting this enables all BAM featurs
  .set_bam_features(bubbles::bam::BamFeatures::CIGAR | bubbles::bam::BamFeatures::READNAMES)
  .fetch_reads(
		|read, header|
    {
			read_count += 1;

			if let Some(rid) = header.ref_name(read.ref_id)
			{
        println!("TID name: {:?}", rid);
			}

      // Only pileup reads if read names contain a string
      if read.read_name.contains("some imaginary requirement")
      {
        return Some(read);
      }
  
			None
	  },
    // Request CIGAR ops and read names.
    // Setting this to None will enable all features to be requested
	  Some(bubbles::bam::BamFeatures::CIGAR | bubbles::bam::BamFeatures::READNAMES),
	)
	.await
	.expect("unable to read sequences");

  // Alternatively, a builder can be opened as follows:
  // let bam_buffer = ... Some aysnc BufReader pointing to memory or a stream
  // let bai_buffer = ... Some aysnc BufReader pointing to memory or a stream
  // let builder = bam::Builder::from_reader(bam_buffer, Some(bai_buffer)).await.expect("unable to open");
  // Or if you don't wish to suppoy BAI data
  // let builder = bam::Builder::from_reader(bam_buffer, None).await.expect("unable to open");
  // Do something with builder (set_features, add_fetch_region, fetch_reads etc)...

	println!("read count = {}", read_count);


  if let Some(pileup) = pileup
  {
  	println!("first 10 pileup positions:");
  	let to_break = 10;
  	for (index, pileup) in pileup.iter().enumerate()
  	{
  		if index == to_break
  		{
  			break;
  		}

  		println!("tid {}:pos {} = {}", pileup.tid, pileup.pos, pileup.score);
  	}
  }

	Ok(())
}
```

# Project aims (Reasons for existing?)
- To be as fast as possible - some room for improvement here!
- To use as little memory as possible
- Support reading from other than a file
- Conditional pileup
- To make use of tokio/async
