use std::sync::Arc;
use std::env;

use std::time::Instant;

use bubbles::bam::FetchRegion;

use bubbles::pufferfish::pool::BgzfBlockPool;
use bubbles::prelude::*;

use tracing_subscriber::EnvFilter;

#[tokio::main]
async fn main() -> anyhow::Result<()>
{
	tracing_subscriber::fmt()
		.with_env_filter(EnvFilter::from_default_env())
		.init();

	let args: Vec<String> = env::args().collect();

	let mut read_count = 0;

	let bam_file = &args[1];

	println!("Input file: {}", &bam_file);

	let block_pool = Arc::new(BgzfBlockPool::new(1000, 64 * 1024));

	let start = Instant::now();

	let mut bam = bubbles::bam::AsyncBamQuery::from_path(bam_file, None, block_pool)
		.await
		.expect("not bam format");

	let pileup = bam
		// add BamFeatures::PILEUP to return pileup data
		.set_features(BamFeatures::CIGAR | BamFeatures::READNAMES)
		.oneshot_query(&vec![FetchRegion::from("chr9")], |read, _| {
			read_count += 1;

			/*let len = read.sequence.len();
			read_lengths.push(len);

			let ref_id = read.ref_id;
			if let Some(rid) = header.ref_name(ref_id)
			{
				*ref_counts.entry(rid.name.clone()).or_insert(0) += 1;
			}*/

			Some(read)
		})
		.await
		.expect("unable to read sequences");

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

	let elapsed = start.elapsed();
	println!("Reader creation + processing took {:?}", elapsed);

	Ok(())
}
