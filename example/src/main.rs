use std::io::Cursor;

use std::sync::Arc;

use bubbles::pufferfish::pool::BgzfBlockPool;
use bubbles::prelude::*;

use tracing_subscriber::EnvFilter;

#[tokio::main]
async fn main() -> anyhow::Result<()>
{
	tracing_subscriber::fmt()
		.with_env_filter(EnvFilter::from_default_env())
		.init();

	let mut read_count = 0;

	let data = std::fs::read("sample.bam")?;

	let block_pool = Arc::new(BgzfBlockPool::new(1000, 64 * 1024));
	let cursor = Cursor::new(data);
	let mut bam = bubbles::bam::AsyncBamQuery::from_reader(cursor, None, block_pool)
		.await
		.expect("not bam format");

	let pileup = bam
		// add BamFeatures::PILEUP to return pileup data
		.set_features(BamFeatures::CIGAR | BamFeatures::READNAMES)
		.query_reads(|read, _| {
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

	Ok(())
}
