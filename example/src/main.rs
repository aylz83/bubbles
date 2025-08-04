use std::io::Cursor;

#[tokio::main]
async fn main() -> anyhow::Result<()>
{
	env_logger::init();

	let mut read_count = 0;

	let data = std::fs::read("sample.bam")?;

	let cursor = Cursor::new(data);
	let mut bam = bubbles::bam::Builder::from_reader(cursor, None)
		.await
		.expect("not bam format");

	let pileup = bam
		// add BamFeatures::PILEUP to return pileup data
		.set_features(bubbles::bam::BamFeatures::CIGAR | bubbles::bam::BamFeatures::READNAMES)
		.fetch_reads(|read, _| {
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
