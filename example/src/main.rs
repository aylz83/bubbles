use bubbles::bam;

use std::io::Cursor;

#[tokio::main]
async fn main() -> anyhow::Result<()>
{
	env_logger::init();

	let mut read_count = 0;

	// tid name in field
	// use tid in field and pileup
	// add star juncrion reading
	//
	let data = std::fs::read("sample.bam")?;

	let cursor = Cursor::new(data);
	let mut bam = bubbles::bam::Builder::from_reader(cursor, None)
		.await
		.expect("not bam format");

	let pileup = bam
		.fetch_reads(|read, header| {
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

	/*let pileup = bam::Builder::from_path("sample.bam")
	.await?
	//.add_fetch_region(bam::FetchRegion::NamedTid("chr2"))?
	//.add_fetch_region(bam::FetchRegion::NamedTidRegion(
	//	"chr2", 119998389, 201064349,
	//))?
	.fetch_reads(|mut read, header| {
		read_count += 1;
		println!("seqname = {:?}", read.ref_id);
		//println!("read seq = {:?}", unsafe {
		//	std::str::from_utf8_unchecked(read.sequence).to_string()
		//});
		//println!("read quality = {:?}", unsafe {
		//std::str::from_utf8_unchecked(read.sequence_quality).to_string()
		//});
		//
		Some(read)
	})
	.await?;*/

	//let mut bai = bai::Reader::from_path(Path::new("sample.bam.bai")).await?;

	/*let mut reader = bam::reader::Reader::from_path(Path::new("sample.bam")).await?;

	let _header = reader.read_header().await?;

	let mut read_count = 0;

	let pileup_vec = reader
		.read_blocks(|read| {
			read_count += 1;
			//println!("read name = {:?}", read.read_name);
			//println!("read seq = {:?}", unsafe {
			//	std::str::from_utf8_unchecked(read.sequence).to_string()
			//});
			//println!("read quality = {:?}", unsafe {
			//std::str::from_utf8_unchecked(read.sequence_quality).to_string()
			//});
		})
		.await?;*/

	println!("read count = {}", read_count);

	println!("top 10 pileup:");
	let to_break = 10;
	for (index, pileup) in pileup.iter().enumerate()
	{
		if index == to_break
		{
			break;
		}

		println!("tid {}:pos {} = {}", pileup.tid, pileup.pos, pileup.score);
	}

	Ok(())
}
