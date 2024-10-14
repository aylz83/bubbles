use nimblefin::bam::reader::Reader;

use std::path::Path;

#[tokio::main]
async fn main() -> anyhow::Result<()>
{
	env_logger::init();

	let mut reader = Reader::from_path(Path::new("sample.bam")).await?;

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
		.await?;

	println!("read count = {}", read_count);

	println!("top 10 pileup:");
	let to_break = 10;
	for (index, pileup) in pileup_vec.iter().enumerate()
	{
		if index == to_break
		{
			break;
		}

		println!("tid {}:pos {} = {}", pileup.tid, pileup.pos, pileup.score);
	}

	Ok(())
}
