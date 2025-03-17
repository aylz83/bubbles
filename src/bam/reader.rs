#![allow(dead_code)]

use tokio::fs::File as TokioFile;
use tokio::io::BufReader as TokioBufReader;

use std::path::Path;

use crate::error;

pub struct Reader
{
	reader: TokioBufReader<TokioFile>,
}

impl Reader
{
	pub async fn from_path(path: &Path) -> error::Result<Reader>
	{
		let file = TokioFile::open(path)
			.await
			.map_err(|_| error::Error::IOError(path.to_string_lossy().to_string()))?;

		let reader = TokioBufReader::new(file);
		Ok(Reader { reader })
	}

	// 	pub async fn read_blocks<F>(
	// 		&mut self,
	// 		mut read_fn: F,
	// 	) -> error::Result<Vec<bam::pileup::Pileup>>
	// 	where
	// 		F: FnMut(&bam::Field),
	// 	{
	// 		let mut pileup = FxHashMap::<(i32, i32), u64>::default();

	// 		loop
	// 		{
	// 			match bam::read_bam_block(&mut self.reader, &mut pileup, &mut read_fn, None, &None)
	// 				.await?
	// 			{
	// 				Some(_bam_block) =>
	// 				{}
	// 				None =>
	// 				{
	// 					debug!("BAM EOF!");
	// 					break;
	// 				}
	// 			}
	// 		}

	// 		let mut pileup: Vec<_> = pileup
	// 			.into_iter()
	// 			.map(|((tid, pos), score)| bam::pileup::Pileup { tid, pos, score })
	// 			.collect();
	// 		pileup.sort_by(|a, b| (a.tid, a.pos).cmp(&(b.tid, b.pos)));

	// 		Ok(pileup)
	// 	}
}
