use std::path::Path;
use std::ops::Range;

use std::io::{Seek, Cursor};
use std::collections::BTreeSet;

use tokio::fs::File as TokioFile;
use tokio::io::{SeekFrom, AsyncSeekExt, BufReader as TokioBufReader};

use log::debug;

use crate::error;
use crate::bai;
use crate::bam;
use crate::AsyncReadSeek;

use pufferfish::BGZ;

pub enum FetchRegion<'a>
{
	Tid(u32),
	NamedTid(&'a str),
	TidRegion(u32, u64, u64),
	NamedTidRegion(&'a str, u64, u64),
	//All,
	//Mapped,
	//Unmapped,
}

pub(crate) struct RegionWithLimits
{
	region: Vec<Range<u64>>,
	limits: Option<Range<u64>>,
	tid: u32,
}

pub struct Builder<R>
where
	R: AsyncReadSeek + std::marker::Send + std::marker::Unpin,
{
	reader: TokioBufReader<R>,
	regions: Vec<RegionWithLimits>,
	header: crate::bam::Header,
	bai_reader: Option<crate::bai::Reader>,
}

impl Builder<TokioFile>
{
	pub async fn from_path<P>(path: P) -> error::Result<Self>
	where
		P: AsRef<Path> + std::marker::Copy,
	{
		let file = TokioFile::open(path)
			.await
			.map_err(|_| error::Error::IOError(path.as_ref().to_string_lossy().to_string()))?;

		let mut reader = TokioBufReader::new(file);

		let mut bai_file = path.as_ref().to_path_buf();
		bai_file.add_extension("bai");

		//debug!("{:?}", bai_file.as_path().file_name().unwrap());

		let bai_reader = if bai_file.exists()
		{
			debug!("Setting up for indexed reading");
			Some(bai::Reader::from_path(&bai_file).await?)
		}
		else
		{
			debug!("BAI file doesn't exist, disabling BAM search");
			None
		};

		let header = crate::bam::read_bam_header(&mut reader).await?;

		Ok(Builder {
			reader,
			regions: Vec::new(),
			header,
			bai_reader,
		})
	}
}

impl<R> Builder<R>
where
	R: AsyncReadSeek + std::marker::Send + std::marker::Unpin,
{
	pub async fn from_reader(reader: R, bai_reader: Option<R>) -> error::Result<Self>
	{
		let mut async_reader = TokioBufReader::new(reader);

		let bai_reader = match bai_reader
		{
			Some(reader) =>
			{
				debug!("Setting up for indexed reading");
				Some(bai::Reader::from_reader(reader).await?)
			}
			None =>
			{
				debug!("BAI file doesn't exist, disabling BAM search");
				None
			}
		};

		/*let bai_reader = if bai_file.exists()
		{
			debug!("Setting up for indexed reading");
			Some(bai::Reader::from_path(&bai_file).await?)
		}
		else
		{
			debug!("BAI file doesn't exist, disabling BAM search");
			None
		};*/

		let header = crate::bam::read_bam_header(&mut async_reader).await?;

		Ok(Builder {
			reader: async_reader,
			regions: Vec::new(),
			header,
			bai_reader,
		})
	}

	pub fn add_fetch_region(&mut self, fetch: FetchRegion) -> error::Result<&mut Self>
	{
		match &self.bai_reader
		{
			Some(bai_reader) =>
			{
				match fetch
				{
					FetchRegion::Tid(tid) =>
					{
						let chunks = bai_reader.ref_indices[tid as usize].flatten_all_bins();
						self.regions.push(RegionWithLimits {
							region: chunks,
							limits: None,
							tid,
						});
					}
					FetchRegion::NamedTid(seqname) =>
					{
						let tid = self
							.header
							.references
							.iter()
							.position(|name| name.name == seqname)
							.unwrap();

						let chunks = bai_reader.ref_indices[tid as usize].flatten_all_bins();
						self.regions.push(RegionWithLimits {
							region: chunks,
							limits: None,
							tid: tid as u32,
						});
					}
					FetchRegion::TidRegion(tid, start, end) =>
					{
						let bins = bai::region_to_bins(start, end);
						let chunks = bai_reader.ref_indices[tid as usize].flatten_bins(&bins);
						self.regions.push(RegionWithLimits {
							region: chunks,
							limits: Some(Range { start, end }),
							tid: tid as u32,
						});
					}
					FetchRegion::NamedTidRegion(seqname, start, end) =>
					{
						let tid = self
							.header
							.references
							.iter()
							.position(|name| name.name == seqname)
							.unwrap();

						let bins = bai::region_to_bins(start, end);
						let chunks = bai_reader.ref_indices[tid as usize].flatten_bins(&bins);
						self.regions.push(RegionWithLimits {
							region: chunks,
							limits: Some(Range { start, end }),
							tid: tid as u32,
						});
					}
				};
			}
			None =>
			{
				return Err(error::Error::NoIndex);
			}
		}

		Ok(self)
	}

	pub async fn fetch_reads<F>(&mut self, mut read_fn: F) -> error::Result<Vec<bam::Pileup>>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let mut reads = Vec::<bam::Field>::new();

		for region in &self.regions
		{
			let mut coverage = BTreeSet::<u64>::new();

			for offset in &region.region
			{
				if let Some(limits) = &region.limits
				{
					if (limits.start..limits.end).all(|pos| coverage.contains(&pos))
					{
						break;
					}
				}

				let (block_start, virtual_start) = (offset.start >> 16, offset.start & 0xFFFF);
				let (block_end, virtual_end) = (offset.end >> 16, offset.end & 0xFFFF);

				debug!("virtual_start = {}", virtual_start);
				debug!("virtual_end = {}", virtual_end);

				self.reader
					.seek(SeekFrom::Start(block_start))
					.await
					.map_err(|_| error::Error::BamSeek(block_start))?;

				if block_start == block_end
				{
					let mut bytes = match self
						.reader
						.read_bgzf_block(Some(pufferfish::is_bam_eof))
						.await
					{
						Ok(bytes) => match bytes
						{
							Some(bytes) => bytes,
							None => break,
						},
						Err(err) => return Err(err.into()),
					};

					let mut seek_bytes = vec![0u8; (virtual_end - virtual_start) as usize];
					let mut cursor = Cursor::new(&mut bytes);
					Seek::seek(&mut cursor, SeekFrom::Start(virtual_start))
						.map_err(|_| error::Error::BamFormat)?;
					std::io::Read::read_exact(&mut cursor, &mut seek_bytes)
						.map_err(|_| error::Error::BamFormat)?;

					match bam::read_bam_block(
						&mut bytes,
						&mut reads,
						//&mut pileup,
						&mut read_fn,
						Some(region.tid),
						&region.limits,
						&mut coverage,
						&mut self.header,
					)
					.await?
					{
						Some(_bam_block) =>
						{}
						None =>
						{
							debug!("BAM EOF!");
							break;
						}
					}
				}
				else
				{
					let mut current_position = block_start;

					while current_position <= block_end
					{
						let start_offset = if current_position == block_start
						{
							virtual_start
						}
						else
						{
							0
						};

						let end_offset = if current_position == block_end
						{
							virtual_end
						}
						else
						{
							u64::MAX
						};

						let mut bytes = match self
							.reader
							.read_bgzf_block(Some(pufferfish::is_bam_eof))
							.await
						{
							Ok(bytes) => match bytes
							{
								Some(bytes) => bytes,
								None => break,
							},
							Err(err) => return Err(err.into()),
						};

						let block_len = bytes.len() as u64;
						let start_offset = start_offset.min(block_len);
						let end_offset = end_offset.min(block_len);

						if start_offset < end_offset
						{
							let mut seek_bytes = vec![0u8; (end_offset - start_offset) as usize];
							let mut cursor = Cursor::new(&mut bytes);
							Seek::seek(&mut cursor, SeekFrom::Start(start_offset))
								.map_err(|_| error::Error::BamFormat)?;
							std::io::Read::read_exact(&mut cursor, &mut seek_bytes)
								.map_err(|_| error::Error::BamFormat)?;
						}

						match bam::read_bam_block(
							&mut bytes,
							&mut reads,
							//&mut pileup,
							&mut read_fn,
							Some(region.tid),
							&region.limits,
							&mut coverage,
							&mut self.header,
						)
						.await?
						{
							Some(_bam_block) =>
							{}
							None =>
							{
								debug!("BAM EOF!");
								break;
							}
						}

						current_position = self
							.reader
							.seek(SeekFrom::Current(0))
							.await
							.map_err(|_| error::Error::BamFormat)?;
					}
				};
			}
		}

		let pileup = bam::generate_pileup(&reads);
		let mut pileup: Vec<_> = pileup
			.into_iter()
			.map(|((tid, pos), score)| bam::pileup::Pileup { tid, pos, score })
			.collect();
		pileup.sort_by(|a, b| (a.tid, a.pos).cmp(&(b.tid, b.pos)));

		Ok(pileup)
	}
}
