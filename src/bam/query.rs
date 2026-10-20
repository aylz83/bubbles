use std::path::Path;
use std::ops::Range;

// use std::io::{Seek, Cursor};
// use std::collections::BTreeSet;

use tokio::fs::File as TokioFile;
use tokio::io::{SeekFrom, AsyncRead, AsyncSeek};

use std::sync::Arc;

use tracing::{error, debug};

use crate::error;
use crate::bai;
use crate::bam;
// use crate::AsyncReadSeek;

use rayon::prelude::*;
use rayon::yield_now;
use rayon::{ThreadPoolBuilder, ThreadPool};

// use futures::stream::Buffered;
use futures::StreamExt;

use crate::bam::CoverageProfile;

const DEFAULT_READ_AHEAD_SIZE: usize = 2;
const DEFAULT_BLOCK_SIZE: usize = 200;

use pufferfish::prelude::*;

#[derive(Debug, Clone)]
pub enum FetchRegion
{
	Tid(u32),
	NamedTid(String),
	TidRegion(u32, u64, u64),
	NamedTidRegion(String, u64, u64),
	//Mapped,
	//Unmapped,
}

impl From<u32> for FetchRegion
{
	fn from(tid: u32) -> Self
	{
		Self::Tid(tid)
	}
}

impl From<&str> for FetchRegion
{
	fn from(name: &str) -> Self
	{
		Self::NamedTid(name.to_owned())
	}
}

impl From<String> for FetchRegion
{
	fn from(name: String) -> Self
	{
		Self::NamedTid(name)
	}
}

impl<'a> From<std::borrow::Cow<'a, str>> for FetchRegion
{
	fn from(name: std::borrow::Cow<'a, str>) -> Self
	{
		Self::NamedTid(name.into_owned().to_string())
	}
}

impl<'a> From<std::borrow::Cow<'a, bstr::BStr>> for FetchRegion
{
	fn from(name: std::borrow::Cow<'a, bstr::BStr>) -> Self
	{
		Self::NamedTid(name.into_owned().to_string())
	}
}

#[derive(Debug, Clone)]
pub(crate) struct RegionWithLimits
{
	pub(crate) region: Vec<Range<u64>>,
	pub(crate) limits: Option<Range<u64>>,
	pub(crate) tid: u32,
}

pub struct ReaderOptions
{
	pub read_ahead_size: Option<usize>,
	pub n_threads: Option<usize>,
	pub n_blocks: Option<usize>,
}

impl Default for ReaderOptions
{
	fn default() -> Self
	{
		Self {
			read_ahead_size: Some(DEFAULT_READ_AHEAD_SIZE),
			n_threads: None,
			n_blocks: Some(DEFAULT_BLOCK_SIZE),
		}
	}
}

impl ReaderOptions
{
	pub fn with_read_ahead_size(mut self, read_ahead_size: usize) -> Self
	{
		self.read_ahead_size = Some(read_ahead_size);
		self
	}

	pub fn with_n_threads(mut self, n: usize) -> Self
	{
		self.n_threads = Some(n);
		self
	}

	pub fn with_block_read_size(mut self, n: usize) -> Self
	{
		self.n_blocks = Some(n);
		self
	}
}

pub struct AsyncBamQuery<R>
where
	R: AsyncRead + AsyncSeek + std::marker::Send + std::marker::Unpin + 'static,
{
	stream: BgzfBlockStream<R>,
	regions: Vec<RegionWithLimits>,
	header: crate::bam::Header,
	bai_reader: Option<crate::bai::Reader>,
	features: Option<crate::bam::BamFeatures>,
	thread_pool: ThreadPool,
	read_ahead_size: usize,
	block_size: usize,
}

impl AsyncBamQuery<TokioFile>
{
	pub async fn from_path<P>(
		path: P,
		bai_path: impl Into<Option<P>>,
		pool: Arc<pool::BgzfBlockPool>,
	) -> error::Result<Self>
	where
		P: AsRef<Path>,
	{
		Self::_from_path(path, bai_path.into(), pool, ReaderOptions::default()).await
	}

	pub async fn from_path_with_options<P>(
		path: P,
		bai_path: impl Into<Option<P>>,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	where
		P: AsRef<Path>,
	{
		Self::_from_path(path, bai_path.into(), pool, options).await
	}

	async fn _from_path(
		path: impl AsRef<Path>,
		bai_path: Option<impl AsRef<Path>>,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	{
		let path = path.as_ref().to_path_buf();
		let mut file = TokioFile::open(&path).await?;

		let bai_file = bai_path
			.map(|p| p.as_ref().to_path_buf())
			.unwrap_or_else(|| {
				let mut path = path.clone();
				path.add_extension("bai");
				path
			});

		let bai_reader = if bai_file.exists()
		{
			debug!(path = ?bai_file, "Setting up for indexed reading");
			Some(bai::Reader::from_path(&bai_file).await?)
		}
		else
		{
			debug!("BAI file doesn't exist, disabling BAM search");
			None
		};

		let header = crate::bam::read_bam_header(&mut file).await?;
		let stream = BgzfBlockStream::new(file, pool, Some(is_bam_eof));

		let mut builder = ThreadPoolBuilder::new();
		if let Some(n) = options.n_threads
		{
			builder = builder.num_threads(n);
		}

		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			read_ahead_size: options.read_ahead_size.unwrap_or(DEFAULT_READ_AHEAD_SIZE),
			block_size: options.n_blocks.unwrap_or(DEFAULT_BLOCK_SIZE),
		})
	}
}

impl<R> AsyncBamQuery<R>
where
	R: AsyncRead + AsyncSeek + std::marker::Send + std::marker::Unpin + 'static,
{
	pub async fn from_reader(
		reader: R,
		bai_reader: impl Into<Option<R>>,
		pool: Arc<pool::BgzfBlockPool>,
	) -> error::Result<Self>
	{
		Self::_from_reader(reader, bai_reader.into(), pool, ReaderOptions::default()).await
	}

	pub async fn from_reader_with_options(
		reader: R,
		bai_reader: impl Into<Option<R>>,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	{
		Self::_from_reader(reader, bai_reader.into(), pool, options).await
	}

	async fn _from_reader(
		mut reader: R,
		bai_reader: Option<R>,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	{
		let bai_reader = match bai_reader
		{
			Some(reader) =>
			{
				debug!("Setting up for indexed reading");
				Some(bai::Reader::from_reader(reader).await?)
			}
			None =>
			{
				debug!("No BAI reader supplied, disabling indexed search");
				None
			}
		};

		let header = crate::bam::read_bam_header(&mut reader).await?;
		let stream = BgzfBlockStream::new(reader, pool, Some(is_bam_eof));

		let mut builder = ThreadPoolBuilder::new();

		if let Some(n) = options.n_threads
		{
			builder = builder.num_threads(n);
		}

		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			read_ahead_size: options.read_ahead_size.unwrap_or(DEFAULT_READ_AHEAD_SIZE),
			block_size: options.n_blocks.unwrap_or(DEFAULT_BLOCK_SIZE),
		})
	}

	pub fn set_features(&mut self, features: bam::BamFeatures) -> &mut Self
	{
		self.features = Some(features);
		self
	}

	// pub async fn set_bai(&mut self, bai_reader: TokioBufReader<R>) -> error::Result<&mut Self>
	// {
	// 	self.bai_reader = Some(bai::Reader::from_reader(bai_reader).await?);

	// 	Ok(self)
	// }

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
							.position(|name| name.name_as_str() == seqname)
							.ok_or(error::Error::NoRef(seqname.to_string()))?;

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
							.position(|name| name.name_as_str() == seqname)
							.ok_or(error::Error::NoRef(seqname.to_string()))?;

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

	pub async fn oneshot_query<F>(
		&mut self,
		fetch_regions: &[FetchRegion],
		mut read_fn: F,
	) -> error::Result<Option<CoverageProfile>>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let bai_reader = self.bai_reader.as_ref().ok_or(error::Error::NoIndex)?;

		let mut regions = Vec::with_capacity(fetch_regions.len());

		for fetch in fetch_regions
		{
			let (tid, limits, chunks) = match fetch
			{
				FetchRegion::Tid(tid) =>
				{
					let chunks = bai_reader.ref_indices[*tid as usize].flatten_all_bins();

					(*tid, None, chunks)
				}
				FetchRegion::NamedTid(seqname) =>
				{
					let tid = self
						.header
						.references
						.iter()
						.position(|name| name.name_as_str() == seqname)
						.ok_or(error::Error::NoRef(seqname.to_string()))? as u32;

					let chunks = bai_reader.ref_indices[tid as usize].flatten_all_bins();

					(tid, None, chunks)
				}
				FetchRegion::TidRegion(tid, start, end) =>
				{
					let bins = bai::region_to_bins(*start, *end);
					let chunks = bai_reader.ref_indices[*tid as usize].flatten_bins(&bins);

					(
						*tid,
						Some(Range {
							start: *start,
							end: *end,
						}),
						chunks,
					)
				}
				FetchRegion::NamedTidRegion(seqname, start, end) =>
				{
					let tid = self
						.header
						.references
						.iter()
						.position(|name| name.name_as_str() == seqname)
						.ok_or(error::Error::NoRef(seqname.to_string()))? as u32;

					let bins = bai::region_to_bins(*start, *end);
					let chunks = bai_reader.ref_indices[tid as usize].flatten_bins(&bins);

					(
						tid,
						Some(Range {
							start: *start,
							end: *end,
						}),
						chunks,
					)
				}
			};

			regions.push(RegionWithLimits {
				region: chunks,
				limits,
				tid,
			});
		}

		let pileup_enabled = self
			.features
			.is_some_and(|f| f.contains(crate::bam::BamFeatures::PILEUP));

		let mut reads = Vec::<bam::Field>::new();

		let stream = &mut self.stream;
		let header = &mut self.header;
		let thread_pool = &mut self.thread_pool;
		let features = self.features.as_ref();

		Self::pileup_regions(
			stream,
			header,
			thread_pool,
			self.read_ahead_size,
			self.block_size,
			features,
			&mut read_fn,
			if pileup_enabled
			{
				Some(&mut reads)
			}
			else
			{
				None
			},
			&regions,
		)
		.await?;

		if pileup_enabled
		{
			Ok(Some(self.sort_pileup(reads)?))
		}
		else
		{
			Ok(None)
		}
	}

	pub async fn query_all_reads<F>(
		&mut self,
		mut read_fn: F,
	) -> error::Result<Option<CoverageProfile>>
	// Return None when CIGAR ops are turned off since this is required for pileup generation
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field> + Send,
	{
		let pileup_enabled = self
			.features
			.is_some_and(|f| f.contains(crate::bam::BamFeatures::PILEUP));

		let mut reads = Vec::<bam::Field>::new();

		if self.regions.is_empty()
		{
			self.pileup_all(
				&mut read_fn,
				if pileup_enabled
				{
					Some(&mut reads)
				}
				else
				{
					None
				},
			)
			.await?;
		}
		else
		{
			let regions = &self.regions;
			let stream = &mut self.stream;
			let header = &mut self.header;
			let thread_pool = &mut self.thread_pool;
			let features = self.features.as_ref();

			Self::pileup_regions(
				stream,
				header,
				thread_pool,
				self.read_ahead_size,
				self.block_size,
				features,
				&mut read_fn,
				if pileup_enabled
				{
					Some(&mut reads)
				}
				else
				{
					None
				},
				regions,
			)
			.await?;
		}

		if self
			.features
			.map_or(false, |f| f.contains(crate::bam::BamFeatures::PILEUP))
		{
			return Some(self.sort_pileup(reads)).transpose();
		}

		Ok(None)
	}

	async fn pileup_all<F>(
		&mut self,
		mut read_fn: F,
		mut reads: Option<&mut Vec<bam::Field>>,
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut bam::Header) -> Option<bam::Field>,
	{
		let batch_size = self.block_size.max(1);
		let (tx, mut rx) = tokio::sync::mpsc::channel::<error::Result<Vec<BgzfCompressedBlock>>>(
			self.read_ahead_size,
		);

		let mut reader_stream = self.stream.clone();

		let producer = tokio::spawn(async move {
			loop
			{
				let mut batch = Vec::with_capacity(batch_size);
				let mut eof = false;

				for _ in 0..batch_size
				{
					let Some(future) = reader_stream.next().await
					else
					{
						eof = true;
						break;
					};

					match future.await
					{
						Ok(Some(block)) => batch.push(block),

						Ok(None) =>
						{
							eof = true;
							break;
						}

						Err(_e) =>
						{
							let _ = tx.send(Err(error::Error::BamFormat)).await;
							return;
						}
					}
				}

				if !batch.is_empty() && tx.send(Ok(batch)).await.is_err()
				{
					break;
				}

				if eof
				{
					break;
				}
			}
		});

		// parallel decompression, then sequential BAM parsing.
		while let Some(batch_result) = rx.recv().await
		{
			let batch = match batch_result
			{
				Ok(batch) => batch,
				Err(e) =>
				{
					producer.abort();
					let _ = producer.await;
					return Err(e);
				}
			};

			let decompressed: Vec<Vec<u8>> = self.thread_pool.install(|| {
				batch
					.into_par_iter()
					.map(
						|block| match pufferfish::bgzf::stream::decompress_bgzf_block(&block)
						{
							Ok(bytes) =>
							{
								yield_now();
								bytes
							}
							Err(e) =>
							{
								error!(error = %e, "Unable to decompress BGZF block");
								Vec::new()
							}
						},
					)
					.collect()
			});

			let total_size: usize = decompressed.iter().map(|b| b.len()).sum();
			let mut blocks = Vec::with_capacity(total_size);

			for bytes in decompressed
			{
				blocks.extend_from_slice(&bytes);
			}

			debug!(parse_blocks = blocks.len(), "Blocks to parse by BAM reader");

			match bam::read_bam_block(
				&mut blocks,
				reads.as_deref_mut(),
				&mut read_fn,
				self.features.as_ref(),
				None,
				&mut self.header,
			)
			.await
			{
				Ok(Some(_)) =>
				{}

				Ok(None) =>
				{
					debug!("BAM EOF!");
					drop(rx);
					producer.abort();
					break;
				}

				Err(e) =>
				{
					drop(rx);
					producer.abort();
					let _ = producer.await;
					return Err(e);
				}
			}
		}

		if let Err(e) = producer.await
		{
			if !e.is_cancelled()
			{
				error!(error = %e, "BGZF reader task failed");
				return Err(error::Error::BamFormat);
			}
		}

		Ok(())
	}

	async fn pileup_regions<F>(
		stream: &mut BgzfBlockStream<R>,
		header: &mut crate::bam::Header,
		thread_pool: &ThreadPool,
		read_ahead_size: usize,
		block_size: usize,
		features: Option<&crate::bam::BamFeatures>,
		read_fn: &mut F,
		mut reads: Option<&mut Vec<bam::Field>>,
		regions: &[RegionWithLimits],
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let mut merged_regions = regions.to_vec();

		for region in &mut merged_regions
		{
			Self::merge_ranges(&mut region.region);
		}

		let merged_length = merged_regions.len();

		let (tx, mut rx) = tokio::sync::mpsc::channel::<
			error::Result<Vec<(Range<usize>, BgzfCompressedBlock)>>,
		>(read_ahead_size);

		let mut reader_stream = stream.clone();

		let merged = merged_regions.clone();
		let producer = tokio::spawn(async move {
			let mut region_index = 0usize;
			let mut offset_index = 0usize;
			let mut block_position: Option<u64> = None;

			while region_index < merged_length
			{
				let batch = Self::compressed_blocks(
					&mut reader_stream,
					&merged_regions,
					&mut region_index,
					&mut offset_index,
					&mut block_position,
					block_size,
				)
				.await;

				match batch
				{
					Ok(blocks) if blocks.is_empty() => break,

					Ok(blocks) =>
					{
						if tx.send(Ok(blocks)).await.is_err()
						{
							break;
						}
					}

					Err(e) =>
					{
						let _ = tx.send(Err(e)).await;
						break;
					}
				}
			}
		});

		while let Some(batch_result) = rx.recv().await
		{
			let compressed_blocks = match batch_result
			{
				Ok(blocks) => blocks,
				Err(e) =>
				{
					producer.abort();
					let _ = producer.await;
					return Err(e);
				}
			};

			let mut blocks = Self::decompress_blocks(thread_pool, compressed_blocks);

			debug!(parse_blocks = blocks.len(), "Blocks to parse by BAM reader");

			match bam::read_bam_block(
				&mut blocks,
				reads.as_deref_mut(),
				&mut *read_fn,
				features,
				Some(&merged),
				&mut *header,
			)
			.await
			{
				Ok(Some(_)) =>
				{}

				Ok(None) =>
				{
					debug!("BAM EOF!");
					// Stop unnecessary read-ahead.
					drop(rx);
					producer.abort();
					break;
				}

				Err(e) =>
				{
					drop(rx);
					producer.abort();
					let _ = producer.await;
					return Err(e);
				}
			}
		}

		// Wait for the producer to finish on normal channel closure.
		// A producer aborted on BAM EOF is also handled here.
		if let Err(e) = producer.await
		{
			if !e.is_cancelled()
			{
				error!(error = %e, "BGZF reader task failed");
				return Err(error::Error::BamFormat);
			}
		}

		Ok(())
	}

	fn merge_ranges(ranges: &mut Vec<Range<u64>>)
	{
		ranges.sort_unstable_by_key(|range| (range.start, range.end));

		let mut merged: Vec<Range<u64>> = Vec::with_capacity(ranges.len());

		for range in ranges.drain(..)
		{
			if let Some(last) = merged.last_mut()
			{
				if range.start <= last.end
				{
					last.end = last.end.max(range.end);
					continue;
				}
			}

			merged.push(range);
		}

		*ranges = merged;
	}

	async fn compressed_blocks(
		stream: &mut BgzfBlockStream<R>,
		regions: &[RegionWithLimits],
		region_index: &mut usize,
		offset_index: &mut usize,
		block_position: &mut Option<u64>,
		block_size: usize,
	) -> error::Result<Vec<(Range<usize>, BgzfCompressedBlock)>>
	{
		let mut blocks = Vec::with_capacity(block_size);

		while *region_index < regions.len() && blocks.len() < block_size
		{
			let region = &regions[*region_index];

			if *offset_index >= region.region.len()
			{
				*region_index += 1;
				*offset_index = 0;
				*block_position = None;
				continue;
			}

			let offset = &region.region[*offset_index];
			let (block_start, virtual_start) =
				(offset.start >> 16, (offset.start & 0xFFFF) as usize);
			let (block_end, virtual_end) = (offset.end >> 16, (offset.end & 0xFFFF) as usize);

			if offset.start >= offset.end
			{
				*offset_index += 1;
				*block_position = None;
				debug!(offset_index = *offset_index);
				continue;
			}

			let last_block = if virtual_end == 0 && block_end > block_start
			{
				block_end - 1
			}
			else
			{
				block_end
			};

			if block_start > last_block
			{
				*offset_index += 1;
				*block_position = None;
				debug!(
					offset_index = *offset_index,
					"block_start is past last_block"
				);
				continue;
			}

			let position = *block_position.get_or_insert(block_start);

			if position > last_block
			{
				*offset_index += 1;
				*block_position = None;
				debug!(offset_index = *offset_index, "position is past last_block");
				continue;
			}

			stream.seek_to(SeekFrom::Start(position)).await?;

			let mut read_count = 0usize;

			while blocks.len() < block_size
			{
				let file_offset = stream.stream_position().await?;

				if file_offset > last_block
				{
					*offset_index += 1;
					*block_position = None;
					debug!(
						offset_index = *offset_index,
						"Completed virtual-offset range"
					);
					break;
				}

				let Some(future) = stream.next().await
				else
				{
					*region_index = regions.len();
					*block_position = None;
					break;
				};

				debug!(file_offset, "Reading BGZF block");

				match future.await
				{
					Ok(Some(block)) =>
					{
						let start = if file_offset == block_start
						{
							virtual_start
						}
						else
						{
							0
						};
						let end = if file_offset == block_end
						{
							virtual_end
						}
						else
						{
							usize::MAX
						};

						blocks.push((start..end, block));
						read_count += 1;

						// Read the actual position after this block; do not add 1.
						let next_position = stream.stream_position().await?;
						*block_position = Some(next_position);

						debug!(
							file_offset,
							next_position,
							blocks_in_batch = blocks.len(),
							"BGZF block appended"
						);

						if next_position > last_block
						{
							*offset_index += 1;
							*block_position = None;
							debug!(
								offset_index = *offset_index,
								"Completed virtual-offset range"
							);
							break;
						}
					}
					Ok(None) =>
					{
						debug!("BGZF EOF!");
						*region_index = regions.len();
						*block_position = None;
						break;
					}
					Err(e) =>
					{
						error!(error = ?e, "Unable to look at block");
						return Err(error::Error::BamFormat);
					}
				}
			}

			debug!(read_count, "Finished BGZF batch");

			if read_count == 0 && *region_index < regions.len()
			{
				// Avoid an infinite loop if no block could be read.
				*offset_index += 1;
				*block_position = None;
			}
		}

		Ok(blocks)
	}

	fn decompress_blocks(
		thread_pool: &ThreadPool,
		blocks: Vec<(Range<usize>, BgzfCompressedBlock)>,
	) -> Vec<u8>
	{
		thread_pool.install(|| {
			let decompressed: Vec<_> = blocks
				.into_par_iter()
				.map(|(range, block)| {
					match pufferfish::bgzf::stream::decompress_bgzf_block(&block)
					{
						Ok(bytes) =>
						{
							let start = range.start.min(bytes.len());
							let end = range.end.min(bytes.len()).max(start);

							debug!(
								virtual_start = start,
								virtual_end = end,
								"Seeking to virtual offsets"
							);

							(bytes, start..end)
						}
						Err(e) =>
						{
							error!(error = %e, "Unable to decompress BGZF block");
							(Vec::new(), 0..0)
						}
					}
				})
				.collect();

			let total_len: usize = decompressed
				.iter()
				.map(|(_bytes, range)| range.end - range.start)
				.sum();

			let mut output = Vec::with_capacity(total_len);

			for (bytes, range) in decompressed
			{
				output.extend_from_slice(&bytes[range]);
			}

			output
		})
	}

	fn sort_pileup(&self, reads: Vec<bam::Field>) -> error::Result<CoverageProfile>
	{
		let pileup = bam::generate_pileup(&reads);
		let mut pileup: Vec<_> = pileup
			.into_iter()
			.map(|((tid, pos), score)| bam::pileup::Pileup { tid, pos, score })
			.collect();
		pileup.sort_by(|a, b| (a.tid, a.pos).cmp(&(b.tid, b.pos)));

		Ok(CoverageProfile(pileup))
	}
}
