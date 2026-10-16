use std::path::Path;
use std::ops::Range;

// use std::io::{Seek, Cursor};
use std::collections::BTreeSet;

use tokio::fs::File as TokioFile;
use tokio::io::{SeekFrom, AsyncRead, AsyncSeek, BufReader as TokioBufReader};

use std::sync::Arc;

use tracing::{error, debug};

use crate::error;
use crate::bai;
use crate::bam;
// use crate::AsyncReadSeek;

use rayon::prelude::*;
use rayon::yield_now;
use rayon::{ThreadPoolBuilder, ThreadPool};

use futures::stream::Buffered;
use futures::StreamExt;

use crate::bam::CoverageProfile;

const DEFAULT_BUFFER_SIZE: usize = 500;
const DEFAULT_BLOCK_SIZE: usize = 500;

use pufferfish::prelude::*;

#[derive(Debug)]
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

#[derive(Debug)]
pub(crate) struct RegionWithLimits
{
	region: Vec<Range<u64>>,
	limits: Option<Range<u64>>,
	tid: u32,
}

pub struct ReaderOptions
{
	pub buffer_size: Option<usize>,
	pub n_threads: Option<usize>,
	pub n_blocks: Option<usize>,
}

impl Default for ReaderOptions
{
	fn default() -> Self
	{
		Self {
			buffer_size: Some(DEFAULT_BUFFER_SIZE),
			n_threads: None,
			n_blocks: Some(DEFAULT_BLOCK_SIZE),
		}
	}
}

impl ReaderOptions
{
	pub fn with_buffer_size(mut self, buffer_size: usize) -> Self
	{
		self.buffer_size = Some(buffer_size);
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
	// reader: TokioBufReader<R>,
	stream: Buffered<BgzfBlockStream<R>>,
	regions: Vec<RegionWithLimits>,
	header: crate::bam::Header,
	bai_reader: Option<crate::bai::Reader>,
	features: Option<crate::bam::BamFeatures>,
	thread_pool: ThreadPool,
	pending_tail: Option<Vec<u8>>,
	block_size: usize,
}

impl AsyncBamQuery<TokioFile>
{
	pub async fn from_path<P>(path: P, pool: Arc<pool::BgzfBlockPool>) -> error::Result<Self>
	where
		P: AsRef<Path>,
	{
		let path = path.as_ref().to_path_buf();
		let mut file = TokioFile::open(&path).await?;

		// let mut reader = TokioBufReader::new(&mut file);

		let mut bai_file = path.clone();
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

		let header = crate::bam::read_bam_header(&mut file).await?;

		let stream = BgzfBlockStream::new(file, pool.clone(), Some(is_bam_eof))
			.buffered(DEFAULT_BUFFER_SIZE);

		let builder = ThreadPoolBuilder::new();
		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			pending_tail: None,
			block_size: DEFAULT_BLOCK_SIZE,
		})
	}

	pub async fn from_path_with_options<P>(
		path: P,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	where
		P: AsRef<Path>,
	{
		let path = path.as_ref().to_path_buf();
		let mut file = TokioFile::open(&path).await?;

		// let mut reader = TokioBufReader::new(&mut file);

		let mut bai_file = path.clone();
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

		let header = crate::bam::read_bam_header(&mut file).await?;

		let stream = BgzfBlockStream::new(file, pool.clone(), Some(is_bam_eof))
			.buffered(options.buffer_size.unwrap_or(DEFAULT_BUFFER_SIZE));

		let builder = ThreadPoolBuilder::new();
		let builder = if let Some(n) = options.n_threads
		{
			builder.num_threads(n)
		}
		else
		{
			builder
		};

		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			pending_tail: None,
			block_size: options.n_blocks.unwrap_or(DEFAULT_BLOCK_SIZE),
		})
	}
}

impl<R> AsyncBamQuery<R>
where
	R: AsyncRead + AsyncSeek + std::marker::Send + std::marker::Unpin + 'static,
{
	pub async fn from_reader(
		mut reader: R,
		bai_reader: impl Into<Option<R>>,
		pool: Arc<pool::BgzfBlockPool>,
	) -> error::Result<Self>
	{
		let bai_reader = match bai_reader.into()
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

		let header = crate::bam::read_bam_header(&mut reader).await?;

		let stream = BgzfBlockStream::new(reader, pool.clone(), Some(is_bam_eof))
			.buffered(DEFAULT_BUFFER_SIZE);

		let builder = ThreadPoolBuilder::new();
		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			pending_tail: None,
			block_size: DEFAULT_BLOCK_SIZE,
		})
	}

	pub async fn from_reader_with_options(
		mut reader: R,
		bai_reader: impl Into<Option<R>>,
		pool: Arc<pool::BgzfBlockPool>,
		options: ReaderOptions,
	) -> error::Result<Self>
	{
		let bai_reader = match bai_reader.into()
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

		let header = crate::bam::read_bam_header(&mut reader).await?;

		let stream = BgzfBlockStream::new(reader, pool.clone(), Some(is_bam_eof))
			.buffered(options.buffer_size.unwrap_or(DEFAULT_BUFFER_SIZE));

		let builder = ThreadPoolBuilder::new();
		let builder = if let Some(n) = options.n_threads
		{
			builder.num_threads(n)
		}
		else
		{
			builder
		};

		let thread_pool = builder.build()?;

		Ok(AsyncBamQuery {
			stream,
			regions: Vec::new(),
			header,
			bai_reader,
			features: None,
			thread_pool,
			pending_tail: None,
			block_size: options.n_blocks.unwrap_or(DEFAULT_BLOCK_SIZE),
		})
	}

	pub fn set_features(&mut self, features: bam::BamFeatures) -> &mut Self
	{
		self.features = Some(features);
		self
	}

	pub async fn set_bai(&mut self, bai_reader: TokioBufReader<R>) -> error::Result<&mut Self>
	{
		self.bai_reader = Some(bai::Reader::from_reader(bai_reader).await?);

		Ok(self)
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

		debug!("regions = {:?}", regions);
		let mut reads = Vec::new();

		let stream = &mut self.stream;
		let header = &mut self.header;
		let thread_pool = &mut self.thread_pool;
		let features = self.features.as_ref();

		Self::_pileup_regions(
			stream,
			header,
			thread_pool,
			self.block_size,
			features,
			&mut read_fn,
			&mut reads,
			&regions,
		)
		.await?;

		if self
			.features
			.map_or(false, |f| f.contains(crate::bam::BamFeatures::PILEUP))
		{
			return Some(self.sort_pileup(reads)).transpose();
		}

		Ok(None)
	}

	pub async fn query_reads<F>(&mut self, mut read_fn: F) -> error::Result<Option<CoverageProfile>>
	// Return None when CIGAR ops are turned off since this is required for pileup generation
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let mut reads = Vec::<bam::Field>::new();

		if self.regions.is_empty()
		{
			self.pileup_all(&mut read_fn, &mut reads).await?;
		}
		else
		{
			self.pileup_regions(&mut read_fn, &mut reads).await?;
		}

		if self
			.features
			.map_or(false, |f| f.contains(crate::bam::BamFeatures::PILEUP))
		{
			return Some(self.sort_pileup(reads)).transpose();
		}

		Ok(None)
	}

	// async fn decompress_blocks<S>(
	// 	stream: &mut Buffered<BgzfBlockStream<S>>,
	// 	header: &mut crate::bam::Header,
	// 	pending_tail: &mut Option<BgzfCompressedBlock>,
	// 	thread_pool: &ThreadPool,
	// 	block_size: usize,
	// ) -> Option<BgzfCompressedBlock>
	// where
	// 	S: AsyncRead + AsyncSeek + std::marker::Send + std::marker::Unpin + 'static,
	// {
	// 	let batch: pufferfish::error::Result<BgzfCompressedBlock> = stream
	// 		.by_ref()
	// 		.take(block_size)
	// 		.collect::<pufferfish::error::Result<BgzfCompressedBlock>>()
	// 		.await;

	// 	let batch: Vec<Option<BgzfCompressedBlock>> = batch
	// 		.into_iter()
	// 		.filter_map(|res| match res
	// 		{
	// 			Ok(Some(block)) => Some(block), // keep the valid block
	// 			_ => None,                      // skip empty blocks
	// 		})
	// 		.collect();

	// 	let batch = if batch.is_empty()
	// 	{
	// 		match pending_tail.take()
	// 		{
	// 			Some(tail) if !tail.is_empty() => vec![tail],
	// 			_ => Vec::<BgzfCompressedBlock>::new(),
	// 		}
	// 	}
	// 	else
	// 	{
	// 		Vec::<BgzfCompressedBlock>::new()
	// 	};

	// 	if batch.is_empty()
	// 	{
	// 		return None;
	// 	}

	// 	let decompressed: Vec<Vec<u8>> = thread_pool.install(|| {
	// 		batch
	// 			.into_par_iter()
	// 			.map(|block| match decompress_bgzf_block(&block)
	// 			{
	// 				Ok(bytes) =>
	// 				{
	// 					yield_now();
	// 					bytes
	// 				}
	// 				Err(e) =>
	// 				{
	// 					error!(error = %e, "Unable to decompress BGZF block");
	// 					Vec::new()
	// 				}
	// 			})
	// 			.collect()
	// 	});

	// 	let total_size: usize = decompressed.iter().map(|b| b.len()).sum();

	// 	let mut blocks = Vec::with_capacity(total_size);

	// 	for bytes in decompressed
	// 	{
	// 		blocks.extend_from_slice(&bytes);
	// 	}

	// 	if let Some(tail) = pending_tail.as_mut()
	// 	{
	// 		if !tail.is_empty()
	// 		{
	// 			let mut merged = Vec::with_capacity(tail.len() + blocks.len());
	// 			merged.extend_from_slice(tail);
	// 			merged.extend_from_slice(&blocks);
	// 			blocks = merged;
	// 			tail.clear();
	// 		}
	// 	}

	// 	if let Some(idx) = memchr::memrchr(b'\n', &blocks)
	// 	{
	// 		let remainder = blocks.split_off(idx + 1);
	// 		*pending_tail.get_or_insert(Vec::new()) = remainder;
	// 	}
	// 	else
	// 	{
	// 		*pending_tail.get_or_insert(Vec::new()) = blocks;
	// 		blocks = Vec::new();
	// 	}

	// 	Some(blocks)
	// }

	async fn pileup_all<F>(
		&mut self,
		mut read_fn: F,
		reads: &mut Vec<bam::Field>,
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut bam::Header) -> Option<bam::Field>,
	{
		loop
		{
			let batch: Vec<_> = self.stream.by_ref().take(self.block_size).collect().await;

			let batch: Vec<_> = batch
				.into_iter()
				.filter_map(|res| match res
				{
					Ok(Some(block)) => Some(block), // keep the valid block
					_ => None,                      // skip empty blocks
				})
				.collect();

			let batch = if batch.is_empty()
			{
				match self.pending_tail.take()
				{
					Some(tail) if !tail.is_empty() => vec![tail],
					_ => Vec::<BgzfCompressedBlock>::new(),
				}
			}
			else
			{
				batch
			};

			if batch.is_empty()
			{
				break;
			}

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

			if let Some(tail) = self.pending_tail.as_mut()
			{
				if !tail.is_empty()
				{
					let mut merged = Vec::with_capacity(tail.len() + blocks.len());
					merged.extend_from_slice(tail);
					merged.extend_from_slice(&blocks);
					blocks = merged;
					tail.clear();
				}
			}

			// if let Some(idx) = memchr::memrchr(b'\n', &blocks)
			// {
			// 	let remainder = blocks.split_off(idx + 1);
			// 	*self.pending_tail.get_or_insert(Vec::new()) = remainder;
			// }
			// else
			// {
			// *self.pending_tail.get_or_insert(Vec::new()) = blocks;
			// blocks = Vec::new();
			// }

			match bam::read_bam_block(
				&mut blocks,
				reads,
				&mut read_fn,
				self.features.as_ref(),
				None,  // No tid filter
				&None, // No position limits
				&mut None,
				&mut self.header,
			)
			.await?
			{
				Some(_) =>
				{}
				None => break, // EOF
			}
		}
		// loop
		// {
		// 	let mut bytes = match self
		// 		.reader
		// 		.read_and_decompress_bgzf_block(Some(is_bam_eof))
		// 		.await
		// 	{
		// 		Ok(Some(bytes)) => bytes,
		// 		Ok(None) => break,
		// 		Err(err) => return Err(err.into()),
		// 	};

		// 	match bam::read_bam_block(
		// 		&mut bytes,
		// 		reads,
		// 		&mut read_fn,
		// 		self.features.as_ref(),
		// 		None,  // No tid filter
		// 		&None, // No position limits
		// 		&mut None,
		// 		&mut self.header,
		// 	)
		// 	.await?
		// 	{
		// 		Some(_) =>
		// 		{}
		// 		None => break, // EOF
		// 	}
		// }

		Ok(())
	}

	async fn pileup_regions<F>(
		&mut self,
		mut read_fn: F,
		reads: &mut Vec<bam::Field>,
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let regions = &self.regions;
		let stream = &mut self.stream;
		let header = &mut self.header;
		let thread_pool = &mut self.thread_pool;
		let features = self.features.as_ref();

		Self::_pileup_regions(
			stream,
			header,
			thread_pool,
			self.block_size,
			features,
			&mut read_fn,
			reads,
			regions,
		)
		.await?;

		Ok(())
	}

	async fn _pileup_regions<F>(
		stream: &mut Buffered<BgzfBlockStream<R>>,
		header: &mut crate::bam::Header,
		thread_pool: &ThreadPool,
		block_size: usize,
		features: Option<&crate::bam::BamFeatures>,
		read_fn: &mut F,
		reads: &mut Vec<bam::Field>,
		regions: &[RegionWithLimits],
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		for region in regions
		{
			Self::_pileup_region(
				stream,
				header,
				thread_pool,
				block_size,
				features,
				&mut *read_fn,
				reads,
				region,
			)
			.await?;
		}

		Ok(())
	}

	async fn _pileup_region<F>(
		stream: &mut Buffered<BgzfBlockStream<R>>,
		header: &mut crate::bam::Header,
		thread_pool: &ThreadPool,
		_block_size: usize,
		features: Option<&crate::bam::BamFeatures>,
		read_fn: &mut F,
		reads: &mut Vec<bam::Field>,
		region: &RegionWithLimits,
	) -> error::Result<()>
	where
		F: FnMut(bam::Field, &mut crate::bam::Header) -> Option<bam::Field>,
	{
		let mut coverage = Some(BTreeSet::<u64>::new());

		for offset in &region.region
		{
			if let Some(limits) = &region.limits
			{
				let coverage_set = coverage.as_ref().expect("unable to unwrap coverage");

				if (limits.start..limits.end).all(|pos| coverage_set.contains(&pos))
				{
					break;
				}
			}

			let (block_start, virtual_start) = (offset.start >> 16, offset.start & 0xFFFF);
			let (block_end, virtual_end) = (offset.end >> 16, offset.end & 0xFFFF);

			debug!(
				virtual_start = virtual_start,
				virtual_end = virtual_end,
				"Seeking to virtual offsets"
			);

			stream
				.get_mut()
				.seek_to(SeekFrom::Start(block_start))
				.await?;

			// if block_start == block_end
			// {
			// 	let batch: Vec<_> = stream.by_ref().take(1).collect().await;

			// 	let batch: Vec<_> = batch
			// 		.into_iter()
			// 		.filter_map(|res| match res
			// 		{
			// 			Ok(Some(block)) => Some(block), // keep the valid block
			// 			_ => None,                      // skip empty blocks
			// 		})
			// 		.collect();

			// 	if batch.is_empty()
			// 	{
			// 		break;
			// 	}

			// 	let decompressed: Vec<Vec<u8>> = thread_pool.install(|| {
			// 		batch
			// 			.into_par_iter()
			// 			.map(|block| {
			// 				match pufferfish::bgzf::stream::decompress_bgzf_block(&block)
			// 				{
			// 					Ok(bytes) =>
			// 					{
			// 						yield_now();
			// 						bytes
			// 					}
			// 					Err(e) =>
			// 					{
			// 						error!(error = %e, "Unable to decompress BGZF block");
			// 						Vec::new()
			// 					}
			// 				}
			// 			})
			// 			.collect()
			// 	});

			// 	let total_size: usize = decompressed.iter().map(|b| b.len()).sum();

			// 	let mut blocks = Vec::with_capacity(total_size);

			// 	for bytes in decompressed
			// 	{
			// 		blocks.extend_from_slice(&bytes);
			// 	}

			// 	// let mut bytes = match reader
			// 	// 	.read_and_decompress_bgzf_block(Some(is_bam_eof))
			// 	// 	.await
			// 	// {
			// 	// 	Ok(Some(bytes)) => bytes,
			// 	// 	Ok(None) => break,
			// 	// 	Err(err) => return Err(err.into()),
			// 	// };

			// 	match bam::read_bam_block(
			// 		&mut blocks[virtual_start as usize..virtual_end as usize],
			// 		reads,
			// 		//&mut pileup,
			// 		&mut *read_fn,
			// 		features,
			// 		Some(region.tid),
			// 		&region.limits,
			// 		&mut coverage,
			// 		&mut *header,
			// 	)
			// 	.await?
			// 	{
			// 		Some(_bam_block) =>
			// 		{}
			// 		None =>
			// 		{
			// 			debug!("BAM EOF!");
			// 			break;
			// 		}
			// 	}
			// }
			// else
			// {
			// let mut current_position = block_start;

			let block_count = (block_end - block_start + 1) as usize;

			let batch: Vec<_> = stream.by_ref().take(block_count).collect().await;

			let batch: Vec<_> = batch
				.into_iter()
				.filter_map(|res| match res
				{
					Ok(Some(block)) => Some(block), // keep the valid block
					_ => None,                      // skip empty blocks
				})
				.collect();

			if batch.is_empty()
			{
				break;
			}

			let decompressed: Vec<Vec<u8>> = thread_pool.install(|| {
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

			// let mut bytes = match reader
			// 	.read_and_decompress_bgzf_block(Some(is_bam_eof))
			// 	.await
			// {
			// 	Ok(Some(bytes)) => bytes,
			// 	Ok(None) => break,
			// 	Err(err) => return Err(err.into()),
			// };

			let block_len = blocks.len() as u64;
			let start_offset = virtual_start.min(block_len);
			let end_offset = virtual_end.min(block_len);

			blocks = if start_offset < end_offset
			{
				blocks[start_offset as usize..end_offset as usize].to_vec()
			}
			else
			{
				Vec::new()
			};

			// if start_offset < end_offset
			// {
			// 	let mut seek_bytes = vec![0u8; (end_offset - start_offset) as usize];
			// 	let mut cursor = Cursor::new(&mut blocks);
			// 	Seek::seek(&mut cursor, SeekFrom::Start(start_offset))
			// 		.map_err(|_| error::Error::BamFormat)?;
			// 	std::io::Read::read_exact(&mut cursor, &mut seek_bytes)
			// 		.map_err(|_| error::Error::BamFormat)?;

			// 	blocks = seek_bytes;
			// }

			match bam::read_bam_block(
				&mut blocks,
				reads,
				//&mut pileup,
				&mut *read_fn,
				features,
				Some(region.tid),
				&region.limits,
				&mut coverage,
				&mut *header,
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

			// current_position = stream.get_mut().seek_to(SeekFrom::Current(0)).await?;
			// current_position = reader
			// 	.seek(SeekFrom::Current(0))
			// 	.await
			// 	.map_err(|_| error::Error::BamFormat)?;
			// };
		}

		Ok(())
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
