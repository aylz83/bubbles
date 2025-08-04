#![allow(dead_code)]

mod cigar;
mod header;
mod pileup;
mod query;
mod tags;

pub use crate::bam::query::*;
pub use crate::bam::pileup::*;
pub use crate::bam::cigar::*;
pub use crate::bam::tags::*;
pub use crate::bam::header::*;

use crate::error;

use log::debug;

use std::simd::Simd;
use std::ops::Range;
use std::collections::BTreeSet;

const ALPHABET: [u8; 16] = [
	b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

use bitflags::bitflags;

bitflags! {
	#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
	pub struct BamFeatures: u32 {
		const READNAMES = 0b0001;
		const SEQUENCES = 0b0010;
		const CIGAR = 0b0100;
		const TAGS = 0b1000;
		const PILEUP = 0b10000 | Self::CIGAR.bits();
		const EMPTY = 0b0000;
	}
}

pub struct Field
{
	pub ref_id: i32,
	pub pos: i32,
	pub mapq: u8,
	pub bin: u16,
	pub flags: u16,
	pub next_ref_id: i32,
	pub next_pos: i32,
	pub tlen: i32,

	pub read_name: Option<Box<[u8]>>,
	pub sequence: Option<Box<[u8]>>,
	pub sequence_quality: Option<Box<[u8]>>,
	pub cigar: Option<Vec<Cigar>>,
	pub tags: Option<Vec<Tag>>,
}

impl Field
{
	pub fn read_name_as_str(&self) -> Option<&str>
	{
		self.read_name.as_ref().map(|name| unsafe {
			std::str::from_utf8_unchecked(&name[..name.len().saturating_sub(1)])
		})
	}

	pub fn sequence_as_str(&self) -> Option<&str>
	{
		self.sequence.as_ref().map(|name| unsafe {
			std::str::from_utf8_unchecked(&name[..name.len().saturating_sub(1)])
		})
	}

	pub fn sequence_quality_as_str(&self) -> Option<&str>
	{
		self.sequence_quality.as_ref().map(|name| unsafe {
			std::str::from_utf8_unchecked(&name[..name.len().saturating_sub(1)])
		})
	}
}

pub(crate) async fn read_bam_block<F>(
	bytes: &mut [u8],
	reads: &mut Vec<Field>,
	mut read_fn: F,
	features: Option<&BamFeatures>,
	tid: Option<u32>,
	region: &Option<Range<u64>>,
	coverage: &mut Option<BTreeSet<u64>>,
	header: &mut Header,
) -> error::Result<Option<()>>
// Result due to some functions have the possibility to fail, and Option so we can detect None = EOF
where
	F: FnMut(Field, &mut Header) -> Option<Field>,
{
	let mut offset: usize = 0;
	let mut start_block_offset: usize;
	while offset < bytes.len()
	{
		start_block_offset = offset;

		let block_size = u32::from_le_bytes([
			bytes[offset],
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
		]);

		debug!("block_size = {}", block_size);

		offset += 3;

		// read_id - int32_t
		let ref_id = i32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);

		debug!("ref_id = {}", ref_id);

		offset += 4;

		// pos - int32_t
		let pos = i32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);

		debug!("pos = {}", pos);

		if let (Some(tid), Some(region)) = (tid, region)
		{
			if !(tid == ref_id as u32 && pos as u64 >= region.start && pos as u64 <= region.end)
			{
				debug!(
					"breaking, tid = {}, ref_tid = {}, region.start = {}, region.end = {}, pos = \
					 {}, block_size = {}",
					tid, ref_id, region.start, region.end, pos, block_size
				);
				offset += (block_size - 3 - 4 - 4) as usize;
				continue;
			}
		}
		else if let Some(tid) = tid
		{
			if !(tid == ref_id as u32)
			{
				offset += (block_size - 3 - 4 - 4) as usize;
				continue;
			}
		}

		offset += 4;

		// l_read_name - uint8_t
		let l_read_name = u8::from_le_bytes([bytes[offset + 1]]);
		debug!("l_read_name = {}", l_read_name);

		offset += 1;

		// mapq - uint8_t
		let mapq = u8::from_le_bytes([bytes[offset + 1]]);
		debug!("mapq = {}", mapq);

		offset += 1;

		// bin - uint16_t
		let bin = u16::from_le_bytes([bytes[offset + 1], bytes[offset + 2]]);
		debug!("bin = {}", bin);

		offset += 2;

		// n_cigar_op - uint16_t
		let n_cigar_op = u16::from_le_bytes([bytes[offset + 1], bytes[offset + 2]]);
		debug!("n_cigar_op = {}", n_cigar_op);

		offset += 2;

		// flag - uint16_t
		let flag = u16::from_le_bytes([bytes[offset + 1], bytes[offset + 2]]);
		debug!("flag = {}", flag);

		offset += 2;

		// l_seq - uint32_t
		let l_seq = u32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);
		debug!("l_seq = {}", l_seq);

		offset += 4;

		// next_ref_id - int32_t
		let next_ref_id = i32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);
		debug!("next_ref_id = {}", next_ref_id);

		offset += 4;

		// next_pos - int32_t
		let next_pos = i32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);
		debug!("next_pos = {}", next_pos);

		offset += 4;

		// tlen - int32_t
		let tlen = i32::from_le_bytes([
			bytes[offset + 1],
			bytes[offset + 2],
			bytes[offset + 3],
			bytes[offset + 4],
		]);
		debug!("tlen = {}", tlen);

		offset += 4;

		// read_name - char[l_read_name]
		let read_name = if features.map_or(false, |f| f.contains(BamFeatures::READNAMES))
		{
			let r_name = &bytes[offset + 1..offset + 1 + l_read_name as usize];
			Some(Box::from(r_name))
		}
		else
		{
			None
		};

		debug!("read_name = {:?}", read_name);

		offset += 1 + l_read_name as usize;

		// cigar - uint32_t[n_cigar_op]
		let mut ref_index = pos;
		let cigar = if features.map_or(false, |f| f.contains(BamFeatures::CIGAR))
		{
			Some(process_cigar(
				&bytes,
				&mut offset,
				n_cigar_op as usize,
				ref_id,
				&mut ref_index,
				//pileup_map,
				coverage,
			))
		}
		else
		{
			offset += n_cigar_op as usize * 4;
			None
		};

		debug!("cigar = {:?}", cigar);

		// seq - uint8_t[(l_seq + 1) / 2]
		let seq = if features.map_or(false, |f| f.contains(BamFeatures::SEQUENCES))
		{
			Some(process_sequence(
				&bytes[offset..(offset + l_seq as usize)],
				l_seq as usize,
			))
		}
		else
		{
			None
		};

		debug!("seq = {:?}", &seq);

		offset += (l_seq as usize + 1) / 2;

		// qual - char[l_seq]
		let qual = if features.map_or(false, |f| f.contains(BamFeatures::SEQUENCES))
		{
			Some(process_quality_scores(
				&bytes[offset..offset + l_seq as usize],
				l_seq as usize,
			))
		}
		else
		{
			None
		};

		debug!("qual = {:?}", qual);

		offset += l_seq as usize;

		let tags_vec = if features.map_or(false, |f| f.contains(BamFeatures::TAGS))
		{
			let tag_bytes = start_block_offset + block_size as usize - offset;
			let estimated_num_tags = tag_bytes / 8;

			let mut tags = Vec::<Tag>::with_capacity(estimated_num_tags);

			read_tags(
				&bytes,
				start_block_offset + block_size as usize,
				&mut offset,
				&mut tags,
			)?;

			Some(tags)
		}
		else
		{
			offset = start_block_offset + block_size as usize + 4;
			None
		};

		let read = Field {
			ref_id,
			pos,
			mapq,
			bin,
			flags: flag,
			next_ref_id,
			next_pos,
			tlen,
			read_name: read_name,
			sequence: seq,
			sequence_quality: qual,
			cigar: cigar,
			tags: tags_vec,
		};

		if let Some(keep_read) = read_fn(read, header)
		{
			reads.push(keep_read);
		}
	}

	Ok(Some(()))
}

fn process_sequence(bytes: &[u8], l_seq: usize) -> Box<[u8]>
{
	let mut seq = Vec::with_capacity(l_seq); // Allocate enough space for the sequence

	const CHUNK_SIZE: usize = 16; // SIMD width for u8

	let mut offset = 0; // Initialize the byte offset

	// Process the sequence in chunks
	while offset < (l_seq + 1) / 2
	{
		let remaining = (l_seq + 1) / 2 - offset; // Remaining bytes to process
		let to_process = remaining.min(CHUNK_SIZE);

		let mut byte_chunk = [0u8; CHUNK_SIZE];
		for i in 0..to_process
		{
			byte_chunk[i] = bytes[offset + i];
		}

		let bytes_simd = Simd::<u8, CHUNK_SIZE>::from_array(byte_chunk);

		// Extract high and low nibbles using SIMD
		let high_nibbles = (bytes_simd >> 4) & Simd::splat(0x0F); // High nibbles
		let low_nibbles = bytes_simd & Simd::splat(0x0F); // Low nibbles

		// Create an array to hold the corresponding bases
		let mut bases = Vec::with_capacity(to_process * 2);

		// Map the high nibbles to bases
		for i in 0..to_process
		{
			bases.push(ALPHABET[high_nibbles[i] as usize]);
			if 2 * i + 1 < l_seq
			{
				bases.push(ALPHABET[low_nibbles[i] as usize]);
			}
		}

		seq.extend_from_slice(&bases);

		// Update the offset
		offset += to_process;
	}

	// Trim the sequence to the exact length if necessary
	seq.truncate(l_seq);

	seq.into_boxed_slice()
}

fn process_quality_scores(bytes: &[u8], l_seq: usize) -> Box<[u8]>
{
	let mut qual = Vec::with_capacity(l_seq);

	// Process the quality scores in chunks
	const CHUNK_SIZE: usize = 16; // This may vary based on your architecture
	let mut i = 0;

	// Process the sequence in chunks of SIMD width
	while i < l_seq
	{
		// Determine how many bytes to process in this iteration
		let remaining = l_seq - i;
		let to_process = remaining.min(CHUNK_SIZE);

		// Load bytes into a SIMD register
		let mut byte_chunk = [0u8; CHUNK_SIZE]; // Create an array to hold the current chunk of bytes
		byte_chunk[..to_process].copy_from_slice(&bytes[i..i + to_process]);

		let bytes_simd = Simd::<u8, CHUNK_SIZE>::from_array(byte_chunk);

		// Add 33 to each byte to convert to ASCII
		let adjusted_simd = bytes_simd + Simd::splat(33);

		for j in 0..to_process
		{
			qual.push(adjusted_simd[j]);
		}

		// Update the index
		i += to_process;
	}

	qual.into_boxed_slice()
}

// fn read_byte_array_simd<T, F>(
// 	bytes: &[u8],
// 	offset: &mut usize,
// 	array_len: usize,
// 	convert: F,
// 	to_bam_val: fn(T) -> TagValueType,
// ) -> Vec<TagValueType>
// where
// 	F: Fn(&[u8]) -> T,
// 	T: Copy + Default + std::simd::SimdElement, // `T` must implement `Copy` and `Default` for the SIMD logic
// {
// 	const CHUNK_SIZE: usize = 8;

// 	let mut result = Vec::with_capacity(array_len);
// 	let mut i = 0;

// 	// SIMD processing
// 	while i + CHUNK_SIZE * std::mem::size_of::<T>() <= array_len * std::mem::size_of::<T>()
// 	{
// 		let mut chunk_values = [T::default(); CHUNK_SIZE];

// 		for j in 0..CHUNK_SIZE
// 		{
// 			let byte_index = i + j * std::mem::size_of::<T>();
// 			chunk_values[j] = convert(&bytes[byte_index..byte_index + std::mem::size_of::<T>()]);
// 		}

// 		let simd_chunk = Simd::from_array(chunk_values); // Creates a SIMD vector from the array
// 		for value in simd_chunk.to_array().iter()
// 		{
// 			result.push(to_bam_val(*value)); // Convert and push the result
// 		}

// 		i += CHUNK_SIZE * std::mem::size_of::<T>(); // Move forward by chunk size
// 	}

// 	// Scalar processing for remaining elements
// 	while i + std::mem::size_of::<T>() <= array_len * std::mem::size_of::<T>()
// 	{
// 		let chunk = &bytes[i..i + std::mem::size_of::<T>()];
// 		result.push(to_bam_val(convert(chunk)));
// 		i += std::mem::size_of::<T>();
// 	}

// 	*offset += i;
// 	result
// }
