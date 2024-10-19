#![allow(dead_code)]

use tokio::fs::File as TokioFile;
use tokio::io::{AsyncReadExt, BufReader as TokioBufReader};
use anyhow::bail;

use std::path::Path;

use log::debug;

use async_compression::tokio::bufread::GzipDecoder;
use rustc_hash::FxHashMap;

use std::simd::{Simd, Mask};
use std::simd::cmp::SimdPartialEq;

const CIGAR_OPS: [char; 9] = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];

const ALPHABET: [char; 16] = [
	'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
];

use crate::bam;

pub struct Reader
{
	reader: TokioBufReader<TokioFile>,
}

impl Reader
{
	pub async fn from_path(path: &Path) -> anyhow::Result<Reader>
	{
		let file = TokioFile::open(path).await?;

		let reader = TokioBufReader::new(file);
		Ok(Reader { reader })
	}

	pub async fn read_header(&mut self) -> anyhow::Result<bam::header::Header>
	{
		let bytes = match Reader::read_bgzf_block(&mut self.reader).await?
		{
			Some(bytes) => bytes,
			None => bail!("BAM EOF before BAM header read?! Invalid BAM perhaps?"),
		};

		if !Reader::is_valid_bam(&bytes)
		{
			bail!("Input is not a valid BAM");
		}

		// obtain header text length
		let l_text = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize + 8;

		let n_ref = u32::from_le_bytes([
			bytes[l_text],
			bytes[l_text + 1],
			bytes[l_text + 2],
			bytes[l_text + 3],
		]);

		debug!("header text (l_text: {}) = {:?}", l_text, unsafe {
			std::str::from_utf8_unchecked(&bytes[8..l_text]).to_string()
		});

		debug!("n_ref: {}", n_ref);

		let mut offset = l_text + 3;

		let mut references = Vec::<bam::header::TID>::new();

		for _ in 1..n_ref
		{
			let l_name = u32::from_le_bytes([
				bytes[offset + 1],
				bytes[offset + 2],
				bytes[offset + 3],
				bytes[offset + 4],
			]) as usize;

			offset += 4;

			let name = &bytes[offset + 1..(offset + 1 + l_name)];

			offset += l_name;

			let l_ref = u32::from_le_bytes([
				bytes[offset + 1],
				bytes[offset + 2],
				bytes[offset + 3],
				bytes[offset + 4],
			]);

			offset += 4;

			references.push(bam::header::TID {
				name: unsafe { std::str::from_utf8_unchecked(name).to_string() },
				length: l_ref,
			});

			debug!(
				"l_name = {}, name = {}, l_ref = {}",
				l_name,
				unsafe { std::str::from_utf8_unchecked(name).to_string() },
				l_ref
			);
		}

		Ok(bam::header::Header {
			header: unsafe { std::str::from_utf8_unchecked(&bytes[8..l_text]).to_string() },
			references,
		})
	}

	pub async fn read_blocks<F>(
		&mut self,
		mut read_fn: F,
	) -> anyhow::Result<Vec<bam::pileup::Pileup>>
	where
		F: FnMut(&bam::Field),
	{
		let mut pileup = FxHashMap::<(i32, i32), u64>::default();

		loop
		{
			match self.read_block(&mut pileup, &mut read_fn).await?
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

		let mut pileup: Vec<_> = pileup
			.into_iter()
			.map(|((tid, pos), score)| bam::pileup::Pileup { tid, pos, score })
			.collect();
		pileup.sort_by(|a, b| (a.tid, a.pos).cmp(&(b.tid, b.pos)));

		Ok(pileup)
	}

	async fn read_block<F>(
		&mut self,
		pileup_map: &mut FxHashMap<(i32, i32), u64>,
		mut read_fn: F,
	) -> anyhow::Result<Option<()>>
	// Result due to some functions have the possibility to fail, and Option so we can detect None = EOF
	where
		F: FnMut(&bam::Field),
	{
		let bytes = match Reader::read_bgzf_block(&mut self.reader).await
		{
			Ok(bytes) => match bytes
			{
				Some(bytes) => bytes,
				None => return Ok(None),
			},
			Err(err) => bail!(err),
		};

		let mut offset: usize = 0;
		let mut start_block_offset: usize;
		while offset < bytes.len()
		{
			start_block_offset = offset;
			debug!("offset = {}", offset);
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
			let read_name = &bytes[offset + 1..offset + 1 + l_read_name as usize];
			debug!("read_name = {}", unsafe {
				std::str::from_utf8_unchecked(read_name).to_string()
			});

			offset += 1 + l_read_name as usize;

			// cigar - uint32_t[n_cigar_op]
			let mut ref_index = pos;
			let cigar = Self::process_cigar(
				&bytes,
				&mut offset,
				n_cigar_op as usize,
				ref_id,
				&mut ref_index,
				pileup_map,
			);

			debug!("cigar = {:?}", cigar);

			// seq - uint8_t[(l_seq + 1) / 2]
			let seq =
				Self::process_sequence(&bytes[offset..(offset + l_seq as usize)], l_seq as usize);

			debug!("seq = {:?}", &seq);

			offset += (l_seq as usize + 1) / 2;

			// qual - char[l_seq]
			let qual = Self::process_quality_scores(
				&bytes[offset..offset + l_seq as usize],
				l_seq as usize,
			);

			debug!("qual = {:?}", qual);

			offset += l_seq as usize;

			let mut tags_vec = Vec::<bam::Tag>::with_capacity(10);
			Reader::read_tags(
				&bytes,
				start_block_offset + block_size as usize,
				&mut offset,
				&mut tags_vec,
			)?;

			let read = bam::Field {
				ref_id,
				pos,
				mapq,
				bin,
				flags: flag,
				next_ref_id,
				next_pos,
				tlen,
				read_name: unsafe { std::str::from_utf8_unchecked(read_name).to_string() },
				sequence: seq.as_bytes(),
				sequence_quality: qual.as_bytes(),
				cigar: &cigar,
				tags: &tags_vec,
			};

			read_fn(&read);
		}

		Ok(Some(()))
	}

	fn process_cigar(
		bytes: &[u8],
		offset: &mut usize,
		n_cigar_op: usize,
		ref_id: i32,
		ref_index: &mut i32,
		pileup_map: &mut FxHashMap<(i32, i32), u64>,
	) -> Vec<bam::Cigar>
	{
		let mut cigar = Vec::<bam::Cigar>::with_capacity(n_cigar_op as usize);
		let chunk_size = 4; // Each CIGAR operation is 4 bytes
		let simd_width = 4; // Using SIMD with 4 `u32` at a time

		while cigar.len() < n_cigar_op
		{
			// Calculate remaining operations
			let remaining_ops = n_cigar_op - cigar.len();

			if remaining_ops >= simd_width && *offset + (simd_width * chunk_size) <= bytes.len()
			{
				// Process 4 CIGAR operations with SIMD
				let cigar_chunk = Simd::<u32, 4>::from_array([
					u32::from_le_bytes([
						bytes[*offset],
						bytes[*offset + 1],
						bytes[*offset + 2],
						bytes[*offset + 3],
					]),
					u32::from_le_bytes([
						bytes[*offset + 4],
						bytes[*offset + 5],
						bytes[*offset + 6],
						bytes[*offset + 7],
					]),
					u32::from_le_bytes([
						bytes[*offset + 8],
						bytes[*offset + 9],
						bytes[*offset + 10],
						bytes[*offset + 11],
					]),
					u32::from_le_bytes([
						bytes[*offset + 12],
						bytes[*offset + 13],
						bytes[*offset + 14],
						bytes[*offset + 15],
					]),
				]);

				// Extract operation codes and lengths in parallel
				let opcodes = cigar_chunk & Simd::splat(0xF); // Last 4 bits for opcode
				let lengths = cigar_chunk >> 4; // Remaining bits for length

				for i in 0..simd_width
				{
					let op = CIGAR_OPS[opcodes[i] as usize]; // Extract opcode
					let length = lengths[i]; // Extract length

					match op
					{
						'M' | '=' | 'X' =>
						{
							for j in 0..length as i32
							{
								let ref_pos = *ref_index + j;
								*pileup_map.entry((ref_id, ref_pos)).or_insert(0) += 1;
							}
							*ref_index += length as i32;
						}
						'D' =>
						{
							*ref_index += length as i32;
						}
						_ =>
						{}
					}

					cigar.push(bam::Cigar { length, opcode: op });
				}

				*offset += simd_width * chunk_size; // Move to the next set of CIGAR operations
			}
			else
			{
				// Fallback to scalar processing for remaining CIGAR operations
				while cigar.len() < n_cigar_op && *offset + chunk_size <= bytes.len()
				{
					let cigar_enc = u32::from_le_bytes([
						bytes[*offset],
						bytes[*offset + 1],
						bytes[*offset + 2],
						bytes[*offset + 3],
					]);
					let op = CIGAR_OPS[(cigar_enc & 0xF) as usize]; // Extract operation code
					let length = cigar_enc >> 4; // Extract length

					match op
					{
						'M' | '=' | 'X' =>
						{
							for j in 0..length as i32
							{
								let ref_pos = *ref_index + j;
								*pileup_map.entry((ref_id, ref_pos)).or_insert(0) += 1;
							}
							*ref_index += length as i32;
						}
						'D' =>
						{
							*ref_index += length as i32;
						}
						_ =>
						{}
					}

					cigar.push(bam::Cigar { length, opcode: op });
					*offset += chunk_size; // Move to the next CIGAR operation
				}
			}
		}

		cigar
	}

	fn process_sequence(bytes: &[u8], l_seq: usize) -> String
	{
		let mut seq = String::with_capacity(l_seq); // Allocate enough space for the sequence

		const CHUNK_SIZE: usize = 16; // SIMD width for u8

		let mut offset = 0; // Initialize the byte offset

		// Process the sequence in chunks
		while offset < (l_seq + 1) / 2
		{
			let remaining = (l_seq + 1) / 2 - offset; // Remaining bytes to process
			let to_process = remaining.min(CHUNK_SIZE);

			// Load bytes into a SIMD register
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

			// Update the sequence with the new bases
			seq.push_str(&bases.iter().collect::<String>());

			// Update the offset
			offset += to_process;
		}

		// Trim the sequence to the exact length if necessary
		seq.truncate(l_seq);

		seq
	}

	fn process_quality_scores(bytes: &[u8], l_seq: usize) -> String
	{
		let mut qual = String::with_capacity(l_seq); // Pre-allocate the string capacity

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

			// Create a SIMD register from the byte chunk
			let bytes_simd = Simd::<u8, CHUNK_SIZE>::from_array(byte_chunk);

			// Add 33 to each byte to convert to ASCII
			let adjusted_simd = bytes_simd + Simd::splat(33);

			// Convert to characters and append to the string
			for j in 0..to_process
			{
				qual.push(adjusted_simd[j] as char);
			}

			// Update the index
			i += to_process;
		}

		qual
	}

	fn read_tags(
		bytes: &[u8],
		block_size: usize,
		offset: &mut usize,
		tags: &mut Vec<bam::Tag>,
	) -> anyhow::Result<()>
	{
		while *offset <= block_size
		{
			// tag (char[2])
			// val_type (char)
			// value
			let tag =
				unsafe { std::str::from_utf8_unchecked(&bytes[*offset..*offset + 2]).to_string() };

			*offset += 2;

			//println!("tag = {}", tag);
			// Type (1 byte)
			let val_type = bytes[*offset] as char;
			*offset += 1;

			//println!("val_type = {}", val_type);
			// Parse value based on type
			let value = match val_type
			{
				'A' => Self::read_value(
					bytes,
					offset,
					|chunk| chunk[0] as char,
					bam::TagValueType::Char,
				),
				'c' =>
				{
					Self::read_value(bytes, offset, |chunk| chunk[0] as i8, bam::TagValueType::I8)
				}
				'C' => Self::read_value(bytes, offset, |chunk| chunk[0], bam::TagValueType::U8),
				's' => Self::read_value(
					bytes,
					offset,
					|chunk| i16::from_le_bytes([chunk[0], chunk[1]]),
					bam::TagValueType::I16,
				),
				'S' => Self::read_value(
					bytes,
					offset,
					|chunk| u16::from_le_bytes([chunk[0], chunk[1]]),
					bam::TagValueType::U16,
				),
				'i' => Self::read_value(
					bytes,
					offset,
					|chunk| i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
					bam::TagValueType::I32,
				),
				'I' => Self::read_value(
					bytes,
					offset,
					|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
					bam::TagValueType::U32,
				),
				'f' => Self::read_value(
					bytes,
					offset,
					|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
					bam::TagValueType::F32,
				),
				'Z' =>
				{
					match Self::find_null_terminator_simd(&bytes[*offset..], offset)
					{
						Some(null_offset) =>
						{
							let string = &bytes[*offset..*offset + null_offset];
							*offset += null_offset + 1; // Skip null terminator
							vec![bam::TagValueType::String(unsafe {
								std::str::from_utf8_unchecked(string).to_string()
							})]
						}
						None =>
						{
							vec![bam::TagValueType::String(String::new())]
						}
					}
					//match bytes[*offset..].iter().position(|&b| b == 0)
					//{
					//	Some(null_offset) =>
					//	{
					//		let string = &bytes[*offset..*offset + null_offset];
					//		*offset += null_offset + 1; // Skip null terminator
					//		vec![bam::TagValueType::String(unsafe {
					//			std::str::from_utf8_unchecked(string).to_string()
					//		})]
					//	}
					//	None =>
					//	{
					//		vec![bam::TagValueType::String(String::new())]
					//	}
					//}
				}
				'B' =>
				{
					// Byte array (Array of typed values)
					let array_type = bytes[*offset] as char;
					*offset += 1;
					let array_len = u32::from_le_bytes([
						bytes[*offset],
						bytes[*offset + 1],
						bytes[*offset + 2],
						bytes[*offset + 3],
					]) as usize;
					*offset += 4;

					match array_type
					{
						'c' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| i8::from_le_bytes([chunk[0]]),
							bam::TagValueType::I8,
						),
						'C' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| u8::from_le_bytes([chunk[0]]),
							bam::TagValueType::U8,
						),
						's' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| i16::from_le_bytes([chunk[0], chunk[1]]),
							bam::TagValueType::I16,
						),
						'S' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| u16::from_le_bytes([chunk[0], chunk[1]]),
							bam::TagValueType::U16,
						),
						'i' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
							bam::TagValueType::I32,
						),
						'I' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
							bam::TagValueType::U32,
						),
						'f' => Self::read_byte_array(
							bytes,
							offset,
							array_len,
							|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
							bam::TagValueType::F32,
						),
						_ => bail!("Unsupported array type: {}", array_type),
					}
				}

				_ => bail!("Unsupported tag type: {}", val_type),
			};

			// 10-15 seconds slower in some cases when adding to tags
			// find more efficient way, or add option to not parse tags
			tags.push(bam::Tag {
				name: tag,
				val_type,
				value,
			});

			//debug!("tag = {}:{}:{:?}", &tag, val_type, &value);
		}

		Ok(())
	}

	fn read_value<T, F>(
		bytes: &[u8],
		offset: &mut usize,
		convert: F,
		to_bam_val: fn(T) -> bam::TagValueType,
	) -> Vec<bam::TagValueType>
	where
		F: Fn(&[u8]) -> T,
	{
		let value = convert(&bytes[*offset..*offset + std::mem::size_of::<T>()]);
		*offset += std::mem::size_of::<T>();
		vec![to_bam_val(value)]
	}

	fn read_byte_array<T, F>(
		bytes: &[u8],
		offset: &mut usize,
		array_len: usize,
		convert: F,
		to_bam_val: fn(T) -> bam::TagValueType,
	) -> Vec<bam::TagValueType>
	where
		F: Fn(&[u8]) -> T,
	{
		let array = bytes[*offset..*offset + (array_len * std::mem::size_of::<T>())]
			.chunks_exact(std::mem::size_of::<T>())
			.map(convert)
			.map(to_bam_val)
			.collect::<Vec<bam::TagValueType>>();

		*offset += array_len * std::mem::size_of::<T>();
		array
	}

	fn read_byte_array_simd<T, F>(
		bytes: &[u8],
		offset: &mut usize,
		array_len: usize,
		convert: F,
		to_bam_val: fn(T) -> bam::TagValueType,
	) -> Vec<bam::TagValueType>
	where
		F: Fn(&[u8]) -> T,
		T: Copy + Default + std::simd::SimdElement, // `T` must implement `Copy` and `Default` for the SIMD logic
	{
		const CHUNK_SIZE: usize = 8; // Simd<T, 8> processes 8 elements at once (adjust size for T)

		let mut result = Vec::with_capacity(array_len);
		let mut i = 0;

		// SIMD processing
		while i + CHUNK_SIZE * std::mem::size_of::<T>() <= array_len * std::mem::size_of::<T>()
		{
			let mut chunk_values = [T::default(); CHUNK_SIZE]; // Create an array of default `T`

			for j in 0..CHUNK_SIZE
			{
				let byte_index = i + j * std::mem::size_of::<T>();
				chunk_values[j] =
					convert(&bytes[byte_index..byte_index + std::mem::size_of::<T>()]);
			}

			// Use SIMD to process these values
			let simd_chunk = Simd::from_array(chunk_values); // Creates a SIMD vector from the array
			for value in simd_chunk.to_array().iter()
			{
				result.push(to_bam_val(*value)); // Convert and push the result
			}

			i += CHUNK_SIZE * std::mem::size_of::<T>(); // Move forward by chunk size
		}

		// Scalar processing for remaining elements
		while i + std::mem::size_of::<T>() <= array_len * std::mem::size_of::<T>()
		{
			let chunk = &bytes[i..i + std::mem::size_of::<T>()];
			result.push(to_bam_val(convert(chunk)));
			i += std::mem::size_of::<T>();
		}

		*offset += i;
		result
	}

	fn find_null_terminator_simd(bytes: &[u8], offset: &mut usize) -> Option<usize>
	{
		const CHUNK_SIZE: usize = 16; // SIMD width: process 16 bytes at a time

		let mut i = 0;

		// SIMD processing loop
		while i + CHUNK_SIZE <= bytes.len()
		{
			// Load 16 bytes at once using SIMD
			let chunk = Simd::from_slice(&bytes[i..i + CHUNK_SIZE]);

			// Compare against zero (null terminator)
			let zeroes: Mask<_, 16> = chunk.simd_eq(Simd::splat(0));

			// Check if any of the elements are zero using Mask::any
			if zeroes.any()
			{
				// Convert the mask to a bitmask, then find the first `true` lane
				let bitmask = zeroes.to_bitmask();
				for j in 0..CHUNK_SIZE
				{
					if bitmask & (1 << j) != 0
					{
						return Some(i + j); // Return the position of the null byte
					}
				}
			}

			i += CHUNK_SIZE;
		}

		// Scalar processing for any remaining bytes
		while i < bytes.len()
		{
			if bytes[i] == 0
			{
				return Some(i);
			}
			i += 1;
		}

		None // Null terminator not found
	}

	async fn read_bgzf_block(
		reader: &mut TokioBufReader<TokioFile>,
	) -> anyhow::Result<Option<Vec<u8>>>
	{
		let mut header = [0; 18];

		match reader.read_exact(&mut header).await
		{
			Ok(_) =>
			{
				if !Reader::is_valid_bgzf_header(&header)
				{
					bail!("Invalid BGZF header"); // Skip to the next block if header is invalid.
				}

				// Calculate the size of the compressed block using BSIZE.
				let bsize = u16::from_le_bytes([header[16], header[17]]) as usize + 1;

				//println!("Header block size {} with contents: {:?}", bsize, header);

				// Read the rest of the BGZF block (bsize - 18 bytes).
				let mut compressed_block = vec![0; bsize];
				compressed_block[..18].copy_from_slice(&header);
				reader.read_exact(&mut compressed_block[18..]).await?;

				if bsize == 28 && Reader::is_bam_eof(&compressed_block)
				{
					debug!("EOF header = {:?}", compressed_block);
					return Ok(None);
				}

				let decompressed_block = Reader::decompress_block(&compressed_block).await?;
				Ok(Some(decompressed_block))
			}
			Err(e) =>
			{
				bail!("Failed to read BGZF header: {:?}", e);
			}
		}
	}

	fn is_bam_eof(bytes: &[u8]) -> bool
	{
		bytes[16..=27]
			== [
				0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
			]
	}

	fn is_valid_bam(bytes: &[u8]) -> bool
	{
		// check for magic BAM string ('BAM\1')
		bytes[0] == b'B' && bytes[1] == b'A' && bytes[2] == b'M' && bytes[3] == 1
	}

	fn is_valid_bgzf_header(header: &[u8]) -> bool
	{
		if header.len() != 18
		{
			return false; // Invalid header size.
		}

		// Check the fixed values in the header (no need for endianness consideration here).
		if header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08
		{
			return false; // Not a valid GZIP header.
		}

		// Check the subfield identifiers and length.
		if header[10] != 0x06 || header[11] != 0x00 || // XLEN = 6
		   header[12] != 0x42 || header[13] != 0x43 || // SI1 = 'B', SI2 = 'C'
	       header[14] != 0x02 || header[15] != 0x00
		{
			// SLEN = 2
			return false; // Not a valid BGZF header.
		}

		// Interpret BSIZE as a little-endian 16-bit integer.
		let bsize = u16::from_le_bytes([header[16], header[17]]);
		if bsize < 18
		{
			return false; // BSIZE should be at least 18 for a valid BGZF block.
		}

		true
	}

	async fn decompress_block(compressed_block: &[u8]) -> anyhow::Result<Vec<u8>>
	{
		let mut bytes: Vec<u8> = Vec::new();
		let mut decoder = GzipDecoder::new(compressed_block); // Skip the header for decompression
		decoder.read_to_end(&mut bytes).await?; // Collect the decompressed bytes into a Vec<u8>
		Ok(bytes)
	}
}
