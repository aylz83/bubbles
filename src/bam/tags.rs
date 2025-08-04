use log::debug;

use std::simd::{Simd, Mask};
use std::simd::cmp::SimdPartialEq;

use crate::error;

#[derive(Debug)]
pub enum TagValueType
{
	Char(char),
	I8(i8),
	U8(u8),
	I16(i16),
	U16(u16),
	I32(i32),
	U32(u32),
	F32(f32),
	String(Box<[u8]>),
}

pub struct Tag
{
	pub name: Box<[u8]>,
	pub val_type: char,
	pub value: Vec<TagValueType>,
}

impl Tag
{
	pub fn name_as_str(&self) -> &str
	{
		unsafe { std::str::from_utf8_unchecked(&self.name[..self.name.len() - 1]) }
	}
}

pub(crate) fn read_tags(
	bytes: &[u8],
	block_size: usize,
	offset: &mut usize,
	tags: &mut Vec<Tag>,
) -> error::Result<()>
{
	while *offset <= block_size
	{
		// tag (char[2])
		let tag = Box::from(&bytes[*offset..*offset + 2]);

		*offset += 2;

		// val_type (char)
		let val_type = bytes[*offset] as char;
		*offset += 1;

		// Parse value based on type
		let value = match val_type
		{
			'A' => read_value(
				bytes,
				offset,
				|chunk| chunk[0] as char,
				1,
				TagValueType::Char,
			),
			'c' => read_value(
				bytes,
				offset,
				|chunk| chunk[0] as i8,
				std::mem::size_of::<i8>(),
				TagValueType::I8,
			),
			'C' => read_value(
				bytes,
				offset,
				|chunk| chunk[0],
				std::mem::size_of::<u8>(),
				TagValueType::U8,
			),
			's' => read_value(
				bytes,
				offset,
				|chunk| i16::from_le_bytes([chunk[0], chunk[1]]),
				std::mem::size_of::<i16>(),
				TagValueType::I16,
			),
			'S' => read_value(
				bytes,
				offset,
				|chunk| u16::from_le_bytes([chunk[0], chunk[1]]),
				std::mem::size_of::<u16>(),
				TagValueType::U16,
			),
			'i' => read_value(
				bytes,
				offset,
				|chunk| i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
				std::mem::size_of::<i32>(),
				TagValueType::I32,
			),
			'I' => read_value(
				bytes,
				offset,
				|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
				std::mem::size_of::<u32>(),
				TagValueType::U32,
			),
			'f' => read_value(
				bytes,
				offset,
				|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
				std::mem::size_of::<f32>(),
				TagValueType::F32,
			),
			'Z' =>
			{
				match find_null_terminator_simd(&bytes[*offset..], offset)
				{
					Some(null_offset) =>
					{
						let string = &bytes[*offset..*offset + null_offset];
						*offset += null_offset + 1; // Skip null terminator
						vec![TagValueType::String(Box::from(string))]
					}
					None =>
					{
						vec![TagValueType::String(Box::from(&[][..]))]
					}
				}
				//match bytes[*offset..].iter().position(|&b| b == 0)
				//{
				//	Some(null_offset) =>
				//	{
				//		let string = &bytes[*offset..*offset + null_offset];
				//		*offset += null_offset + 1; // Skip null terminator
				//		vec![TagValueType::String(unsafe {
				//			std::str::from_utf8_unchecked(string).to_string()
				//		})]
				//	}
				//	None =>
				//	{
				//		vec![TagValueType::String(String::new())]
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
					'c' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| i8::from_le_bytes([chunk[0]]),
						TagValueType::I8,
					),
					'C' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| u8::from_le_bytes([chunk[0]]),
						TagValueType::U8,
					),
					's' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| i16::from_le_bytes([chunk[0], chunk[1]]),
						TagValueType::I16,
					),
					'S' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| u16::from_le_bytes([chunk[0], chunk[1]]),
						TagValueType::U16,
					),
					'i' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
						TagValueType::I32,
					),
					'I' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
						TagValueType::U32,
					),
					'f' => read_byte_array(
						bytes,
						offset,
						array_len,
						|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
						TagValueType::F32,
					),
					_ => return Err(error::Error::BamArrayType(array_type)),
				}
			}

			_ =>
			{
				debug!("tag error = {:?}:{}", &tag, val_type);
				return Err(error::Error::BamArrayTag(val_type));
			}
		};

		debug!("tag = {:?}:{}:{:?}", &tag, val_type, &value);
		// 10-15 seconds slower in some cases when adding to tags
		// find more efficient way
		tags.push(Tag {
			name: tag,
			val_type,
			value,
		});
	}

	Ok(())
}

pub(crate) fn read_value<T, F>(
	bytes: &[u8],
	offset: &mut usize,
	convert: F,
	size: usize,
	to_bam_val: fn(T) -> TagValueType,
) -> Vec<TagValueType>
where
	F: Fn(&[u8]) -> T,
{
	let value = convert(&bytes[*offset..*offset + size]);
	*offset += size;
	vec![to_bam_val(value)]
}

pub(crate) fn read_byte_array<T, F>(
	bytes: &[u8],
	offset: &mut usize,
	array_len: usize,
	convert: F,
	to_bam_val: fn(T) -> TagValueType,
) -> Vec<TagValueType>
where
	F: Fn(&[u8]) -> T,
{
	let array = bytes[*offset..*offset + (array_len * std::mem::size_of::<T>())]
		.chunks_exact(std::mem::size_of::<T>())
		.map(convert)
		.map(to_bam_val)
		.collect::<Vec<TagValueType>>();

	*offset += array_len * std::mem::size_of::<T>();
	array
}

pub(crate) fn find_null_terminator_simd(bytes: &[u8], _offset: &mut usize) -> Option<usize>
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
