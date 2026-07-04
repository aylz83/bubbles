use std::simd::Simd;
use std::collections::BTreeSet;

use std::fmt;

use crate::bam::blocks::AlignedBlocks;
use crate::splicing::SplicedSegments;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchKind
{
	Match = b'=' as u8,
	Mismatch = b'X' as u8,
	Legacy = b'M' as u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Cigar
{
	Match(u32, MatchKind), // M, = or X
	Insertion(u32),        // I
	Deletion(u32),         // D
	Softclip(u32),         // S
	Hardclip(u32),         // H
	Skip(u32),             // N
	Pad(u32),              // P
	Unknown,
}

#[derive(Debug, Clone, Default)]
pub struct CigarString
{
	pub(crate) ops: Vec<Cigar>,
}

impl Cigar
{
	fn from(opcode: u8, length: u32) -> Self
	{
		match opcode
		{
			b'=' => Cigar::Match(length, MatchKind::Match),
			b'X' => Cigar::Match(length, MatchKind::Mismatch),
			b'M' => Cigar::Match(length, MatchKind::Legacy),
			b'I' => Cigar::Insertion(length),
			b'D' => Cigar::Deletion(length),
			b'S' => Cigar::Softclip(length),
			b'H' => Cigar::Hardclip(length),
			b'N' => Cigar::Skip(length),
			b'P' => Cigar::Pad(length),
			_ => Cigar::Unknown,
		}
	}

	pub fn len(&self) -> u32
	{
		match *self
		{
			Cigar::Match(l, _)
			| Cigar::Skip(l)
			| Cigar::Pad(l)
			| Cigar::Deletion(l)
			| Cigar::Insertion(l)
			| Cigar::Hardclip(l)
			| Cigar::Softclip(l) => l,
			Cigar::Unknown => 0,
		}
	}
}

impl CigarString
{
	pub fn new(ops: Vec<Cigar>) -> Self
	{
		Self { ops }
	}

	pub fn read_len(&self) -> u32
	{
		self.ops
			.iter()
			.map(|op| match *op
			{
				Cigar::Match(l, _) | Cigar::Insertion(l) | Cigar::Softclip(l) => l,
				_ => 0,
			})
			.sum()
	}

	pub fn ref_len(&self) -> u32
	{
		self.ops
			.iter()
			.map(|op| match *op
			{
				Cigar::Match(l, _) | Cigar::Deletion(l) | Cigar::Skip(l) => l,
				_ => 0,
			})
			.sum()
	}

	pub fn total_len(&self) -> u32
	{
		self.ops.iter().map(Cigar::len).sum()
	}

	pub fn is_empty(&self) -> bool
	{
		self.ops.is_empty()
	}

	pub fn is_spliced(&self) -> bool
	{
		self.ops.iter().any(|op| matches!(op, Cigar::Skip(_)))
	}

	pub fn is_contiguous(&self) -> bool
	{
		!self.ops.iter().any(|op| matches!(op, Cigar::Skip(_)))
	}

	pub fn has_indels(&self) -> bool
	{
		self.ops
			.iter()
			.any(|op| matches!(op, Cigar::Insertion(_) | Cigar::Deletion(_)))
	}

	pub fn splice_count(&self) -> usize
	{
		self.ops
			.iter()
			.filter(|op| matches!(op, Cigar::Skip(_)))
			.count()
	}

	pub fn splice_lengths(&self) -> impl Iterator<Item = u32> + '_
	{
		self.ops.iter().filter_map(|op| {
			if let Cigar::Skip(l) = op
			{
				Some(*l)
			}
			else
			{
				None
			}
		})
	}

	pub fn infer_introns(&self) -> impl Iterator<Item = (u32, u32)> + '_
	{
		let mut refp = 0u32;

		self.ops.iter().filter_map(move |op| match *op
		{
			Cigar::Match(l, _) =>
			{
				refp += l;
				None
			}

			Cigar::Insertion(_) | Cigar::Softclip(_) => None,

			Cigar::Deletion(l) =>
			{
				refp += l;
				None
			}

			Cigar::Skip(l) =>
			{
				let j = (refp, refp + l);
				refp += l;
				Some(j)
			}

			_ => None,
		})
	}

	pub fn ref_coverage(&self) -> impl Iterator<Item = (u32, u32)> + '_
	{
		self.aligned_blocks().map(|b| (b.ref_start, b.ref_end))
	}

	pub fn leading_softclip(&self) -> u32
	{
		match self.ops.first()
		{
			Some(Cigar::Softclip(l)) => *l,
			_ => 0,
		}
	}

	pub fn trailing_softclip(&self) -> u32
	{
		match self.ops.last()
		{
			Some(Cigar::Softclip(l)) => *l,
			_ => 0,
		}
	}

	pub fn spliced_segments(&self) -> SplicedSegments<'_>
	{
		SplicedSegments {
			ops: &self.ops,
			idx: 0,
			read_pos: 0,
			ref_pos: 0,
		}
	}

	pub fn aligned_blocks(&self) -> AlignedBlocks<'_>
	{
		AlignedBlocks {
			ops: &self.ops,
			idx: 0,
			read_pos: 0,
			ref_pos: 0,
		}
	}

	pub fn ops(&self) -> &[Cigar]
	{
		&self.ops
	}

	pub fn advance_ref(&self, start_ref: u32) -> u32
	{
		let delta: u32 = self
			.ops
			.iter()
			.map(|op| match *op
			{
				Cigar::Match(l, _) | Cigar::Deletion(l) | Cigar::Skip(l) => l,
				_ => 0,
			})
			.sum();

		start_ref + delta
	}

	pub fn read_to_ref(&self, read_coord: u32) -> Option<u32>
	{
		for block in self.aligned_blocks()
		{
			if read_coord >= block.read_start && read_coord < block.read_end
			{
				return Some(block.ref_start + (read_coord - block.read_start));
			}
		}
		None
	}

	pub fn ref_to_read(&self, ref_coord: u32) -> Option<u32>
	{
		for block in self.aligned_blocks()
		{
			if ref_coord >= block.ref_start && ref_coord < block.ref_end
			{
				return Some(block.read_start + (ref_coord - block.ref_start));
			}
		}
		None
	}
}

impl fmt::Display for CigarString
{
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result
	{
		for op in &self.ops
		{
			match *op
			{
				Cigar::Match(l, MatchKind::Legacy) => write!(f, "{}M", l)?,
				Cigar::Match(l, MatchKind::Match) => write!(f, "{}=", l)?,
				Cigar::Match(l, MatchKind::Mismatch) => write!(f, "{}X", l)?,
				Cigar::Skip(l) => write!(f, "{}N", l)?,
				Cigar::Pad(l) => write!(f, "{}P", l)?,
				Cigar::Deletion(l) => write!(f, "{}D", l)?,
				Cigar::Insertion(l) => write!(f, "{}I", l)?,
				Cigar::Softclip(l) => write!(f, "{}S", l)?,
				Cigar::Hardclip(l) => write!(f, "{}H", l)?,
				Cigar::Unknown => write!(f, "*")?,
			}
		}
		Ok(())
	}
}

const CIGAR_OPS: [u8; 9] = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];

const CIGAR_INDEX_LOOKUP: [bool; 256] = {
	let mut t = [false; 256];
	t[b'M' as usize] = true;
	t[b'=' as usize] = true;
	t[b'X' as usize] = true;
	t[b'D' as usize] = true;
	t[b'N' as usize] = true;
	t
};

pub(crate) fn process_cigar(
	bytes: &[u8],
	offset: &mut usize,
	n_cigar_op: usize,
	_ref_id: i32,
	ref_index: &mut i32,
	//pileup_map: &mut FxHashMap<(i32, i32), u64>,
	coverage_set: &mut Option<BTreeSet<u64>>,
) -> Vec<Cigar>
{
	let mut cigar = vec![Cigar::Unknown; n_cigar_op as usize];
	let chunk_size = 4; // Each CIGAR operation is 4 bytes
	let simd_width = 4; // Using SIMD with 4 `u32` at a time

	let mut cigar_pos = 0;

	while cigar_pos < n_cigar_op
	{
		// Calculate remaining operations
		let remaining_ops = n_cigar_op - cigar_pos;

		if remaining_ops >= simd_width && *offset + (simd_width * chunk_size) <= bytes.len()
		{
			// Process 4 CIGAR operations with SIMD
			let cigar_chunk = unsafe {
				core::slice::from_raw_parts(bytes.as_ptr().add(*offset) as *const u32, simd_width)
			};

			let cigar_chunk = Simd::<u32, 4>::from_array([
				cigar_chunk[0].to_le(),
				cigar_chunk[1].to_le(),
				cigar_chunk[2].to_le(),
				cigar_chunk[3].to_le(),
			]);

			// Extract operation codes and lengths in parallel
			let opcodes = cigar_chunk & Simd::splat(0xF); // Last 4 bits for opcode
			let lengths = cigar_chunk >> 4; // Remaining bits for length

			for i in 0..simd_width
			{
				let op = CIGAR_OPS[opcodes[i] as usize];
				let length = lengths[i];

				if op == b'M' || op == b'=' || op == b'X'
				{
					if let Some(coverage) = coverage_set
					{
						let mut cov_buf = [0u64; 1024];
						let mut cov_count = 0;

						for j in 0..length
						{
							cov_buf[cov_count] = (*ref_index + j as i32) as u64;
							cov_count += 1;
						}

						coverage.extend(&cov_buf[..cov_count]);
					}
				}

				if CIGAR_INDEX_LOOKUP[op as usize]
				{
					*ref_index += length as i32;
				}

				cigar[cigar_pos] = Cigar::from(op, length);
				cigar_pos = cigar_pos + 1;
			}

			*offset += simd_width * chunk_size; // Move to the next set of CIGAR operations
		}
		else
		{
			// Fallback to scalar processing for remaining CIGAR operations
			while cigar_pos < n_cigar_op && *offset + chunk_size <= bytes.len()
			{
				let cigar_enc = u32::from_le_bytes([
					bytes[*offset],
					bytes[*offset + 1],
					bytes[*offset + 2],
					bytes[*offset + 3],
				]);

				let op = CIGAR_OPS[(cigar_enc & 0xF) as usize]; // Extract operation code
				let length = cigar_enc >> 4; // Extract length

				if op == b'M' || op == b'=' || op == b'X'
				{
					if let Some(coverage) = coverage_set
					{
						let mut cov_buf = [0u64; 1024];
						let mut cov_count = 0;

						for j in 0..length
						{
							cov_buf[cov_count] = (*ref_index + j as i32) as u64;
							cov_count += 1;
						}

						coverage.extend(&cov_buf[..cov_count]);
					}
				}

				if CIGAR_INDEX_LOOKUP[op as usize]
				{
					*ref_index += length as i32;
				}

				cigar[cigar_pos] = Cigar::from(op, length);
				cigar_pos = cigar_pos + 1;
				*offset += chunk_size; // Move to the next CIGAR operation
			}
		}
	}

	cigar
}
