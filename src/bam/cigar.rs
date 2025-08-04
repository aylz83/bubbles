use std::simd::Simd;
use std::collections::BTreeSet;

#[derive(Debug, Clone)]
pub enum Cigar
{
	Match(u32),
	Deletion(u32),
	Insertion(u32),
	Softclip(u32),
	Unknown,
}

impl Cigar
{
	fn from(opcode: u8, length: u32) -> Self
	{
		match opcode
		{
			b'M' | b'=' | b'X' => Cigar::Match(length),
			b'D' => Cigar::Deletion(length),
			b'I' => Cigar::Insertion(length),
			b'S' => Cigar::Softclip(length),
			_ => Cigar::Unknown,
		}
	}
}

const CIGAR_OPS: [u8; 9] = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];

const CIGAR_INDEX_LOOKUP: [bool; 256] = {
	let mut t = [false; 256];
	t[b'M' as usize] = true;
	t[b'=' as usize] = true;
	t[b'X' as usize] = true;
	t[b'D' as usize] = true;
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
