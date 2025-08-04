use hashbrown::hash_map::RawEntryMut;
use hashbrown::HashMap;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

use crate::bam::{Field, Cigar};

type FxRawMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

pub struct Pileup
{
	pub tid: i32,
	pub pos: i32,
	pub score: u64,
}

pub(crate) fn generate_pileup(reads: &Vec<Field>) -> FxRawMap<(i32, i32), u64>
{
	let mut pileup: FxRawMap<(i32, i32), u64> = FxRawMap::default();

	for read in reads.iter()
	{
		let mut ref_pos = read.pos;
		//let mut seq_pos = 0;

		if let Some(cigar_ops) = &read.cigar
		{
			for cigar in cigar_ops
			{
				match cigar
				{
					Cigar::Match(length) =>
					{
						for i in 0..*length
						{
							let key = (read.ref_id, ref_pos + i as i32);
							match pileup.raw_entry_mut().from_key(&key)
							{
								RawEntryMut::Occupied(mut entry) =>
								{
									*entry.get_mut() += 1;
								}
								RawEntryMut::Vacant(entry) =>
								{
									entry.insert(key, 1);
								}
							}
						}
						ref_pos += *length as i32;
					}
					Cigar::Deletion(length) =>
					{
						// Move reference position forward (gap in query)
						ref_pos += *length as i32;
					}
					_ =>
					{}
				}
			}
		}
	}

	pileup
}
