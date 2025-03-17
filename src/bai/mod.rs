#![allow(dead_code)]

mod reader;

use std::ops::Range;
use std::collections::HashMap;

pub(crate) use crate::bai::reader::*;

pub(crate) struct Header
{
	pub(crate) n_refs: u32,
}

pub(crate) type Region = Vec<Range<u64>>;

pub(crate) struct Reference
{
	pub(crate) bins: HashMap<u32, Region>,
}

impl Reference
{
	pub(crate) fn flatten_bins(&self, search_bins: &[u32]) -> Vec<Range<u64>>
	{
		self.bins
			.iter()
			.filter(|(bins, _)| search_bins.contains(*bins))
			.flat_map(|(_, region)| region.clone())
			.collect()
	}

	pub(crate) fn flatten_all_bins(&self) -> Vec<Range<u64>>
	{
		self.bins
			.iter()
			.flat_map(|(_, region)| region.clone())
			.collect()
	}
}

pub(crate) fn region_to_bin(beg: u32, end: u32) -> u32
{
	let end = end - 1;
	if beg >> 14 == end >> 14
	{
		return ((1 << 15) - 1) / 7 + (beg >> 14);
	}

	if beg >> 17 == end >> 17
	{
		return ((1 << 12) - 1) / 7 + (beg >> 17);
	}

	if beg >> 20 == end >> 20
	{
		return ((1 << 9) - 1) / 7 + (beg >> 20);
	}

	if beg >> 23 == end >> 23
	{
		return ((1 << 6) - 1) / 7 + (beg >> 23);
	}

	if beg >> 26 == end >> 26
	{
		return ((1 << 3) - 1) / 7 + (beg >> 26);
	}

	0
}

pub(crate) fn region_to_bins(start: u64, stop: u64) -> Vec<u32>
{
	let mut bins = Vec::new();
	let stop = stop - 1;

	bins.push(0);

	for k in (1 + (start >> 26))..=(1 + (stop >> 26))
	{
		bins.push(k as u32);
	}

	for k in (9 + (start >> 23))..=(9 + (stop >> 23))
	{
		bins.push(k as u32);
	}

	for k in (73 + (start >> 20))..=(73 + (stop >> 20))
	{
		bins.push(k as u32);
	}

	for k in (585 + (start >> 17))..=(585 + (stop >> 17))
	{
		bins.push(k as u32);
	}

	for k in (4681 + (start >> 14))..=(4681 + (stop >> 14))
	{
		bins.push(k as u32);
	}

	bins.sort_by(|a, b| b.cmp(a));
	bins
}
