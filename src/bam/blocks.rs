use crate::bam::Cigar;

#[derive(Debug, Clone, Copy)]
pub struct AlignedBlock
{
	pub read_start: u32,
	pub read_end: u32, // exclusive
	pub ref_start: u32,
	pub ref_end: u32, // exclusive
}

pub struct AlignedBlocks<'a>
{
	pub(crate) ops: &'a [Cigar],
	pub(crate) idx: usize,
	pub(crate) read_pos: u32,
	pub(crate) ref_pos: u32,
}

impl<'a> Iterator for AlignedBlocks<'a>
{
	type Item = AlignedBlock;

	fn next(&mut self) -> Option<Self::Item>
	{
		while let Some(op) = self.ops.get(self.idx)
		{
			self.idx += 1;

			match *op
			{
				Cigar::Match(len, _) =>
				{
					let block = AlignedBlock {
						read_start: self.read_pos,
						read_end: self.read_pos + len,
						ref_start: self.ref_pos,
						ref_end: self.ref_pos + len,
					};

					self.read_pos += len;
					self.ref_pos += len;

					return Some(block);
				}

				Cigar::Insertion(len) | Cigar::Softclip(len) =>
				{
					self.read_pos += len;
				}

				Cigar::Skip(len) =>
				{
					self.ref_pos += len;
				}

				Cigar::Deletion(len) =>
				{
					self.ref_pos += len;
				}

				Cigar::Hardclip(_) | Cigar::Pad(_) =>
				{}

				Cigar::Unknown =>
				{}
			}
		}

		None
	}
}
