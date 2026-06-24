use crate::bam::Cigar;

#[derive(Debug, Clone, Copy)]
pub struct SpliceSegment
{
	pub read_start: u32,
	pub read_end: u32,
	pub ref_start: u32,
	pub ref_end: u32,
}

pub struct SplicedSegments<'a>
{
	pub(crate) ops: &'a [Cigar],
	pub(crate) idx: usize,
	pub(crate) read_pos: u32,
	pub(crate) ref_pos: u32,
}

impl<'a> Iterator for SplicedSegments<'a>
{
	type Item = SpliceSegment;

	fn next(&mut self) -> Option<Self::Item>
	{
		let mut started = false;
		let mut seg = None;

		while let Some(op) = self.ops.get(self.idx)
		{
			self.idx += 1;

			match *op
			{
				Cigar::Match(len, _) =>
				{
					if !started
					{
						seg = Some((self.read_pos, self.ref_pos));
						started = true;
					}

					self.read_pos += len;
					self.ref_pos += len;
				}

				Cigar::Insertion(len) =>
				{
					self.read_pos += len;
				}

				Cigar::Deletion(len) =>
				{
					self.ref_pos += len;
				}

				Cigar::Softclip(len) =>
				{
					self.read_pos += len;
				}

				Cigar::Skip(len) =>
				{
					// 🚨 splice boundary
					self.ref_pos += len;

					if let Some((rs, rfs)) = seg.take()
					{
						return Some(SpliceSegment {
							read_start: rs,
							read_end: self.read_pos,
							ref_start: rfs,
							ref_end: self.ref_pos - len,
						});
					}
				}

				Cigar::Hardclip(_) | Cigar::Pad(_) =>
				{}
				Cigar::Unknown =>
				{}
			}
		}

		// flush last segment
		seg.map(|(rs, rfs)| SpliceSegment {
			read_start: rs,
			read_end: self.read_pos,
			ref_start: rfs,
			ref_end: self.ref_pos,
		})
	}
}
