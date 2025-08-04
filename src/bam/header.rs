use crate::error;

use pufferfish::BGZ;

use tokio::io::BufReader as TokioBufReader;
use tokio::io::AsyncRead;

use log::debug;

#[derive(Debug)]
pub struct TID
{
	pub name: Box<[u8]>,
	pub length: u32,
}

impl TID
{
	pub fn name_as_str(&self) -> &str
	{
		unsafe { std::str::from_utf8_unchecked(&self.name[..self.name.len() - 1]) }
	}
}

#[derive(Debug)]
pub struct Header
{
	pub header: Box<[u8]>,
	pub references: Vec<TID>,
}

impl Header
{
	pub fn ref_name(&self, ref_id: i32) -> Option<&TID>
	{
		self.references.get(ref_id as usize)
	}

	pub fn header_as_str(&self) -> &str
	{
		unsafe { std::str::from_utf8_unchecked(&self.header[..self.header.len() - 1]) }
	}
}

pub(crate) async fn read_bam_header<R>(reader: &mut TokioBufReader<R>) -> error::Result<Header>
where
	R: AsyncRead + Send + std::marker::Unpin,
{
	let bytes = match reader.read_bgzf_block(Some(pufferfish::is_bam_eof)).await?
	{
		Some(bytes) => bytes,
		None => return Err(error::Error::BamFormat),
	};

	if !is_valid_bam(&bytes)
	{
		return Err(error::Error::BamFormat);
	}

	// obtain header text length
	let l_text = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize + 8;

	let n_ref = u32::from_le_bytes([
		bytes[l_text],
		bytes[l_text + 1],
		bytes[l_text + 2],
		bytes[l_text + 3],
	]);

	//debug!("header text (l_text: {}) = {:?}", l_text, unsafe {
	//	std::str::from_utf8_unchecked(&bytes[8..l_text]).to_string()
	//});

	debug!("n_ref: {}", n_ref);

	let mut offset = l_text + 3;

	let mut references = Vec::<TID>::new();

	for _ in 0..n_ref
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

		references.push(TID {
			name: Box::from(&name[0..l_name - 1]),
			length: l_ref,
		});

		//debug!(
		//	"l_name = {}, name = {}, l_ref = {}",
		//	l_name,
		//	unsafe { std::str::from_utf8_unchecked(name).to_string() },
		//	l_ref
		//);
	}

	Ok(Header {
		header: Box::from(&bytes[8..l_text]),
		references,
	})
}

pub(crate) fn is_valid_bam(bytes: &[u8]) -> bool
{
	// check for magic BAM string ('BAM\1')
	bytes[0] == b'B' && bytes[1] == b'A' && bytes[2] == b'M' && bytes[3] == 1
}
