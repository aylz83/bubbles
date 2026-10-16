use crate::error;

use pufferfish::prelude::*;

use tokio::io::{AsyncRead, AsyncSeek};

// use log::debug;

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
		unsafe { std::str::from_utf8_unchecked(&self.name[..self.name.len()]) }
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

pub(crate) async fn read_bam_header<R>(reader: &mut R) -> error::Result<Header>
where
	R: AsyncRead + AsyncSeek + Send + Unpin + BGZ,
{
	let mut bytes = Vec::new();

	loop
	{
		// Read enough to validate the magic and get l_text.
		while bytes.len() < 8
		{
			let block = reader
				.read_and_decompress_bgzf_block(Some(is_bam_eof))
				.await?
				.ok_or(error::Error::BamFormat)?;

			bytes.extend_from_slice(&block);
		}

		if !is_valid_bam(&bytes)
		{
			return Err(error::Error::BamFormat);
		}

		let l_text = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;

		let text_end = 8usize.checked_add(l_text).ok_or(error::Error::BamFormat)?;

		while bytes.len() < text_end + 4
		{
			let block = reader
				.read_and_decompress_bgzf_block(Some(is_bam_eof))
				.await?
				.ok_or(error::Error::BamFormat)?;

			bytes.extend_from_slice(&block);
		}

		let n_ref = u32::from_le_bytes(bytes[text_end..text_end + 4].try_into()?) as usize;

		let mut offset = text_end + 4;
		let mut complete = true;

		for _ in 0..n_ref
		{
			if bytes.len() < offset + 4
			{
				complete = false;
				break;
			}

			let l_name = u32::from_le_bytes(bytes[offset..offset + 4].try_into()?) as usize;

			let ref_end = offset
				.checked_add(4)
				.and_then(|v| v.checked_add(l_name))
				.and_then(|v| v.checked_add(4))
				.ok_or(error::Error::BamFormat)?;

			if bytes.len() < ref_end
			{
				complete = false;
				break;
			}

			offset = ref_end;
		}

		if complete
		{
			let mut references = Vec::<TID>::with_capacity(n_ref);
			offset = text_end + 4;

			for _ in 0..n_ref
			{
				let l_name = u32::from_le_bytes(bytes[offset..offset + 4].try_into()?) as usize;
				offset += 4;

				let name = &bytes[offset..offset + l_name - 1];
				offset += l_name;

				let l_ref = u32::from_le_bytes(bytes[offset..offset + 4].try_into()?);
				offset += 4;

				references.push(TID {
					name: Box::from(name),
					length: l_ref,
				});
			}

			return Ok(Header {
				header: Box::from(&bytes[8..text_end]),
				references,
			});
		}

		// The header spans more blocks. Read another one and retry.
		let block = reader
			.read_and_decompress_bgzf_block(Some(is_bam_eof))
			.await?
			.ok_or(error::Error::BamFormat)?;

		bytes.extend_from_slice(&block);
	}
}

// pub(crate) async fn read_bam_header<R>(reader: &mut TokioBufReader<R>) -> error::Result<Header>
// where
// 	R: AsyncRead + AsyncSeek + Send + std::marker::Unpin,
// {
// 	debug!("Reading BAM header");

// 	// reader.seek(SeekFrom::Start(0)).await?;

// 	let mut bytes = Vec::new();

// 	while let Some(buffer) = reader
// 		.read_and_decompress_bgzf_block(Some(is_bam_eof))
// 		.await?
// 	{
// 		bytes.extend_from_slice(&buffer);
// 	}

// 	// let bytes = match reader
// 	// 	.read_and_decompress_bgzf_block(Some(is_bam_eof))
// 	// 	.await?
// 	// {
// 	// 	Some(bytes) => bytes,
// 	// 	None => return Err(error::Error::BamFormat),
// 	// };

// 	debug!("Read BAM header block");

// 	if !is_valid_bam(&bytes)
// 	{
// 		debug!("Not a valid BAM: {:?}", bytes);
// 		return Err(error::Error::BamFormat);
// 	}

// 	debug!("bytes = {:?}", bytes);

// 	// obtain header text length
// 	let l_text = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize + 8;

// 	debug!("l_text = {}", l_text);

// 	let n_ref = u32::from_le_bytes([
// 		bytes[l_text],
// 		bytes[l_text + 1],
// 		bytes[l_text + 2],
// 		bytes[l_text + 3],
// 	]);

// 	//debug!("header text (l_text: {}) = {:?}", l_text, unsafe {
// 	//	std::str::from_utf8_unchecked(&bytes[8..l_text]).to_string()
// 	//});

// 	debug!("n_ref: {}", n_ref);

// 	let mut offset = l_text + 3;

// 	let mut references = Vec::<TID>::new();

// 	for _ in 0..n_ref
// 	{
// 		let l_name = u32::from_le_bytes([
// 			bytes[offset + 1],
// 			bytes[offset + 2],
// 			bytes[offset + 3],
// 			bytes[offset + 4],
// 		]) as usize;

// 		offset += 4;

// 		let name = &bytes[offset + 1..(offset + 1 + l_name)];

// 		offset += l_name;

// 		let l_ref = u32::from_le_bytes([
// 			bytes[offset + 1],
// 			bytes[offset + 2],
// 			bytes[offset + 3],
// 			bytes[offset + 4],
// 		]);

// 		offset += 4;

// 		references.push(TID {
// 			name: Box::from(&name[0..l_name - 1]),
// 			length: l_ref,
// 		});

// 		//debug!(
// 		//	"l_name = {}, name = {}, l_ref = {}",
// 		//	l_name,
// 		//	unsafe { std::str::from_utf8_unchecked(name).to_string() },
// 		//	l_ref
// 		//);
// 	}

// 	Ok(Header {
// 		header: Box::from(&bytes[8..l_text]),
// 		references,
// 	})
// }

pub(crate) fn is_valid_bam(bytes: &[u8]) -> bool
{
	// check for magic BAM string ('BAM\1')
	bytes[0] == b'B' && bytes[1] == b'A' && bytes[2] == b'M' && bytes[3] == 1
}
