#![allow(dead_code)]

use std::path::Path;
use std::ops::Range;
use std::collections::HashMap;

use tokio::fs::File as TokioFile;
use tokio::io::{AsyncReadExt, BufReader as TokioBufReader};

use log::debug;

use crate::bai::{Header, Region, Reference};
use crate::error;
use crate::AsyncReadSeek;

pub(crate) struct Reader
{
	pub(crate) ref_indices: Vec<Reference>,
}

impl Reader
{
	pub async fn from_reader<R>(reader: R) -> error::Result<Reader>
	where
		R: AsyncReadSeek + std::marker::Send + std::marker::Unpin,
	{
		let mut reader = TokioBufReader::new(reader);

		let mut bytes = [0u8; 4];
		reader
			.read_exact(&mut bytes)
			.await
			.map_err(|_| error::Error::BaiFormat)?;

		//debug!("{:?}", unsafe { std::str::from_utf8_unchecked(&bytes) });

		reader
			.read_exact(&mut bytes)
			.await
			.map_err(|_| error::Error::BaiFormat)?;

		let n_refs = u32::from_le_bytes(bytes);
		debug!("n_refs = {}", n_refs);

		let _header = Header { n_refs };

		let mut large_bytes = [0u8; 8];

		let mut ref_indices = Vec::with_capacity(n_refs as usize);

		for _ in 0..n_refs
		{
			//debug!("at_ref = {}", at_ref);

			reader
				.read_exact(&mut bytes)
				.await
				.map_err(|_| error::Error::BaiFormat)?;

			let n_bin = u32::from_le_bytes(bytes);
			//debug!("n_bin = {}", n_bin);

			let mut bins_map = HashMap::with_capacity(n_bin as usize);

			for _ in 0..n_bin
			{
				//debug!("at_bin = {}", at_bin);

				reader
					.read_exact(&mut bytes)
					.await
					.map_err(|_| error::Error::BaiFormat)?;

				let bin = u32::from_le_bytes(bytes);
				//debug!("bin = {}", n_bin);

				reader
					.read_exact(&mut bytes)
					.await
					.map_err(|_| error::Error::BaiFormat)?;

				let n_chunk = u32::from_le_bytes(bytes);
				//debug!("n_chunk = {}", n_chunk);

				let mut chunks = Region::with_capacity(n_chunk as usize);

				for _ in 0..n_chunk
				{
					reader
						.read_exact(&mut large_bytes)
						.await
						.map_err(|_| error::Error::BaiFormat)?;

					let chunk_beg = u64::from_le_bytes(large_bytes);
					//debug!("chunk_beg = {}", chunk_beg);

					reader
						.read_exact(&mut large_bytes)
						.await
						.map_err(|_| error::Error::BaiFormat)?;

					let chunk_end = u64::from_le_bytes(large_bytes);
					//debug!("chunk_end = {}", chunk_end);

					chunks.push(Range {
						start: chunk_beg,
						end: chunk_end,
					});
				}

				if bin != 37450
				// skip unmapped for now
				{
					bins_map.insert(bin, chunks);
				}
			}

			reader
				.read_exact(&mut bytes)
				.await
				.map_err(|_| error::Error::BaiFormat)?;

			let n_intv = u32::from_le_bytes(bytes);
			//debug!("n_intv = {}", n_intv);

			for _ in 0..n_intv
			{
				//debug!("at_intv = {}", at_intv);

				reader
					.read_exact(&mut large_bytes)
					.await
					.map_err(|_| error::Error::BaiFormat)?;

				let _ioffset = u64::from_le_bytes(large_bytes);
				//debug!("ioffset = {}", ioffset);
			}

			ref_indices.push(Reference { bins: bins_map });
		}

		Ok(Reader { ref_indices })
	}

	pub async fn from_path(path: &Path) -> error::Result<Reader>
	{
		let bai_file = TokioFile::open(path)
			.await
			.map_err(|_| error::Error::IOError(path.to_string_lossy().to_string()))?;

		Self::from_reader(bai_file).await
	}
}
