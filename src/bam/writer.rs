use std::sync::Arc;

use tokio::io::AsyncWrite;

use crate::error;
use pufferfish::bgzf::pool::BgzfBlockPool;
use pufferfish::bgzf::writer::{BgzfBlockWriter, EncodeBgzfPayload};

use super::{Field, Header};

pub struct AsyncBamWriter<W>
where
	W: AsyncWrite + Send + Unpin,
{
	writer: BgzfBlockWriter<W>,
}

impl<W> AsyncBamWriter<W>
where
	W: AsyncWrite + Send + Unpin,
{
	pub async fn new_with_header(
		writer: W,
		header: Header,
		pool: Arc<rayon::ThreadPool>,
		block_pool: Arc<BgzfBlockPool>,
	) -> error::Result<Self>
	{
		let writer = BgzfBlockWriter::new_with_header(writer, header, pool, block_pool).await?;

		Ok(Self { writer })
	}

	pub async fn write<I, T>(&self, records: I) -> error::Result<()>
	where
		I: IntoIterator<Item = T>,
		T: EncodeBgzfPayload<Error = error::Error>,
	{
		self.writer.write(records).await
	}

	pub async fn write_record(&self, record: Field) -> error::Result<()>
	{
		self.writer.write(std::iter::once(record)).await
	}

	pub async fn flush(&self) -> error::Result<()>
	{
		self.writer.flush().await?;
		Ok(())
	}

	pub async fn finish(&self) -> error::Result<()>
	{
		self.writer.finish().await?;
		Ok(())
	}
}

#[cfg(test)]
mod test
{
	#[cfg(test)]
	mod test
	{
		use std::sync::Arc;

		#[tokio::test]
		async fn test_write_bam() -> crate::error::Result<()>
		{
			let path = std::env::temp_dir().join("test_output.bam");

			let file = tokio::fs::File::create(&path).await?;

			let header = crate::bam::Header {
				header: b"@HD\tVN:1.6\tSO:coordinate\n".to_vec().into_boxed_slice(),
				references: vec![],
			};

			let pool = Arc::new(
				rayon::ThreadPoolBuilder::new()
					.num_threads(2)
					.build()
					.expect("failed to create thread pool"),
			);

			// Replace this with the actual BgzfBlockPool constructor.
			let block_pool = Arc::new(crate::pufferfish::pool::BgzfBlockPool::new(
				10000,
				64 * 1024,
			));

			let writer =
				crate::bam::AsyncBamWriter::new_with_header(file, header, pool, block_pool).await?;

			writer
				.write(std::iter::empty::<crate::bam::Field>())
				.await?;
			writer.finish().await?;

			let bytes = tokio::fs::read(&path).await?;

			assert!(!bytes.is_empty());
			assert!(bytes.ends_with(&pufferfish::bgzf::BGZF_EOF_BLOCK));

			tokio::fs::remove_file(path).await?;

			Ok(())
		}
	}
}
