use thiserror::Error;

use pufferfish::error::Error as PufferfishError;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug)]
pub enum Error
{
	#[error("Not in BAM format")]
	BamFormat,
	#[error("Not in BAI format")]
	BaiFormat,
	#[error("No BAI index found")]
	NoIndex,
	#[error("Unsupported array type '{0}' in BAM")]
	BamArrayType(u8),
	#[error("Unsupported array tag '{0}' in BAM")]
	BamArrayTag(u8),
	#[error("Unable to seek to {0} in BAM")]
	BamSeek(u64),
	#[error("Unable to read BGZ block in BAM")]
	BamBlock,
	#[error(transparent)]
	Pufferfish(#[from] PufferfishError),
	#[error("IO error: {0}")]
	Io(#[from] std::io::Error),
}
