use thiserror::Error;

use pufferfish::error::Error as PufferfishError;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Error, Debug)]
pub enum Error
{
	#[error("Unable to open file {0}")]
	IOError(String),
	#[error("Not in BAM format")]
	BamFormat,
	#[error("Not in BAI format")]
	BaiFormat,
	#[error("Unable to decompress BGZ block")]
	BGZDecompress,
	#[error("Unable to read BGZ block")]
	BGZRead,
	#[error("Invalid BGZ header: {0:?}")]
	BGZInvalidHeader([u8; 18]),
	#[error("No BAI index found")]
	NoIndex,
	#[error("Unsupported array type '{0}' in BAM")]
	BamArrayType(char),
	#[error("Unsupported array tag '{0}' in BAM")]
	BamArrayTag(char),
	#[error("Unable to seek to {0} in BAM")]
	BamSeek(u64),
	#[error("Unable to read BGZ block in BAM")]
	BamBlock,
	#[error(transparent)]
	Pufferfish(#[from] PufferfishError),
}
