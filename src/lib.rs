#![feature(portable_simd)]

mod bai;
pub mod bam;
pub mod error;
pub mod splicing;

pub use pufferfish::prelude as pufferfish;

pub mod prelude
{
	pub use crate::bam::BamFeatures;
}
