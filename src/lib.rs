#![feature(portable_simd)]
#![feature(path_add_extension)]

mod bai;
pub mod bam;

pub mod error;

pub trait AsyncReadSeek: tokio::io::AsyncRead + tokio::io::AsyncSeek {}
impl<T: tokio::io::AsyncRead + tokio::io::AsyncSeek + ?Sized> AsyncReadSeek for T {}

// pub mod bgz;
