pub mod header;
pub mod pileup;
pub mod reader;

#[derive(Debug)]
pub enum TagValueType
{
	Char(char),
	I8(i8),
	U8(u8),
	I16(i16),
	U16(u16),
	I32(i32),
	U32(u32),
	F32(f32),
	String(String),
}

#[derive(Debug)]
pub struct Cigar
{
	pub length: u32,
	pub opcode: char,
}

pub struct Tag
{
	pub name: String,
	pub val_type: char,
	pub value: Vec<TagValueType>,
}

pub struct Field<'a>
{
	pub ref_id: i32,
	pub pos: i32,
	pub mapq: u8,
	pub bin: u16,
	pub flags: u16,
	pub next_ref_id: i32,
	pub next_pos: i32,
	pub tlen: i32,

	pub read_name: String,
	pub sequence: &'a [u8],
	pub sequence_quality: &'a [u8],
	pub cigar: &'a Vec<Cigar>,
	pub tags: &'a Vec<Tag>,
}
