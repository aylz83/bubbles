#[derive(Debug)]
pub struct TID
{
	pub name: String,
	pub length: u32,
}

#[derive(Debug)]
pub struct Header
{
	pub header: String,
	pub references: Vec<TID>,
}
