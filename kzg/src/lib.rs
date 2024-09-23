#![deny(warnings, unused, future_incompatible, nonstandard_style)]
use std::error::Error as ErrorTrait;

pub mod trivial_kzg;
pub mod batch_kzg;
pub mod biv_trivial_kzg;

pub type Error = Box<dyn ErrorTrait>;


