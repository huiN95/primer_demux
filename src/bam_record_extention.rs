use rust_htslib::bam::{record::Aux, record::Record};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ReadRecord {
    pub id: Arc<str>,
    pub sequence: Vec<u8>,
    pub quality: Option<Vec<u8>>,
    pub dw: Option<Vec<u8>>,
    pub cr: Option<Vec<u8>>,
    pub rq: Option<f32>, // Set default to 0.8 to obtain more data for SMC
    pub np: Option<i32>,
    pub cx: Option<i32>,
    pub ch: Option<i32>,
    pub sn: Option<Vec<f32>>,
    pub RG: Option<String>,
    pub be: Option<Vec<u32>>,
    pub cq: Option<f32>,
    pub ar: Option<Vec<u32>>,
    pub nn: Option<Vec<u8>>,
}

impl Default for ReadRecord {
    fn default() -> Self {
        ReadRecord {
            id: Arc::from(""),
            sequence: Vec::new(),
            quality: None,
            dw: None,
            cr: None,
            rq: Some(0.8),
            np: None,
            cx: None,
            ch: None,
            sn: Some(vec![20.0, 20.0, 20.0, 20.0]),
            RG: Some("0425".to_string()),
            be: None,
            cq: None,
            ar: None,
            nn: None,
        }
    }
}

pub trait AuxExt<'a> {
    fn i32(&self, tag: &[u8]) -> Option<i32>; // ch:i, cx:i …
    fn u32(&self, tag: &[u8]) -> Option<u32>; // np:i, cx:i …
    fn float(&self, tag: &[u8]) -> Option<f32>; // rq:f, cq:f …
    fn array_u8(&self, tag: &[u8]) -> Option<Vec<u8>>; // dw:B:C …
    fn array_f32(&self, tag: &[u8]) -> Option<Vec<f32>>; // sn:B:f …
    fn array_u32(&self, tag: &[u8]) -> Option<Vec<u32>>; // sn:B:f …
    fn string(&self, tag: &[u8]) -> Option<&str>; // RG:Z …
}

/// —— ② Implementation for external types ——
/// As long as the trait is defined in your crate, you can implement it for external types.
impl<'a> AuxExt<'a> for Record {
    fn i32(&self, tag: &[u8]) -> Option<i32> {
        match self.aux(tag).ok()? {
            Aux::I32(i) => Some(i),
            Aux::I16(i) => Some(i as i32),
            Aux::I8(i) => Some(i as i32),
            _ => None,
        }
    }
    /* ---------- Pure integers ---------- */
    fn u32(&self, tag: &[u8]) -> Option<u32> {
        match self.aux(tag).ok()? {
            Aux::U32(i) => Some(i),
            Aux::U16(i) => Some(i as u32),
            Aux::U8(i) => Some(i as u32),
            _ => None,
        }
    }

    /* ---------- Single precision floats ---------- */
    fn float(&self, tag: &[u8]) -> Option<f32> {
        match self.aux(tag).ok()? {
            Aux::Float(f) => Some(f),
            _ => None,
        }
    }

    /* ---------- u8 arrays ---------- */
    fn array_u8(&self, tag: &[u8]) -> Option<Vec<u8>> {
        match self.aux(tag).ok()? {
            Aux::ArrayU8(a) => Some(a.iter().collect()),
            Aux::ArrayU16(a) => Some(
                a.iter()
                    .map(|x| if x < 255 { x as u8 } else { 255 })
                    .collect(),
            ),
            Aux::ArrayU32(a) => Some(
                a.iter()
                    .map(|x| if x < 255 { x as u8 } else { 255 })
                    .collect(),
            ),
            Aux::ArrayFloat(a) => Some(
                a.iter()
                    .map(|x| if x < 255.0 { x as u8 } else { 255 })
                    .collect(),
            ),

            _ => None,
        }
    }

    /* ---------- f32 arrays ---------- */
    fn array_f32(&self, tag: &[u8]) -> Option<Vec<f32>> {
        match self.aux(tag).ok()? {
            Aux::ArrayFloat(a) => Some(a.iter().collect()),
            _ => None,
        }
    }
    fn array_u32(&self, tag: &[u8]) -> Option<Vec<u32>> {
        match self.aux(tag).ok()? {
            Aux::ArrayU32(a) => Some(a.iter().collect()),
            Aux::ArrayU8(a) => Some(a.iter().map(|x| x as u32).collect()),
            Aux::ArrayU16(a) => Some(a.iter().map(|x| x as u32).collect()),
            _ => None,
        }
    }
    fn string(&self, tag: &[u8]) -> Option<&str> {
        match self.aux(tag).ok()? {
            Aux::String(s) => Some(s),
            _ => None,
        }
    }
}
