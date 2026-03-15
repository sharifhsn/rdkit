use bitvec::prelude::*;
use cxx::SharedPtr;

#[derive(Clone, Debug)]
pub struct Fingerprint(pub BitVec<u8, bitvec::order::Lsb0>);

impl Fingerprint {
    pub fn new(ptr: SharedPtr<rdkit_sys::fingerprint_ffi::ExplicitBitVect>) -> Self {
        let num_bits = rdkit_sys::fingerprint_ffi::explicit_bit_vect_num_bits(&ptr) as usize;
        let unique_ptr_bytes = rdkit_sys::fingerprint_ffi::explicit_bit_vect_to_u64_vec(&ptr);
        let u8_bytes: Vec<u8> = unique_ptr_bytes
            .iter()
            .flat_map(|&w| w.to_le_bytes())
            .collect();
        let mut bv = BitVec::<u8, Lsb0>::from_vec(u8_bytes);
        bv.truncate(num_bits);
        Fingerprint(bv)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn tanimoto_distance(&self, other: &Fingerprint) -> f32 {
        let a = self.0.as_raw_slice();
        let b = other.0.as_raw_slice();
        let mut and_count: u32 = 0;
        let mut a_count: u32 = 0;
        let mut b_count: u32 = 0;
        for (&x, &y) in a.iter().zip(b.iter()) {
            and_count += (x & y).count_ones();
            a_count += x.count_ones();
            b_count += y.count_ones();
        }
        let union_count = a_count + b_count - and_count;
        and_count as f32 / union_count as f32
    }
}
