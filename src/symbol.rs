//! Nucleotide symbol encoding. Port of `Symbol.hpp` / `Symbol.cpp`.

/// A sequence symbol. Only `N_UNK`, `N_T`, `N_C`, `N_A` and `N_G` are ever
/// produced by [`symbol_for_char`]; the remaining codes exist for parity with
/// the C++ header.
pub type Symbol = u8;

pub const NUM_SYMBOLS: usize = 27;

pub const N_UNK: Symbol = 0;
pub const N_T: Symbol = 1;
pub const N_C: Symbol = 2;
pub const N_A: Symbol = 3;
pub const N_G: Symbol = 4;
pub const N_R: Symbol = 5;
pub const N_Y: Symbol = 6;
pub const N_K: Symbol = 7;
pub const N_M: Symbol = 8;
pub const N_S: Symbol = 9;
pub const N_W: Symbol = 10;
pub const N_B: Symbol = 11;
pub const N_D: Symbol = 12;
pub const N_H: Symbol = 13;
pub const N_V: Symbol = 14;
pub const N_N: Symbol = 15;
pub const N_GAP: Symbol = 16;

/// Maps an input byte to its symbol. Anything other than `ACGTacgt` is
/// `N_UNK`.
#[inline]
pub fn symbol_for_char(c: u8) -> Symbol {
    match c {
        b'a' | b'A' => N_A,
        b't' | b'T' => N_T,
        b'c' | b'C' => N_C,
        b'g' | b'G' => N_G,
        _ => N_UNK,
    }
}

/// Maps a symbol back to an upper-case byte. Anything other than the four
/// nucleotides is `b'N'`.
#[inline]
pub fn char_for_symbol(s: Symbol) -> u8 {
    match s {
        N_A => b'A',
        N_T => b'T',
        N_G => b'G',
        N_C => b'C',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_nucleotides() {
        for (lower, upper) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T')] {
            assert_eq!(symbol_for_char(lower), symbol_for_char(upper));
            assert_eq!(char_for_symbol(symbol_for_char(lower)), upper);
        }
    }

    #[test]
    fn unknown_bytes_become_n() {
        for c in [b'N', b'n', b'R', b'-', b'0', b' ', 0u8, 200u8] {
            assert_eq!(symbol_for_char(c), N_UNK);
        }
        assert_eq!(char_for_symbol(N_UNK), b'N');
        assert_eq!(char_for_symbol(N_N), b'N');
        assert_eq!(char_for_symbol(26), b'N');
    }
}
