//! A port of libstdc++'s `std::sort` (introsort) and `std::partial_sort`.
//!
//! `std::sort` is not stable, and ULTRA sorts repeats with a comparator under
//! which distinct repeats can compare equal (the same repeat found in two
//! overlapping windows). Which of the tied elements ends up first decides
//! which one survives the overlap join, so reproducing the C++ output byte
//! for byte needs the same permutation libstdc++ produces, not just a sorted
//! result. Rust's own sorts order tied elements differently.
//!
//! The functions here follow `bits/stl_algo.h` and `bits/stl_heap.h` step by
//! step. Where libstdc++ lifts a value out and shifts elements over a hole,
//! this code keeps the value in the hole and swaps or rotates, which yields
//! the same sequence of comparisons and the same final arrangement without
//! needing `Clone` or `unsafe`.

/// Ranges at or below this length are left to insertion sort
/// (`_S_threshold` in libstdc++).
const THRESHOLD: usize = 16;

/// Sorts `v` exactly as `std::sort(v.begin(), v.end(), less)` does in
/// libstdc++, including the final order of elements that compare equal.
///
/// `less(a, b)` must be a strict weak ordering returning true when `a` goes
/// before `b`.
pub fn std_sort<T, F>(v: &mut [T], mut less: F)
where
    F: FnMut(&T, &T) -> bool,
{
    let n = v.len();
    if n == 0 {
        return;
    }
    introsort_loop(v, 0, n, lg(n) * 2, &mut less);
    final_insertion_sort(v, 0, n, &mut less);
}

/// Sorts `v` exactly as `std::partial_sort(v.begin(), v.end(), v.end(), less)`
/// does in libstdc++ (a heap sort). This is the fallback `std::sort` uses when
/// quicksort recursion gets too deep.
pub fn std_heap_sort<T, F>(v: &mut [T], mut less: F)
where
    F: FnMut(&T, &T) -> bool,
{
    let n = v.len();
    partial_sort_all(v, 0, n, &mut less);
}

/// `std::__lg`: floor of log2 for a positive value.
#[inline]
fn lg(n: usize) -> usize {
    (usize::BITS - 1 - n.leading_zeros()) as usize
}

/// `std::__introsort_loop` on `v[first..last]`.
fn introsort_loop<T, F>(v: &mut [T], first: usize, mut last: usize, mut depth_limit: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    while last - first > THRESHOLD {
        if depth_limit == 0 {
            partial_sort_all(v, first, last, less);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(v, first, last, less);
        introsort_loop(v, cut, last, depth_limit, less);
        last = cut;
    }
}

/// `std::__unguarded_partition_pivot`.
fn unguarded_partition_pivot<T, F>(v: &mut [T], first: usize, last: usize, less: &mut F) -> usize
where
    F: FnMut(&T, &T) -> bool,
{
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, less);
    unguarded_partition(v, first + 1, last, first, less)
}

/// `std::__move_median_to_first`: swaps the median of `a`, `b`, `c` into
/// `result`.
fn move_median_to_first<T, F>(v: &mut [T], result: usize, a: usize, b: usize, c: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            v.swap(result, b);
        } else if less(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if less(&v[a], &v[c]) {
        v.swap(result, a);
    } else if less(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

/// `std::__unguarded_partition`. The scans rely on the pivot and the
/// median-of-three arrangement to stop them, as in libstdc++; the bounds
/// checks here only matter for a comparator that is not a strict weak order.
fn unguarded_partition<T, F>(v: &mut [T], mut first: usize, mut last: usize, pivot: usize, less: &mut F) -> usize
where
    F: FnMut(&T, &T) -> bool,
{
    loop {
        while less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

/// `std::__final_insertion_sort`.
fn final_insertion_sort<T, F>(v: &mut [T], first: usize, last: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    if last - first > THRESHOLD {
        insertion_sort(v, first, first + THRESHOLD, less);
        // `std::__unguarded_insertion_sort`
        for i in first + THRESHOLD..last {
            unguarded_linear_insert(v, first, i, less);
        }
    } else {
        insertion_sort(v, first, last, less);
    }
}

/// `std::__insertion_sort`.
fn insertion_sort<T, F>(v: &mut [T], first: usize, last: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    if first == last {
        return;
    }
    for i in first + 1..last {
        if less(&v[i], &v[first]) {
            // Lift v[i], shift [first, i) up by one, drop it at the front.
            v[first..=i].rotate_right(1);
        } else {
            unguarded_linear_insert(v, first, i, less);
        }
    }
}

/// `std::__unguarded_linear_insert`: walks the element at `last` down past
/// every predecessor it is less than. `floor` is the lowest index the walk may
/// reach; libstdc++ has no such guard and relies on a sentinel being present.
fn unguarded_linear_insert<T, F>(v: &mut [T], floor: usize, last: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    let mut pos = last;
    while pos > floor && less(&v[last], &v[pos - 1]) {
        pos -= 1;
    }
    if pos != last {
        v[pos..=last].rotate_right(1);
    }
}

/// `std::__partial_sort(first, last, last)`: `__heap_select` (which with
/// `middle == last` is just `__make_heap`) followed by `__sort_heap`.
fn partial_sort_all<T, F>(v: &mut [T], first: usize, last: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    let heap = &mut v[first..last];
    make_heap(heap, less);

    // `std::__sort_heap`
    let mut end = heap.len();
    while end > 1 {
        end -= 1;
        // `std::__pop_heap(first, last, result = last)`: the old top goes to
        // `result`, and the value that was there is sifted in from the root.
        heap.swap(0, end);
        adjust_heap(heap, 0, end, less);
    }
}

/// `std::__make_heap`.
fn make_heap<T, F>(heap: &mut [T], less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    let len = heap.len();
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        adjust_heap(heap, parent, len, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// `std::__adjust_heap` with the value to place sitting at `hole_index`.
///
/// libstdc++ first walks the hole down to a leaf, always promoting the larger
/// child without looking at the value, and then bubbles the value back up
/// with `__push_heap`. The same two phases are done here with swaps.
fn adjust_heap<T, F>(heap: &mut [T], mut hole_index: usize, len: usize, less: &mut F)
where
    F: FnMut(&T, &T) -> bool,
{
    let top_index = hole_index;
    let mut second_child = hole_index;
    while second_child < (len - 1) / 2 {
        second_child = 2 * (second_child + 1);
        if less(&heap[second_child], &heap[second_child - 1]) {
            second_child -= 1;
        }
        heap.swap(hole_index, second_child);
        hole_index = second_child;
    }
    if len & 1 == 0 && second_child == (len - 2) / 2 {
        second_child = 2 * (second_child + 1);
        heap.swap(hole_index, second_child - 1);
        hole_index = second_child - 1;
    }

    // `std::__push_heap`
    while hole_index > top_index {
        let parent = (hole_index - 1) / 2;
        if !less(&heap[parent], &heap[hole_index]) {
            break;
        }
        heap.swap(hole_index, parent);
        hole_index = parent;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Same generator as the C++ harness that produced the expected hashes.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    /// Builds `n` (key, original index) pairs. `pattern` selects the key
    /// distribution; see the C++ harness for the matching definitions.
    fn make_input(n: usize, key_range: u64, pattern: u32, seed: u64) -> Vec<(u64, u32)> {
        let mut rng = Lcg(seed);
        (0..n)
            .map(|i| {
                let r = rng.next();
                let key = match pattern {
                    0 => r % key_range,
                    1 => i as u64 % key_range,
                    2 => (n - i) as u64 % key_range,
                    3 => {
                        // organ pipe
                        let half = n / 2;
                        (if i < half { i } else { n - i }) as u64 % key_range
                    }
                    _ => (i as u64 / 3 + r % 2) % key_range,
                };
                (key, i as u32)
            })
            .collect()
    }

    /// FNV-style hash of the output order of the original indices.
    fn order_hash(v: &[(u64, u32)]) -> u64 {
        let mut h: u64 = 14695981039346656037;
        for &(_, tag) in v {
            h = h.wrapping_mul(1099511628211) ^ tag as u64;
        }
        h
    }

    const SIZES: [usize; 22] = [
        0, 1, 2, 3, 5, 15, 16, 17, 18, 31, 32, 33, 34, 64, 100, 257, 1000, 1023, 1024, 5000, 20000, 65537,
    ];
    const RANGES: [u64; 7] = [1, 2, 3, 10, 100, 1000, 1 << 30];

    /// Combined hash over every (size, range, pattern, direction) case, in the
    /// order the C++ harness visits them. `kind` 0 is `std::sort`, 1 is
    /// `std::partial_sort(first, last, last)`.
    fn combined_hash(kind: u32) -> u64 {
        let mut total: u64 = 0;
        let mut seed = 1u64;
        for &n in &SIZES {
            for &range in &RANGES {
                for pattern in 0..5u32 {
                    for descending in [false, true] {
                        seed += 1;
                        let mut v = make_input(n, range, pattern, seed);
                        let less = |a: &(u64, u32), b: &(u64, u32)| {
                            if descending {
                                a.0 > b.0
                            } else {
                                a.0 < b.0
                            }
                        };
                        if kind == 0 {
                            std_sort(&mut v, less);
                        } else {
                            std_heap_sort(&mut v, less);
                        }
                        for w in v.windows(2) {
                            assert!(!less(&w[1], &w[0]), "not sorted: n={n} range={range} pattern={pattern}");
                        }
                        total = total.rotate_left(7) ^ order_hash(&v);
                    }
                }
            }
        }
        total
    }

    #[test]
    fn sort_matches_libstdcxx_permutation() {
        assert_eq!(combined_hash(0), EXPECTED_SORT_HASH);
    }

    #[test]
    fn heap_sort_matches_libstdcxx_permutation() {
        assert_eq!(combined_hash(1), EXPECTED_HEAP_HASH);
    }

    #[test]
    fn small_cases_spelled_out() {
        // Hand-checkable: all keys equal, so the result is whatever
        // permutation the algorithm leaves, which for insertion sort is the
        // input order.
        let mut v: Vec<(u64, u32)> = (0..10).map(|i| (7, i)).collect();
        std_sort(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v.iter().map(|p| p.1).collect::<Vec<_>>(), (0..10).collect::<Vec<_>>());

        let mut v = vec![(3u64, 0u32), (1, 1), (2, 2), (1, 3), (3, 4)];
        std_sort(&mut v, |a, b| a.0 < b.0);
        assert_eq!(v, vec![(1, 1), (1, 3), (2, 2), (3, 0), (3, 4)]);
    }

    // Produced by the g++ 13 harness (libstdc++), see the module docs.
    const EXPECTED_SORT_HASH: u64 = 10351184119196295753;
    const EXPECTED_HEAP_HASH: u64 = 16425017049699444346;
}
