use pyo3::prelude::*;


fn cardinality(seq: &str, k: usize) -> usize {
    let substrings: std::collections::HashSet<&str> = (0..=seq.len() - k).map(|i| &seq[i..i + k]).collect();
    substrings.len()
}

#[pyfunction]
fn delta(_p: Python, seq: &str) -> PyResult<f64> {
    if seq.is_empty() {
        Ok(0.0)
    } else {
        let mut max_ratio = 0.0;
        for k in 1..=seq.len() {
            let ratio = cardinality(seq, k) as f64 / k as f64;
             if max_ratio < ratio {
                max_ratio = ratio;
        // max_ratio = max_ratio.max(ratio);
    }

        }
        Ok(max_ratio)
    }
}


/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "delta_rust")]
fn delta_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(delta, m)?)?;
    Ok(())
}

