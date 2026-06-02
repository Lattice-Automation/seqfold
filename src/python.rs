//! PyO3 bindings: the `seqfold._core` extension module.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use crate::core::fold as cfold;
use crate::core::pyfloat::pyfloat_str;
use crate::core::tm as ctm;

const FMT: &str = "{:>4} {:>4} {:>6}  {:<15}";

/// A single structure with a free energy, description, and inward children.
#[pyclass(name = "Struct")]
#[derive(Clone)]
struct Struct {
    #[pyo3(get)]
    e: f64,
    #[pyo3(get)]
    desc: String,
    ij: Vec<(i64, i64)>,
}

impl Struct {
    fn from_core(s: cfold::Struct) -> Struct {
        Struct {
            e: s.e,
            desc: s.desc,
            // core stores 32-bit indices; Python exposes ints (i64 tuples)
            ij: s.ij.iter().map(|&(a, b)| (a as i64, b as i64)).collect(),
        }
    }
}

#[pymethods]
impl Struct {
    #[new]
    #[pyo3(signature = (e = f64::NEG_INFINITY, desc = String::new(), ij = Vec::new()))]
    fn new(e: f64, desc: String, ij: Vec<(i64, i64)>) -> Struct {
        Struct { e, desc, ij }
    }

    #[classattr]
    fn fmt() -> &'static str {
        FMT
    }

    #[getter]
    fn ij(&self) -> Vec<(i64, i64)> {
        self.ij.clone()
    }

    fn with_ij(&self, ij: Vec<(i64, i64)>) -> Struct {
        Struct {
            e: self.e,
            desc: self.desc.clone(),
            ij,
        }
    }

    fn __eq__(&self, other: &Struct) -> bool {
        self.e == other.e && self.ij == other.ij
    }

    fn __bool__(&self) -> bool {
        self.e != f64::INFINITY && self.e != f64::NEG_INFINITY
    }

    fn __str__(&self) -> String {
        let (i, j) = if self.ij.is_empty() {
            (String::new(), String::new())
        } else {
            (self.ij[0].0.to_string(), self.ij[0].1.to_string())
        };
        let e = pyfloat_str(self.e);
        // Mirror Python's str.format width/alignment for the four columns.
        format!(
            "{:>4} {:>4} {:>6}  {:<15}",
            i, j, e, self.desc
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

fn runtime(e: cfold::FoldError) -> PyErr {
    PyRuntimeError::new_err(e.0)
}

fn value(e: ctm::TmError) -> PyErr {
    PyValueError::new_err(e.0)
}

#[pyfunction]
#[pyo3(signature = (seq, temp = 37.0))]
fn fold(seq: &str, temp: f64) -> PyResult<Vec<Struct>> {
    cfold::fold(seq, temp)
        .map(|v| v.into_iter().map(Struct::from_core).collect())
        .map_err(runtime)
}

#[pyfunction]
#[pyo3(signature = (seq, temp = 37.0))]
fn dg(seq: &str, temp: f64) -> PyResult<f64> {
    cfold::dg(seq, temp).map_err(runtime)
}

#[pyfunction]
#[pyo3(signature = (seq, temp = 37.0))]
fn dg_cache(seq: &str, temp: f64) -> PyResult<Vec<Vec<f64>>> {
    cfold::dg_cache(seq, temp).map_err(runtime)
}

#[pyfunction]
fn dot_bracket(seq: &str, structs: Vec<Struct>) -> String {
    let mut result = vec![b'.'; seq.len()];
    for s in &structs {
        if s.ij.len() == 1 {
            let (i, j) = s.ij[0];
            result[i as usize] = b'(';
            result[j as usize] = b')';
        }
    }
    String::from_utf8(result).expect("ascii")
}

#[pyfunction]
#[pyo3(signature = (seq1, seq2 = String::new(), pcr = true))]
fn tm(seq1: &str, seq2: String, pcr: bool) -> PyResult<f64> {
    ctm::tm(seq1, &seq2, pcr).map_err(value)
}

#[pyfunction]
#[pyo3(signature = (seq1, seq2 = String::new(), pcr = true))]
fn tm_cache(seq1: &str, seq2: String, pcr: bool) -> PyResult<Vec<Vec<f64>>> {
    ctm::tm_cache(seq1, &seq2, pcr).map_err(value)
}

#[pyfunction]
fn gc_cache(seq: &str) -> Vec<Vec<f64>> {
    ctm::gc_cache(seq)
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Struct>()?;
    m.add_function(wrap_pyfunction!(fold, m)?)?;
    m.add_function(wrap_pyfunction!(dg, m)?)?;
    m.add_function(wrap_pyfunction!(dg_cache, m)?)?;
    m.add_function(wrap_pyfunction!(dot_bracket, m)?)?;
    m.add_function(wrap_pyfunction!(tm, m)?)?;
    m.add_function(wrap_pyfunction!(tm_cache, m)?)?;
    m.add_function(wrap_pyfunction!(gc_cache, m)?)?;
    Ok(())
}
