/*
 * Quantum Operations C Extension for FoTChemistry
 * High-performance matrix operations for 8096-dimensional vQbit states
 * 
 * Critical paths that need C acceleration:
 * - Matrix-vector multiplication for large quantum states
 * - Eigenvalue decomposition for virtue operators  
 * - Quantum Fourier Transform operations
 * - State normalization and fidelity calculations
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>
#include <math.h>

// Use Apple's Accelerate framework on macOS
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#define HAS_BLAS 1
#else
// Try to include BLAS if available, otherwise use naive implementation
#ifdef HAVE_CBLAS
#include <cblas.h>
#define HAS_BLAS 1
#else
#define HAS_BLAS 0
#endif
#endif

// Fast matrix-vector multiplication for quantum state evolution
static PyObject* fast_matvec(PyObject* self, PyObject* args) {
    PyArrayObject *matrix, *vector, *result;
    
    if (!PyArg_ParseTuple(args, "OO", &matrix, &vector)) {
        return NULL;
    }
    
    // Get array dimensions
    npy_intp n = PyArray_DIM(matrix, 0);
    npy_intp m = PyArray_DIM(matrix, 1);
    
    if (PyArray_DIM(vector, 0) != m) {
        PyErr_SetString(PyExc_ValueError, "Matrix and vector dimensions don't match");
        return NULL;
    }
    
    // Create result array
    npy_intp dims[1] = {n};
    result = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_COMPLEX128, 0);
    
    // Get data pointers
    double complex* mat_data = (double complex*)PyArray_DATA(matrix);
    double complex* vec_data = (double complex*)PyArray_DATA(vector);
    double complex* res_data = (double complex*)PyArray_DATA(result);
    
    // Perform optimized matrix-vector multiplication
#if HAS_BLAS
    // Using BLAS for maximum performance
    double complex alpha = 1.0 + 0.0*I;
    double complex beta = 0.0 + 0.0*I;
    cblas_zgemv(CblasRowMajor, CblasNoTrans, n, m, 
                &alpha, mat_data, m, vec_data, 1, &beta, res_data, 1);
#else
    // Naive implementation without BLAS
    for (npy_intp i = 0; i < n; i++) {
        res_data[i] = 0.0;
        for (npy_intp j = 0; j < m; j++) {
            res_data[i] += mat_data[i * m + j] * vec_data[j];
        }
    }
#endif
    
    return (PyObject*)result;
}

// Fast quantum state normalization
static PyObject* normalize_state(PyObject* self, PyObject* args) {
    PyArrayObject *state;
    
    if (!PyArg_ParseTuple(args, "O", &state)) {
        return NULL;
    }
    
    npy_intp n = PyArray_DIM(state, 0);
    double complex* state_data = (double complex*)PyArray_DATA(state);
    
    // Calculate norm: ||ψ||² = Σ|αᵢ|²
    double norm_squared = 0.0;
    for (npy_intp i = 0; i < n; i++) {
        double complex amp = state_data[i];
        norm_squared += creal(amp * conj(amp));
    }
    
    double norm = sqrt(norm_squared);
    
    // Normalize in-place: ψ → ψ/||ψ||
    if (norm > 1e-15) {  // Avoid division by zero
        for (npy_intp i = 0; i < n; i++) {
            state_data[i] /= norm;
        }
    }
    
    return PyFloat_FromDouble(norm);
}

// Fast quantum fidelity calculation F = |⟨ψ₁|ψ₂⟩|²
static PyObject* quantum_fidelity(PyObject* self, PyObject* args) {
    PyArrayObject *state1, *state2;
    
    if (!PyArg_ParseTuple(args, "OO", &state1, &state2)) {
        return NULL;
    }
    
    npy_intp n = PyArray_DIM(state1, 0);
    if (PyArray_DIM(state2, 0) != n) {
        PyErr_SetString(PyExc_ValueError, "State dimensions don't match");
        return NULL;
    }
    
    double complex* data1 = (double complex*)PyArray_DATA(state1);
    double complex* data2 = (double complex*)PyArray_DATA(state2);
    
    // Calculate inner product ⟨ψ₁|ψ₂⟩ = Σ ψ₁*[i] * ψ₂[i]
    double complex inner_product = 0.0;
    for (npy_intp i = 0; i < n; i++) {
        inner_product += conj(data1[i]) * data2[i];
    }
    
    // Fidelity is |⟨ψ₁|ψ₂⟩|²
    double fidelity = creal(inner_product * conj(inner_product));
    
    return PyFloat_FromDouble(fidelity);
}

// Fast virtue operator expectation value ⟨ψ|V̂|ψ⟩
static PyObject* virtue_expectation(PyObject* self, PyObject* args) {
    PyArrayObject *state, *virtue_operator;
    
    if (!PyArg_ParseTuple(args, "OO", &state, &virtue_operator)) {
        return NULL;
    }
    
    npy_intp n = PyArray_DIM(state, 0);
    
    // For diagonal virtue operators (memory efficient)
    if (PyArray_NDIM(virtue_operator) == 1) {
        double complex* state_data = (double complex*)PyArray_DATA(state);
        double* virtue_data = (double*)PyArray_DATA(virtue_operator);
        
        double expectation = 0.0;
        for (npy_intp i = 0; i < n; i++) {
            double complex amp = state_data[i];
            double prob = creal(amp * conj(amp));  // |αᵢ|²
            expectation += prob * virtue_data[i];
        }
        
        return PyFloat_FromDouble(expectation);
    }
    
    // For full matrix virtue operators
    // V̂_exp = ⟨ψ|V̂|ψ⟩ = ψ† V̂ ψ
    // This would require more complex implementation
    Py_RETURN_NONE;
}

// Method definitions
static PyMethodDef QuantumOpsMethods[] = {
    {"fast_matvec", fast_matvec, METH_VARARGS, "Fast matrix-vector multiplication"},
    {"normalize_state", normalize_state, METH_VARARGS, "Normalize quantum state"},
    {"quantum_fidelity", quantum_fidelity, METH_VARARGS, "Calculate quantum fidelity"},
    {"virtue_expectation", virtue_expectation, METH_VARARGS, "Calculate virtue expectation value"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef quantum_ops_module = {
    PyModuleDef_HEAD_INIT,
    "quantum_ops",
    "High-performance quantum operations for FoTChemistry",
    -1,
    QuantumOpsMethods
};

// Module initialization
PyMODINIT_FUNC PyInit_quantum_ops(void) {
    import_array();  // Initialize NumPy C API
    return PyModule_Create(&quantum_ops_module);
}
