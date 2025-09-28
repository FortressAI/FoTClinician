"""
Setup script for compiling the quantum_ops C extension
Optimized for performance-critical quantum operations

Usage:
    python core/setup_quantum_ops.py build_ext --inplace
"""

from setuptools import setup, Extension
import numpy as np
import platform
import os

# Platform-specific optimizations
extra_compile_args = ["-O3", "-ffast-math"]
extra_link_args = []

if platform.system() == "Darwin":  # macOS
    # Apple Silicon optimizations (avoid -march=native issues)
    if platform.machine() == "arm64":
        extra_compile_args.extend(["-mcpu=apple-a14"])  # Apple Silicon
    else:
        extra_compile_args.extend(["-march=x86-64"])    # Intel Mac
    # Link against Accelerate framework for optimized BLAS
    extra_link_args.extend(["-framework", "Accelerate"])
    
elif platform.system() == "Linux":
    # Linux optimizations
    extra_compile_args.extend([
        "-march=native",
        "-mtune=native",
        "-fopenmp"  # OpenMP for parallel operations
    ])
    extra_link_args.extend(["-lgomp", "-lopenblas"])

# Define the extension module
quantum_ops_ext = Extension(
    name='quantum_ops',
    sources=['core/quantum_ops.c'],
    include_dirs=[
        np.get_include(),
        '/usr/local/include',  # For BLAS headers
        '/opt/homebrew/include'  # For Apple Silicon homebrew
    ],
    library_dirs=[
        '/usr/local/lib',
        '/opt/homebrew/lib'
    ],
    libraries=['blas', 'lapack'],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    language='c'
)

setup(
    name='quantum_ops',
    version='1.0.0',
    description='High-performance quantum operations for FoTChemistry',
    ext_modules=[quantum_ops_ext],
    zip_safe=False,
    python_requires='>=3.8'
)

if __name__ == "__main__":
    print("🔨 Compiling quantum operations C extension...")
    print("💡 This will provide significant speedup for 8096-dimensional quantum operations")
    print()
    print("Platform optimizations:")
    print(f"  - System: {platform.system()}")
    print(f"  - Architecture: {platform.machine()}")
    print(f"  - Compile flags: {' '.join(extra_compile_args)}")
    print(f"  - Link flags: {' '.join(extra_link_args)}")
