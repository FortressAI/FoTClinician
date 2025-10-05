#!/usr/bin/env python3
"""
FoT Synthetic Medical Data Generator

Generates synthetic medical images and audio for testing and demonstration.
Creates placeholder files with realistic metadata for clinical pipeline validation.
"""

import json
import os
import wave
import struct
import math
import datetime
from typing import Dict, Any

def create_directories():
    """Create necessary directories for synthetic data"""
    os.makedirs("data/synth/images", exist_ok=True)
    os.makedirs("data/synth/audio", exist_ok=True)
    os.makedirs("data/synth/metadata", exist_ok=True)

def write_synthetic_wav(path: str, seconds: float = 12.0, fs: int = 8000, 
                       freq: float = 180.0, amp: float = 0.25, noise: float = 0.02) -> None:
    """
    Generate synthetic WAV file with specified parameters
    
    Args:
        path: Output file path
        seconds: Duration in seconds
        fs: Sample rate in Hz
        freq: Primary frequency component
        amp: Signal amplitude
        noise: Noise level
    """
    n = int(seconds * fs)
    
    with wave.open(path, "w") as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(fs)
        
        for i in range(n):
            t = i / fs
            # Generate signal with primary frequency and noise
            s = amp * math.sin(2 * math.pi * freq * t) + noise * (2 * (i % 17) / 17.0 - 1.0)
            s = max(-0.99, min(0.99, s))
            wf.writeframesraw(struct.pack("<h", int(s * 32767)))

def generate_lung_audio() -> Dict[str, Any]:
    """Generate synthetic lung audio with good quality metrics"""
    audio_path = "data/synth/audio/lung_good.wav"
    write_synthetic_wav(audio_path, seconds=20, fs=8000, freq=180, amp=0.25, noise=0.02)
    
    metadata = {
        "wave_file": audio_path,
        "metadata_file": "tests/fot_audio/fixtures/lung_good.json",
        "description": "Synthetic lung audio with good quality metrics",
        "clinical_context": "Normal breath sounds, posterior lower left",
        "quality_assessment": "High SNR, low noise floor, minimal artifacts"
    }
    
    return metadata

def generate_heart_audio() -> Dict[str, Any]:
    """Generate synthetic heart audio with poor quality metrics"""
    audio_path = "data/synth/audio/heart_nearmiss.wav"
    write_synthetic_wav(audio_path, seconds=8, fs=4000, freq=100, amp=0.20, noise=0.10)
    
    metadata = {
        "wave_file": audio_path,
        "metadata_file": "tests/fot_audio/fixtures/heart_nearmiss.json",
        "description": "Synthetic heart audio with poor quality metrics",
        "clinical_context": "Possible murmur, cardiac apex",
        "quality_assessment": "Low SNR, high noise floor, significant artifacts"
    }
    
    return metadata

def generate_retina_image_metadata() -> Dict[str, Any]:
    """Generate synthetic retina image metadata"""
    return {
        "image_file": "data/synth/images/retina_good.png",
        "metadata_file": "tests/fot_image/fixtures/retina_good.json",
        "description": "Synthetic retina fundus image with good quality",
        "clinical_context": "Diabetic retinopathy screening, left macula",
        "quality_assessment": "High focus score, proper exposure, good SNR"
    }

def generate_derm_image_metadata() -> Dict[str, Any]:
    """Generate synthetic dermatology image metadata"""
    return {
        "image_file": "data/synth/images/derm_nearmiss.png",
        "metadata_file": "tests/fot_image/fixtures/derm_nearmiss.json",
        "description": "Synthetic dermatology image with poor quality",
        "clinical_context": "Mole screening, left forearm",
        "quality_assessment": "Poor focus, overexposed, low SNR"
    }

def create_synthetic_image_placeholders():
    """Create placeholder PNG files for synthetic images"""
    # Note: In a real implementation, you would generate actual PNG files
    # For now, we'll create metadata files that reference expected image locations
    
    retina_placeholder = {
        "type": "synthetic_image",
        "modality": "fundus",
        "dimensions": [1024, 1024],
        "description": "Synthetic retina fundus image placeholder",
        "note": "Replace with actual PNG generation in production"
    }
    
    derm_placeholder = {
        "type": "synthetic_image", 
        "modality": "derm",
        "dimensions": [480, 360],
        "description": "Synthetic dermatology image placeholder",
        "note": "Replace with actual PNG generation in production"
    }
    
    with open("data/synth/images/retina_good.json", "w") as f:
        json.dump(retina_placeholder, f, indent=2)
    
    with open("data/synth/images/derm_nearmiss.json", "w") as f:
        json.dump(derm_placeholder, f, indent=2)

def main():
    """Main function to generate all synthetic test data"""
    print("🔬 FoT Synthetic Medical Data Generator")
    print("======================================")
    
    # Create directories
    create_directories()
    
    # Generate audio files
    print("🎵 Generating synthetic audio files...")
    lung_meta = generate_lung_audio()
    heart_meta = generate_heart_audio()
    
    # Generate image metadata
    print("🖼️  Generating synthetic image metadata...")
    retina_meta = generate_retina_image_metadata()
    derm_meta = generate_derm_image_metadata()
    
    # Create image placeholders
    create_synthetic_image_placeholders()
    
    # Create master metadata file
    master_metadata = {
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "generator_version": "1.0",
        "description": "Synthetic medical data for FoT clinical testing",
        "audio_files": [lung_meta, heart_meta],
        "image_files": [retina_meta, derm_meta],
        "usage_notes": [
            "These are synthetic files for testing purposes only",
            "Do not use for actual clinical diagnosis",
            "Replace with real medical data for production use",
            "Respect all dataset licenses and terms of use"
        ]
    }
    
    with open("data/synth/metadata/synthetic_data_manifest.json", "w") as f:
        json.dump(master_metadata, f, indent=2)
    
    print("")
    print("✅ Synthetic data generation complete!")
    print("")
    print("📁 Generated files:")
    print("   🎵 data/synth/audio/lung_good.wav")
    print("   🎵 data/synth/audio/heart_nearmiss.wav")
    print("   🖼️  data/synth/images/retina_good.json")
    print("   🖼️  data/synth/images/derm_nearmiss.json")
    print("   📋 data/synth/metadata/synthetic_data_manifest.json")
    print("")
    print("🧪 Test with readiness checkers:")
    print("   python tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/lung_good.json")
    print("   python tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/heart_nearmiss.json")
    print("   python tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/retina_good.json")
    print("   python tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/derm_nearmiss.json")

if __name__ == "__main__":
    main()
