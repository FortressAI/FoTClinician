#!/usr/bin/env python3
"""
FoT Audio Readiness Checker

Validates medical audio data against Field of Truth clinical readiness requirements.
Checks for calibration, quality metrics, and technical specifications.
"""

import json
import sys
from typing import Dict, List, Any, Optional

def readiness_report(sig: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate comprehensive readiness report for medical audio signal.
    
    Args:
        sig: Medical audio signal metadata
        
    Returns:
        Dictionary with readiness status, missing fields, and warnings
    """
    missing, warn = [], []

    def need(k: str) -> None:
        """Check if required field is present and non-empty"""
        if sig.get(k) in (None, "", []):
            missing.append(k)

    # Required technical fields
    for k in ("bodySite", "sampleRateHz", "bitDepth", "channels", "durationSec", "deviceModel", "acquiredAt"):
        need(k)

    # Calibration check
    if sig.get("calibrationPassed") not in (1, True):
        warn.append("calibration not passed")

    # Quality measurements analysis
    q = {m.get("hasMetric"): m for m in sig.get("qualityMeasurements", [])} if sig.get("qualityMeasurements") else {}
    
    def qv(metric: str) -> Optional[float]:
        """Extract quality metric value"""
        m = q.get(metric)
        return None if m is None else m.get("value")

    snr = qv("faud:Quality_SNR_dB")
    nf = qv("faud:Quality_NoiseFloor_dBFS")
    art = qv("faud:Quality_ArtifactScore")

    quality_ok = True
    
    # Technical specifications check
    try:
        if int(sig.get("sampleRateHz", 0)) < 4000:
            quality_ok = False
            warn.append("sampleRateHz < 4000")
        
        if int(sig.get("bitDepth", 0)) < 16:
            quality_ok = False
            warn.append("bitDepth < 16")
        
        if int(sig.get("channels", 0)) != 1:
            warn.append("channels != 1 (not blocking)")
        
        if float(sig.get("durationSec", 0)) < 10.0:
            warn.append("duration < 10s (may harm performance)")
            
    except Exception:
        missing.append("numeric audio attributes malformed")

    # Quality metrics check
    if snr is not None and snr < 20:
        quality_ok = False
        warn.append("SNR < 20 dB")
    
    if nf is not None and nf > -35:
        quality_ok = False
        warn.append("noise floor > -35 dBFS")
    
    if art is not None and art > 0.4:
        quality_ok = False
        warn.append("artifact score > 0.4")

    # Overall readiness determination
    ready = (len(missing) == 0 and quality_ok)
    
    return {
        "ready": ready,
        "missing": missing,
        "warnings": warn,
        "quality_metrics": {
            "snr_db": snr,
            "noise_floor_dbfs": nf,
            "artifact_score": art,
            "calibration_passed": sig.get("calibrationPassed") in (1, True)
        },
        "technical_specs": {
            "sample_rate_hz": sig.get("sampleRateHz"),
            "bit_depth": sig.get("bitDepth"),
            "channels": sig.get("channels"),
            "duration_sec": sig.get("durationSec")
        }
    }

def main():
    """Command-line interface for audio readiness checking"""
    if len(sys.argv) != 2:
        print("Usage: python audio_readiness_checker.py <audio_metadata.json>")
        sys.exit(1)
    
    try:
        with open(sys.argv[1], 'r') as f:
            sig = json.load(f)
        
        report = readiness_report(sig)
        print(json.dumps(report, indent=2))
        
    except FileNotFoundError:
        print(f"Error: File {sys.argv[1]} not found")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {sys.argv[1]}: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
