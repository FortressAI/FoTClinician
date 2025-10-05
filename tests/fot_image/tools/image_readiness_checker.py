#!/usr/bin/env python3
"""
FoT Image Readiness Checker

Validates medical imaging data against Field of Truth clinical readiness requirements.
Checks for PHI compliance, quality metrics, and technical specifications.
"""

import json
import sys
from typing import Dict, List, Any, Optional

def readiness_report(study: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate comprehensive readiness report for medical imaging study.
    
    Args:
        study: Medical imaging study metadata
        
    Returns:
        Dictionary with readiness status, missing fields, and warnings
    """
    missing, warn = [], []

    def need(k: str) -> None:
        """Check if required field is present and non-empty"""
        if study.get(k) in (None, "", []):
            missing.append(k)

    # Required technical fields
    for k in ("modality", "bodySite", "acquiredAt", "deviceModel", "widthPx", "heightPx"):
        need(k)

    # Need either calibrated spacing or a scale reference
    if study.get("pixelSpacingMm") is None and not study.get("scaleRef"):
        missing.append("pixelSpacingMm OR scaleRef")

    # PHI compliance check
    if study.get("phiBurninFlag") not in (0, None):
        warn.append("phiBurninFlag indicates possible PHI; block until redacted")

    # Quality measurements analysis
    q = {m.get("hasMetric"): m for m in study.get("qualityMeasurements", [])} if study.get("qualityMeasurements") else {}
    
    def qv(metric: str) -> Optional[float]:
        """Extract quality metric value"""
        m = q.get(metric)
        return None if m is None else m.get("value")

    focus = qv("fimg:Quality_FocusScore")
    expos = qv("fimg:Quality_ExposureScore")
    snr = qv("fimg:Quality_SNR_dB")

    quality_ok = True
    
    # Focus quality check
    if focus is not None and focus < 0.6:
        quality_ok = False
        warn.append("focus < 0.6")
    
    # Exposure quality check
    if expos is not None and (expos < 0.3 or expos > 0.9):
        quality_ok = False
        warn.append("exposure out of 0.3–0.9")
    
    # Signal-to-noise ratio check
    if snr is not None and snr < 20:
        quality_ok = False
        warn.append("SNR < 20 dB")

    # Resolution check
    try:
        if min(int(study.get("widthPx", 0)), int(study.get("heightPx", 0))) < 512:
            warn.append("shortest side < 512px")
            quality_ok = False
    except Exception:
        missing.append("widthPx/heightPx integers")

    # Overall readiness determination
    ready = (len(missing) == 0 and quality_ok and study.get("phiBurninFlag") in (0, None))
    
    return {
        "ready": ready,
        "missing": missing,
        "warnings": warn,
        "quality_metrics": {
            "focus_score": focus,
            "exposure_score": expos,
            "snr_db": snr,
            "resolution_ok": quality_ok
        },
        "phi_compliant": study.get("phiBurninFlag") in (0, None)
    }

def main():
    """Command-line interface for image readiness checking"""
    if len(sys.argv) != 2:
        print("Usage: python image_readiness_checker.py <study_metadata.json>")
        sys.exit(1)
    
    try:
        with open(sys.argv[1], 'r') as f:
            study = json.load(f)
        
        report = readiness_report(study)
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
