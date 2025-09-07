"""
Centralized configuration for classification confidence thresholds and CNV analysis rules.

This module provides a single source of truth for confidence thresholds
and CNV classification rules used across the robin application (GUI, reporting, etc.).
"""

from typing import Dict, Tuple

# Confidence thresholds for different classifiers
CLASSIFIER_CONFIDENCE_THRESHOLDS: Dict[str, Dict[str, float]] = {
    "sturgeon": {
        "high": 85.0,
        "medium": 65.0,
        "low": 0.0,
    },
    "nanodx": {
        "high": 50.0,
        "medium": 25.0,
        "low": 0.0,
    },
    "pannanodx": {
        "high": 50.0,
        "medium": 25.0,
        "low": 0.0,
    },
    "random_forest": {
        "high": 85.0,
        "medium": 65.0,
        "low": 0.0,
    },
}

# Default thresholds (fallback)
DEFAULT_CONFIDENCE_THRESHOLDS = {
    "high": 80.0,
    "medium": 50.0,
    "low": 20.0,
}

def get_confidence_level(classifier: str, confidence: float) -> str:
    """
    Get confidence level based on classifier-specific thresholds.
    
    Args:
        classifier: Name of the classifier (e.g., 'sturgeon', 'nanodx')
        confidence: Confidence value as percentage (0-100)
    
    Returns:
        Confidence level string ('High confidence', 'Medium confidence', etc.)
    """
    thresholds = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )
    
    if confidence >= thresholds["high"]:
        return "High confidence"
    elif confidence >= thresholds["medium"]:
        return "Medium confidence"
    elif confidence >= thresholds["low"]:
        return "Low confidence"
    else:
        return "Very low confidence"

def get_confidence_status(classifier: str, confidence: float) -> tuple[str, str]:
    """
    Get confidence status and color for reporting.
    
    Args:
        classifier: Name of the classifier (e.g., 'sturgeon', 'nanodx')
        confidence: Confidence value as decimal (0-1)
    
    Returns:
        Tuple of (status, color_hex) where status is 'High', 'Medium', or 'Low'
    """
    # Convert to percentage for threshold comparison
    confidence_percent = confidence * 100
    
    thresholds = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )
    
    if confidence_percent >= thresholds["high"]:
        return "High", "#059669"  # Green
    elif confidence_percent >= thresholds["medium"]:
        return "Medium", "#D97706"  # Amber
    else:
        return "Low", "#DC2626"  # Red

# CNV Analysis Configuration
CNV_THRESHOLDS = {
    "autosomes": {
        "gain": 0.4,
        "loss": -0.4,
    },
    "chrX": {
        "male": {
            "gain": 0.3,
            "loss": -0.3,
        },
        "female": {
            "gain": 0.75,
            "loss": -0.75,
        }
    },
    "chrY": {
        "male": {
            "gain": 0.5,
            "loss": -0.5,
        },
        "female": {
            "gain": 0.2,  # Fixed: positive threshold for gains
            "loss": -1.0,
        }
    }
}

# CNV Event Detection Rules
CNV_EVENT_RULES = {
    "whole_chromosome": {
        "min_proportion_affected": 0.7,  # 70% of bins must support the event
        "min_arm_proportion": 0.4,  # Each arm must show at least 40% effect
        "single_arm_multiplier": 1.5,  # For single-arm chromosomes, use 1.5x threshold
    },
    "arm_specific": {
        "min_proportion_affected": 0.7,  # 70% of arm must be affected
    },
    "resolution": {
        "max_bin_width": 10_000_000,  # 10Mb - resolution too low for accurate CNV calling
    }
}

def get_cnv_thresholds(chromosome: str, sex_estimate: str) -> Tuple[float, float]:
    """
    Get CNV gain/loss thresholds for a specific chromosome and sex.
    
    Args:
        chromosome: Chromosome name (e.g., 'chr1', 'chrX', 'chrY')
        sex_estimate: Sex estimate ('XY', 'XX', 'Male', 'Female')
    
    Returns:
        Tuple of (gain_threshold, loss_threshold)
    """
    # Normalize sex estimate
    sex_key = "male" if sex_estimate.upper() in ("XY", "MALE") else "female"
    
    if chromosome == "chrX":
        thresholds = CNV_THRESHOLDS["chrX"][sex_key]
    elif chromosome == "chrY":
        thresholds = CNV_THRESHOLDS["chrY"][sex_key]
    else:
        # Autosomes
        thresholds = CNV_THRESHOLDS["autosomes"]
    
    return thresholds["gain"], thresholds["loss"]

def is_whole_chromosome_event(
    p_arm_mean: float, 
    q_arm_mean: float, 
    p_arm_proportion: float, 
    q_arm_proportion: float,
    gain_threshold: float,
    loss_threshold: float
) -> Tuple[bool, str]:
    """
    Determine if a chromosome shows a whole chromosome event.
    
    Args:
        p_arm_mean: Mean CNV value for p-arm
        q_arm_mean: Mean CNV value for q-arm  
        p_arm_proportion: Proportion of p-arm affected
        q_arm_proportion: Proportion of q-arm affected
        gain_threshold: Gain threshold for this chromosome
        loss_threshold: Loss threshold for this chromosome
    
    Returns:
        Tuple of (is_whole_chr_event, event_type) where event_type is 'GAIN', 'LOSS', or 'NORMAL'
    """
    rules = CNV_EVENT_RULES["whole_chromosome"]
    
    # Check for gains - both arms must show gains, but with more nuanced proportion requirements
    both_arms_gained = (
        p_arm_mean > gain_threshold
        and q_arm_mean > gain_threshold
        and (p_arm_proportion > rules["min_proportion_affected"] or q_arm_proportion > rules["min_proportion_affected"])  # At least one arm > 70%
        and (p_arm_proportion > rules["min_arm_proportion"] and q_arm_proportion > rules["min_arm_proportion"])  # Both arms > 40%
    )
    
    # Check for losses - both arms must show losses, but with more nuanced proportion requirements
    both_arms_lost = (
        p_arm_mean < loss_threshold
        and q_arm_mean < loss_threshold
        and (p_arm_proportion > rules["min_proportion_affected"] or q_arm_proportion > rules["min_proportion_affected"])  # At least one arm > 70%
        and (p_arm_proportion > rules["min_arm_proportion"] and q_arm_proportion > rules["min_arm_proportion"])  # Both arms > 40%
    )
    
    if both_arms_gained:
        return True, "GAIN"
    elif both_arms_lost:
        return True, "LOSS"
    else:
        return False, "NORMAL"

def is_arm_event(
    arm_mean: float,
    arm_proportion: float,
    gain_threshold: float,
    loss_threshold: float
) -> Tuple[bool, str]:
    """
    Determine if a chromosome arm shows a significant event.
    
    Args:
        arm_mean: Mean CNV value for the arm
        arm_proportion: Proportion of arm affected
        gain_threshold: Gain threshold for this chromosome
        loss_threshold: Loss threshold for this chromosome
    
    Returns:
        Tuple of (is_arm_event, event_type) where event_type is 'GAIN', 'LOSS', or 'NORMAL'
    """
    rules = CNV_EVENT_RULES["arm_specific"]
    
    if arm_proportion > rules["min_proportion_affected"]:
        if arm_mean > gain_threshold:
            return True, "GAIN"
        elif arm_mean < loss_threshold:
            return True, "LOSS"
    
    return False, "NORMAL"

def is_resolution_sufficient(bin_width: int) -> bool:
    """
    Check if the resolution is sufficient for CNV calling.
    
    Args:
        bin_width: Bin width in base pairs
    
    Returns:
        True if resolution is sufficient for CNV calling
    """
    return bin_width <= CNV_EVENT_RULES["resolution"]["max_bin_width"]
