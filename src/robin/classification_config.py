"""
Centralized configuration for classification confidence thresholds.

This module provides a single source of truth for confidence thresholds
used across the robin application (GUI, reporting, etc.).
"""

from typing import Dict

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

def get_classifier_thresholds(classifier: str) -> Dict[str, float]:
    """
    Get the confidence thresholds for a specific classifier.
    
    Args:
        classifier: Name of the classifier
    
    Returns:
        Dictionary with 'high', 'medium', 'low' threshold values
    """
    return CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier, DEFAULT_CONFIDENCE_THRESHOLDS
    )
