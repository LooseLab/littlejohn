"""
Configuration and helper functions for GUI section visibility.
"""

from typing import List, Set, Optional

# Re-export get_confidence_level from classification_config for backward compatibility
try:
    from robin.classification_config import get_confidence_level
except ImportError:
    # Fallback if classification_config is not available
    def get_confidence_level(classifier: str, confidence: float) -> str:
        """Fallback confidence level function."""
        if confidence >= 80:
            return "High confidence"
        elif confidence >= 50:
            return "Medium confidence"
        elif confidence >= 20:
            return "Low confidence"
        else:
            return "Very low confidence"


# Map workflow step names to GUI section names
WORKFLOW_STEP_TO_SECTION = {
    "target": "target",
    "mgmt": "mgmt",
    "fusion": "fusion",
    "cnv": "cnv",
    "sturgeon": "sturgeon",
    "nanodx": "nanodx",
    "random_forest": "random_forest",
    "pannanodx": "pannanodx",
}

# Map workflow steps to their display names in classification section
CLASSIFICATION_STEPS = {
    "sturgeon": "Sturgeon",
    "nanodx": "NanoDX",
    "random_forest": "Random Forest",
    "pannanodx": "PanNanoDX",
}


def get_enabled_sections(workflow_steps: Optional[List[str]]) -> Set[str]:
    """
    Determine which sections should be enabled based on workflow steps.
    
    Args:
        workflow_steps: List of workflow step names (e.g., ['target', 'mgmt', 'sturgeon'])
        
    Returns:
        Set of enabled section names
    """
    if not workflow_steps:
        # If no workflow steps specified, show all sections (backward compatibility)
        return set(WORKFLOW_STEP_TO_SECTION.values())
    
    enabled = set()
    for step in workflow_steps:
        # Handle workflow steps that might have queue prefixes (e.g., "classification:sturgeon")
        step_name = step.split(":")[-1] if ":" in step else step
        if step_name in WORKFLOW_STEP_TO_SECTION:
            enabled.add(WORKFLOW_STEP_TO_SECTION[step_name])
    
    return enabled


def is_section_enabled(section_name: str, workflow_steps: Optional[List[str]]) -> bool:
    """
    Check if a specific section should be enabled.
    
    Args:
        section_name: Name of the section to check
        workflow_steps: List of workflow step names
        
    Returns:
        True if section should be enabled, False otherwise
    """
    enabled_sections = get_enabled_sections(workflow_steps)
    
    # If no workflow steps specified, show all sections (backward compatibility)
    if not workflow_steps:
        return True
    
    return section_name in enabled_sections


def get_enabled_classification_steps(workflow_steps: Optional[List[str]]) -> Set[str]:
    """
    Get the set of enabled classification steps.
    
    Args:
        workflow_steps: List of workflow step names
        
    Returns:
        Set of enabled classification step names (e.g., {'sturgeon', 'nanodx'})
    """
    enabled_sections = get_enabled_sections(workflow_steps)
    return enabled_sections.intersection(set(CLASSIFICATION_STEPS.keys()))
