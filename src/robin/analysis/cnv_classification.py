"""
Centralized CNV classification and analysis functions.

This module provides CNV analysis functions that use the centralized
classification rules from classification_config.py.
"""

import logging
import numpy as np
import pandas as pd
import natsort
from typing import Dict, List, Tuple, Optional, Any

from robin.classification_config import (
    get_cnv_thresholds,
    is_whole_chromosome_event,
    is_arm_event,
    is_resolution_sufficient
)

logger = logging.getLogger(__name__)


class CNVEvent:
    """Represents a CNV event with metadata."""
    
    def __init__(
        self,
        chromosome: str,
        event_type: str,
        mean_cnv: float,
        start_pos: int,
        end_pos: int,
        length: int,
        genes: List[str] = None,
        confidence: str = "Unknown",
        arm: Optional[str] = None,
        proportion_affected: float = 0.0
    ):
        self.chromosome = chromosome
        self.event_type = event_type  # 'GAIN', 'LOSS', 'WHOLE_CHR_GAIN', 'WHOLE_CHR_LOSS'
        self.mean_cnv = mean_cnv
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.length = length
        self.genes = genes or []
        self.confidence = confidence
        self.arm = arm
        self.proportion_affected = proportion_affected
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for GUI/reporting."""
        return {
            "chromosome": self.chromosome,
            "event_type": self.event_type,
            "mean_cnv": self.mean_cnv,
            "start_pos": self.start_pos,
            "end_pos": self.end_pos,
            "length": self.length,
            "genes": self.genes,
            "confidence": self.confidence,
            "arm": self.arm,
            "proportion_affected": self.proportion_affected,
            "start_mb": f"{self.start_pos/1e6:.2f}",
            "end_mb": f"{self.end_pos/1e6:.2f}",
            "length_mb": f"{self.length/1e6:.2f}",
            "mean_cnv_str": f"{self.mean_cnv:.3f}",
            "genes_str": ", ".join(self.genes) if self.genes else "",
        }


def analyze_chromosome_arms(
    cnv_data: Dict[str, np.ndarray],
    chromosome: str,
    bin_width: int,
    sex_estimate: str,
    cytobands_df: pd.DataFrame
) -> Tuple[Optional[float], Optional[float], float, float]:
    """
    Analyze p and q arms of a chromosome for CNV events.
    
    Args:
        cnv_data: CNV data dictionary
        chromosome: Chromosome to analyze
        bin_width: Bin width in base pairs
        sex_estimate: Sex estimate
        cytobands_df: Cytobands dataframe
    
    Returns:
        Tuple of (p_arm_mean, q_arm_mean, p_arm_proportion, q_arm_proportion)
    """
    if chromosome not in cnv_data:
        return None, None, 0.0, 0.0
    
    # Get thresholds for this chromosome
    gain_threshold, loss_threshold = get_cnv_thresholds(chromosome, sex_estimate)
    
    # Get chromosome cytobands
    chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
    if chr_cytobands.empty:
        return None, None, 0.0, 0.0
    
    # Debug: log cytoband names for troubleshooting
    logger.debug(f"{chromosome} cytoband names: {chr_cytobands['name'].tolist()}")
    
    # Analyze p arm - cytoband names start with 'p' (e.g., p36.33, p36.32)
    p_arm_cytobands = chr_cytobands[
        chr_cytobands["name"].str.startswith("p", na=False)
    ]
    
    p_arm_mean = None
    p_arm_proportion = 0.0
    
    if not p_arm_cytobands.empty:
        p_arm_values = []
        for _, band in p_arm_cytobands.iterrows():
            start_bin = max(0, int(band["start"] // bin_width))
            end_bin = min(len(cnv_data[chromosome]) - 1, int(band["end"] // bin_width))
            if end_bin >= start_bin:
                region_values = cnv_data[chromosome][start_bin:end_bin + 1]
                p_arm_values.extend(region_values)
        
        if p_arm_values:
            p_arm_mean = float(np.mean(p_arm_values))
            p_arm_proportion_gain = sum(1 for v in p_arm_values if v > gain_threshold) / len(p_arm_values)
            p_arm_proportion_loss = sum(1 for v in p_arm_values if v < loss_threshold) / len(p_arm_values)
            p_arm_proportion = max(p_arm_proportion_gain, p_arm_proportion_loss)
    
    logger.debug(f"{chromosome} p-arm: {len(p_arm_cytobands)} bands found, mean={p_arm_mean}, prop={p_arm_proportion:.3f}")
    
    # Analyze q arm - cytoband names start with 'q' (e.g., q11, q12, q21.1)
    q_arm_cytobands = chr_cytobands[
        chr_cytobands["name"].str.startswith("q", na=False)
    ]
    
    q_arm_mean = None
    q_arm_proportion = 0.0
    
    if not q_arm_cytobands.empty:
        q_arm_values = []
        for _, band in q_arm_cytobands.iterrows():
            start_bin = max(0, int(band["start"] // bin_width))
            end_bin = min(len(cnv_data[chromosome]) - 1, int(band["end"] // bin_width))
            if end_bin >= start_bin:
                region_values = cnv_data[chromosome][start_bin:end_bin + 1]
                q_arm_values.extend(region_values)
        
        if q_arm_values:
            q_arm_mean = float(np.mean(q_arm_values))
            q_arm_proportion_gain = sum(1 for v in q_arm_values if v > gain_threshold) / len(q_arm_values)
            q_arm_proportion_loss = sum(1 for v in q_arm_values if v < loss_threshold) / len(q_arm_values)
            q_arm_proportion = max(q_arm_proportion_gain, q_arm_proportion_loss)
    
    logger.debug(f"{chromosome} q-arm: {len(q_arm_cytobands)} bands found, mean={q_arm_mean}, prop={q_arm_proportion:.3f}")
    
    return p_arm_mean, q_arm_mean, p_arm_proportion, q_arm_proportion


def detect_cnv_events(
    cnv_data: Dict[str, np.ndarray],
    bin_width: int,
    sex_estimate: str,
    cytobands_df: pd.DataFrame,
    gene_df: Optional[pd.DataFrame] = None
) -> List[CNVEvent]:
    """
    Detect CNV events using centralized classification rules.
    
    Args:
        cnv_data: CNV data dictionary
        bin_width: Bin width in base pairs
        sex_estimate: Sex estimate
        cytobands_df: Cytobands dataframe
        gene_df: Optional gene dataframe
    
    Returns:
        List of CNVEvent objects
    """
    events = []
    
    logger.debug(f"Detecting CNV events with sex_estimate='{sex_estimate}', bin_width={bin_width}")
    
    # Check if resolution is sufficient
    if not is_resolution_sufficient(bin_width):
        logger.warning(f"Resolution insufficient for CNV calling: bin_width={bin_width}")
        return events
    
    # Analyze each chromosome
    logger.debug(f"Available chromosomes: {list(cnv_data.keys())}")
    for chromosome in natsort.natsorted(cnv_data.keys()):
        if chromosome == "chrM" or not chromosome.startswith("chr"):
            continue
        
        # Skip Y chromosome for male samples (expected absence)
        if chromosome == "chrY" and sex_estimate.upper() in ("XY", "MALE"):
            logger.debug(f"Skipping {chromosome} for male sample")
            continue
        
        logger.debug(f"Analyzing chromosome {chromosome} for CNV events")
        
        # Get thresholds
        gain_threshold, loss_threshold = get_cnv_thresholds(chromosome, sex_estimate)
        
        # Analyze arms
        p_arm_mean, q_arm_mean, p_arm_proportion, q_arm_proportion = analyze_chromosome_arms(
            cnv_data, chromosome, bin_width, sex_estimate, cytobands_df
        )
        
        # Check for whole chromosome events
        if p_arm_mean is not None and q_arm_mean is not None:
            logger.debug(f"{chromosome}: p_arm_mean={p_arm_mean:.3f}, q_arm_mean={q_arm_mean:.3f}, p_prop={p_arm_proportion:.3f}, q_prop={q_arm_proportion:.3f}")
            is_whole_chr, event_type = is_whole_chromosome_event(
                p_arm_mean, q_arm_mean, p_arm_proportion, q_arm_proportion,
                gain_threshold, loss_threshold
            )
            
            if is_whole_chr:
                # Create whole chromosome event
                chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                if not chr_cytobands.empty:
                    start_pos = int(chr_cytobands["start"].min())
                    end_pos = int(chr_cytobands["end"].max())
                    length = end_pos - start_pos
                    
                    # Get genes in chromosome
                    genes = []
                    if gene_df is not None:
                        genes = gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                    
                    # Calculate overall chromosome mean
                    chr_mean = float(np.mean(cnv_data[chromosome]))
                    
                    event = CNVEvent(
                        chromosome=chromosome,
                        event_type=f"WHOLE_CHR_{event_type}",
                        mean_cnv=chr_mean,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        length=length,
                        genes=genes,
                        confidence="High" if min(p_arm_proportion, q_arm_proportion) > 0.8 else "Medium",
                        proportion_affected=max(p_arm_proportion, q_arm_proportion)
                    )
                    events.append(event)
                    logger.info(f"Detected whole chromosome {event_type} for {chromosome}")
            
            # Always check for individual arm events, even if whole chromosome event was detected
            # This ensures we capture significant regional variations
            if chromosome != "chrY":  # Skip Y chromosome arm events
                # Check p arm
                if p_arm_mean is not None:
                    is_p_event, p_event_type = is_arm_event(
                        p_arm_mean, p_arm_proportion, gain_threshold, loss_threshold
                    )
                    if is_p_event:
                        chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                        p_bands = chr_cytobands[chr_cytobands["name"].str.startswith("p", na=False)]
                        if not p_bands.empty:
                            start_pos = int(p_bands["start"].min())
                            end_pos = int(p_bands["end"].max())
                            length = end_pos - start_pos
                            
                            # Get genes in p arm
                            genes = []
                            if gene_df is not None:
                                genes = gene_df[
                                    (gene_df["chrom"] == chromosome) &
                                    (gene_df["start"] <= end_pos) &
                                    (gene_df["end"] >= start_pos)
                                ]["gene"].astype(str).tolist()
                            
                            event = CNVEvent(
                                chromosome=chromosome,
                                event_type=p_event_type,
                                mean_cnv=p_arm_mean,
                                start_pos=start_pos,
                                end_pos=end_pos,
                                length=length,
                                genes=genes,
                                confidence="High" if p_arm_proportion > 0.8 else "Medium",
                                arm="p",
                                proportion_affected=p_arm_proportion
                            )
                            events.append(event)
                            logger.info(f"Detected p-arm {p_event_type} for {chromosome}")
                
                # Check q arm
                if q_arm_mean is not None:
                    is_q_event, q_event_type = is_arm_event(
                        q_arm_mean, q_arm_proportion, gain_threshold, loss_threshold
                    )
                    if is_q_event:
                        chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                        q_bands = chr_cytobands[chr_cytobands["name"].str.startswith("q", na=False)]
                        if not q_bands.empty:
                            start_pos = int(q_bands["start"].min())
                            end_pos = int(q_bands["end"].max())
                            length = end_pos - start_pos
                            
                            # Get genes in q arm
                            genes = []
                            if gene_df is not None:
                                genes = gene_df[
                                    (gene_df["chrom"] == chromosome) &
                                    (gene_df["start"] <= end_pos) &
                                    (gene_df["end"] >= start_pos)
                                ]["gene"].astype(str).tolist()
                            
                            event = CNVEvent(
                                chromosome=chromosome,
                                event_type=q_event_type,
                                mean_cnv=q_arm_mean,
                                start_pos=start_pos,
                                end_pos=end_pos,
                                length=length,
                                genes=genes,
                                confidence="High" if q_arm_proportion > 0.8 else "Medium",
                                arm="q",
                                proportion_affected=q_arm_proportion
                            )
                            events.append(event)
                            logger.info(f"Detected q-arm {q_event_type} for {chromosome}")
        else:
            # Single arm chromosome - use stricter threshold
            # Only use this logic if we truly have only one arm (like chrY in some cases)
            # For chromosomes that should have both arms, this indicates a cytoband parsing issue
            logger.warning(f"{chromosome}: Only one arm detected - this may indicate a cytoband parsing issue")
            
            whole_chr_mean = float(np.mean(cnv_data[chromosome]))
            single_arm_multiplier = 1.5  # From CNV_EVENT_RULES
            
            if abs(whole_chr_mean) > abs(gain_threshold) * single_arm_multiplier:
                chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                if not chr_cytobands.empty:
                    start_pos = int(chr_cytobands["start"].min())
                    end_pos = int(chr_cytobands["end"].max())
                    length = end_pos - start_pos
                    
                    # Get genes in chromosome
                    genes = []
                    if gene_df is not None:
                        genes = gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                    
                    event_type = "GAIN" if whole_chr_mean > gain_threshold else "LOSS"
                    event = CNVEvent(
                        chromosome=chromosome,
                        event_type=f"WHOLE_CHR_{event_type}",
                        mean_cnv=whole_chr_mean,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        length=length,
                        genes=genes,
                        confidence="Medium",  # Single arm events are less certain
                        proportion_affected=1.0  # Entire chromosome
                    )
                    events.append(event)
                    logger.info(f"Detected single-arm whole chromosome {event_type} for {chromosome}")
    
    return events


def get_cnv_summary(events: List[CNVEvent]) -> Dict[str, Any]:
    """
    Generate a summary of CNV events.
    
    Args:
        events: List of CNVEvent objects
    
    Returns:
        Dictionary with summary statistics
    """
    summary = {
        "total_events": len(events),
        "whole_chromosome_events": [],
        "arm_events": [],
        "gene_containing_events": [],
        "total_genes_affected": set(),
        "high_confidence_events": 0,
        "medium_confidence_events": 0,
    }
    
    for event in events:
        if event.event_type.startswith("WHOLE_CHR_"):
            summary["whole_chromosome_events"].append(event)
        else:
            summary["arm_events"].append(event)
        
        if event.genes:
            summary["gene_containing_events"].append(event)
            summary["total_genes_affected"].update(event.genes)
        
        if event.confidence == "High":
            summary["high_confidence_events"] += 1
        elif event.confidence == "Medium":
            summary["medium_confidence_events"] += 1
    
    summary["total_genes_affected"] = len(summary["total_genes_affected"])
    
    return summary
