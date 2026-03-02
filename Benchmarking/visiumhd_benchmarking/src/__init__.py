# Visium HD H&E Morphology Benchmarking
"""
Benchmarking framework for H&E morphology-based cell type assignment
using pre-trained ViT features and proportion-guided MIL.

Modules:
    create_pseudo_visium: Generate pseudo-Visium spots from Visium HD data
    run_cellpose_he: Cellpose segmentation wrapper for H&E images
    extract_patches_he: Extract 224x224 H&E patches per nucleus
    vit_extractor: Pre-trained ViT feature extraction (UNI model)
    proportion_mil: Proportion-guided MIL aggregation
    train_mil: Training pipeline for MIL head
    evaluate_single_cell: Single-cell assignment evaluation
    run_benchmark: Main benchmark pipeline
"""
