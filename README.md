# pose estimation:
- Lightweight CNN-based approaches to inference of body nodes and connectivity-based identity tracking
- Example analysis of drosophila in close social contact, kinetics, and spatial relationships
- Custom pipeline for parallelization, visualization, statistical validation, and auto-sourcing examples, resulting in >10x efficiency in model development, improved performance (from 80% to >98%), and 3-4x reduction in training set examples
- Systematic model development as above, blending of top-down and bottom-up models, and other innovations enabled markerless disambiguation of occlusions and made model performance robust in interesting but difficult to parse social situations, e.g. aggression
- Enabled high-resolution analysis of behavioral dynamics of mice in close social contact (manuscript in preparation)

## pipeline for validation, visualization, and sourcing of key training examples:
1. use wrapper_validation_pipeline.m to start
2. play with sample datafiles in '/drosophila'
3. various auxiliary functions are distributed depending on class
4. generates plots of spatial statistics
5. interpolates low confidence values
6. generate re-annotated videos for visualization of key mistakes
7. generate heatmaps and trail plots for visualization
8. visualize, identify statistical outliers, identify key error frames, and re-label recommended training examples for recalibrating training set

## example inferenced videos:
<p align="center">
<img src="example_movies_processed/PairS1_20190128_113421.120s_ann_samplecrop2.gif">
</p>

## example analysis tools:
<p align="center">
<img src="example_figures/fig5.png">
</p>

<p align="center">
<img src="example_figures/hm_1-1.gif">
</p>

## acknowledgements:
- sleap: Pereira et al., 2020 (sleap.ai)
- deepposekit: Graving et al., 2019 (deeposekit.org)
- repnan: www.ig.utexas.edu/people/students/cgreene/
- printStruct
- avi2gif (Lindner)
