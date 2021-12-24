# pose_estimation
vision-based approaches to identifying key body parts and behavioral dynamics

# pipeline for validation, visualization, and sourcing of key training examples
and generated video annotations for key mistakes
1. use wrapper_validation_pipeline.m to start
2. play with sample datafiles in '/drosophila'
3. various auxiliary functions are distributed depending on class
4. generates plots of spatial statistics
5. interpolates low confidence values
6. generate re-annotated videos for visualization
7. generate heatmaps and trail plots for visualization
8. can then be used to visualize, identify statistical outliers, identify key error frames, and re-label recommended training examples for recalibrating training set

references for utilized tools:
- sleap: Pereira et al., 2020 (sleap.ai):
- repnan: www.ig.utexas.edu/people/students/cgreene/
- printStruct
