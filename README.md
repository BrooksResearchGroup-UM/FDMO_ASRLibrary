# Flavin Dependent Monooxygenase Ancestral Library Creation 

Computational pipeline supporting the Gen3 FDMO ancestral sequence reconstruction library paper.

## Repository Structure

### pre_ancestors/
Scripts and models for the original FDMO latent space analysis (Narayan Lab VAE, no ancestors).
- `latent_space/` — VAE training, 2D latent space model
- `clustering/` — K-means clustering (k=40), elbow method
- `scripts/` — supporting utilities

### post_ancestors/
Scripts for the ancestral sequence reconstruction pipeline and updated latent space analysis.
- `ASR_pipeline/` — MUSCLE5 alignment, indel filtering, TrimAl, PhyML, PyASR, post-processing, filtering
- `latent_space/` — retrained VAE with 276 Gen3 ancestors injected
- `analysis/` — reactivity data processing
- `data/` — plate layouts, reactivity CSVs

## Dependencies
- Python 3.11
- PyTorch, scikit-learn, BioPython, pandas, numpy, matplotlib, plotly
- MUSCLE5, TrimAl, PhyML, ModelTest-NG, PyASR, PAML, HMMER

## Model
The latent space model used for figures is a 2D VAE trained on 33,972 FDMO sequences (352 aa)
with 276 Gen3 ancestors injected. Model checkpoint in `post_ancestors/latent_space/models/`.

## HMMER Installation 
1. Download HMMER from the official [HMMER website](http://hmmer.org/download.html).
2. Follow the installation instructions provided on the website or in the downloaded package.
3. To verify that HMMER has been installed successfully, run the following command in your terminal or command line:
   ```bash
   hmmalign --version
## MUSCLE5 Installation 
1. Download from the MUSCLE github page (https://github.com/rcedgar/muscle/releases/tag/v5.3)
   ```bash
   MUSCLE5 --version

## PhyML Installation 
1. Download from the PhyML github page (https://github.com/stephaneguindon/phyml)

## PyASR Installation 
1. Go to the PyASR webpage (http://pypi.org/project/pyasr/)
2. Copy the pip install command from the pypi page.
   ```bash
   pip install pyasr

