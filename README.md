# GC-based nucleosome prediction

Uses GC content of a genome to predict nucleosome occupancy. Outputs the positions to a tab-delimited file with GC-smoothed values based on GC counts. This is done by binarizing the genome to 1 where C/G and 0 where A/T and then applying a Gaussian convolution to these binary values.


## Usage

```
gc_nucleosome_predictor.py [-hg] <Genome FASTA File>

Options:
  -h, --help            show this help message and exit
  -g GAUSSIAN_INTERVAL  Interval length for Gaussian smoother; in bp) [138]
  --no-secondary-smoothing
                        Turn off SMA (Simple Moving Average) to smooth
                        afterwards for less jagged graph
```

(Previously hosted at http://baderlab.org/Software/nucleosome-prediction)