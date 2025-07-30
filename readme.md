# Pooling Strategy Calculator for Equimolar Sequencing

A Python tool for calculating optimal sub-pooling strategies for equimolar DNA library sequencing based on tapestation compact region table data.

## Overview

This script processes compact region table CSV files to create optimal sub-pools for equimolar sequencing. It analyzes dimer and target library concentrations to calculate precise volumes needed for each sample, ensuring equimolar contribution while respecting volume and sample count constraints.

## Features

- **Automated data processing**: Handles multiple CSV files with different text encodings
- **Dimer/Library classification**: Automatically identifies dimer (130-160bp) and target library (180-350bp) regions
- **Pool type segregation**: Separates strong (>5 nmol/l) and weak (≤5 nmol/l) samples into different sub-pools
- **Equimolar calculations**: Calculates exact volumes for equimolar contribution across samples
- **Constraint enforcement**: Respects volume limits (100-150μl per pool) and sample count limits
- **Comprehensive output**: Generates detailed CSV with all calculated parameters and flags

## Requirements

- Python 3.6+
- pandas
- argparse (standard library)

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd pooling-tapes

# Install dependencies
pip install pandas
```

## Usage

### Basic Usage

```bash
python pooling.py /path/to/csv/folder
```

### Advanced Usage

```bash
# Specify maximum samples per sub-pool (default: 48)
python pooling.py /path/to/csv/folder --max-samples 24

# Example with full path
python sub-pooling.py "C:/data/tape_results" --max-samples 96
```

## Input Data Format

The script expects D1000 compact region table CSV files with the following columns:

- `FileName`: Name of the D1000 analysis file
- `WellId`: Well identifier (e.g., A1, B2, etc.)
- `Sample Description`: Sample description (optional)
- `From [bp]`: Start position in base pairs
- `To [bp]`: End position in base pairs
- `Average Size [bp]`: Average fragment size
- `Conc. [pg/μl]`: Concentration in picograms per microliter
- `Region Molarity [pmol/l]`: Molarity in picomoles per liter
- `% of Total`: Percentage of total
- `Region Comment`: Comments (optional)

## Algorithm Logic

### Sample Classification

1. **Dimer identification**: Fragments with From ≤ 160bp and To ≤ 200bp
2. **Library identification**: Fragments with From ≥ 160bp
3. **Pool type determination**: 
   - Strong pools: Library molarity > 5 nmol/l (1.5-7μl per sample)
   - Weak pools: Library molarity ≤ 5 nmol/l (7-20μl per sample)

### Pooling Strategy

1. **Sort samples**: By library molarity (strongest to weakest)
2. **Initialize sub-pool**: Start with strongest available sample
3. **Add compatible samples**: Only samples of same pool type (strong/weak)
4. **Volume calculation**: Ensure equimolar contribution based on strongest sample
5. **Constraint checking**:
   - Volume per sample: 1.5-7μl (strong) or 7-20μl (weak)
   - Total pool volume: 100-150μl
   - Maximum samples per pool: Configurable (default 48)

### Quality Controls

- **Error detection**: Flags samples requiring <1μl or >20μl
- **Volume warnings**: Notes pools below 100μl minimum
- **Data validation**: Ensures no duplicate regions per sample

## Output Format

The script generates a CSV file named `YYYY-MM-DD_sub-pooling.csv` with the following columns:

| Column | Description |
|--------|-------------|
| FileName | Original HSD1000 filename |
| Tape Well | Well identifier |
| Dimer Conc. | Dimer concentration [pg/μl] |
| Dimer Molarity | Dimer molarity [pmol/l] |
| Lib Conc. | Library concentration [pg/μl] |
| Lib Molarity | Library molarity [pmol/l] |
| target ratio | Target/dimer molarity ratio |
| sub-pool number | Assigned sub-pool number |
| volume added | Volume to add to sub-pool [μl] |
| target molarity contribution | Molar contribution to pool |
| sub-pool volume | Total sub-pool volume [μl] |
| sub-pool samples | Number of samples in sub-pool |
| notes | Warnings and flags |

## Example Output Summary

```
Processing CSV files in /data/results...
Pooling strategy saved to: /data/results/2025-07-30_sub-pooling.csv
Processed 192 samples into 7 sub-pools

Sub-pool Summary:
Pool 1: 48 samples, 141.7μl total
Pool 2: 48 samples, 85.5μl total
Pool 3: 38 samples, 83.3μl total
Pool 4: 19 samples, 149.9μl total
Pool 5: 19 samples, 145.4μl total
Pool 6: 17 samples, 140.4μl total
Pool 7: 3 samples, 22.2μl total
```

## Troubleshooting

### Common Issues

1. **Encoding errors**: Script automatically tries multiple encodings (UTF-8, Latin-1, CP1252, ISO-8859-1)
2. **Missing columns**: Check that CSV files contain required concentration and molarity columns
3. **Multiple regions**: Error if sample has multiple dimer or library regions
4. **Empty results**: Verify CSV files contain valid data and proper column headers

### Error Messages

- `Multiple dimer regions found`: Sample has more than one dimer region
- `Too strong - requires <1μl`: Sample molarity too high for practical pipetting
- `Too weak - requires >20μl`: Sample molarity too low for efficient pooling
- `Pool below 100μl minimum`: Sub-pool volume insufficient for handling

## Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Strong pool threshold | 5 nmol/l | Molarity cutoff for strong vs weak classification |
| Strong volume range | 1.5-7μl | Volume limits for strong samples |
| Weak volume range | 7-20μl | Volume limits for weak samples |
| Pool volume range | 100-150μl | Total volume limits per sub-pool |
| Max samples per pool | 48 | Maximum number of samples per sub-pool |
| Minimum starting volume | 1.5μl | Minimum volume for strongest sample |

## Use Cases

- **High-throughput sequencing prep**: Optimize library pooling for Illumina sequencing
- **Quality control**: Identify problematic samples before sequencing
- **Resource planning**: Calculate total volumes and pool requirements
- **Automation integration**: Generate pipetting instructions for liquid handlers
