# Pooling Strategy Calculator for Equimolar Sequencing

A Python tool for calculating optimal sub-pooling strategies for equimolar DNA library sequencing based on tapestation compact region table data.

## Recent Updates (August 2025)

- ✅ **Dual Molarity Unit Support**: Automatically handles both nmol/l and pmol/l inputs
- ✅ **Smart Plate Assignment**: Sequential plate numbering with clear console reporting
- ✅ **Organized Output Structure**: Automatic output folder creation with timestamped filenames
- ✅ **Improved Volume Constraints**: Updated minimum volume from 1.5μl to 3μl for strong samples
- ✅ **Enhanced File Management**: Automatic filtering of previous output files during processing

## Overview

This script processes compact region table CSV files to create optimal sub-pools for equimolar sequencing. It analyzes dimer and target library concentrations to calculate precise volumes needed for each sample, ensuring equimolar contribution while respecting volume and sample count constraints.

## Features

- **Automated data processing**: Handles multiple CSV files with different text encodings
- **Dual unit support**: Automatically detects and converts between nmol/l and pmol/l molarity units
- **Smart plate identification**: Sequential plate numbering with clear console reporting
- **Organized output**: Automatic creation of output folders with timestamped filenames
- **Dimer/Library classification**: Automatically identifies dimer (130-160bp) and target library (180-350bp) regions
- **Pool type segregation**: Separates strong (>5 nmol/l) and weak (≤5 nmol/l) samples into different sub-pools
- **Equimolar calculations**: Calculates exact volumes for equimolar contribution across samples
- **Constraint enforcement**: Respects volume limits (100-150μl per pool) and sample count limits
- **TECAN integration**: Direct compatibility with TECAN Freedom EVOware and Fluent Control
- **Liquid handling optimization**: Volume ranges optimized for automated pipetting systems
- **Comprehensive output**: Generates detailed CSV with all calculated parameters and TECAN-ready columns

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
- `Region Molarity [nmol/l]` or `Region Molarity [pmol/l]`: Molarity in nanomoles or picomoles per liter
- `% of Total`: Percentage of total
- `Region Comment`: Comments (optional)

### Molarity Unit Support

The script automatically detects and handles both molarity unit formats:
- **Primary**: `nmol/l` (nanomoles per liter) - used directly
- **Legacy**: `pmol/l` (picomoles per liter) - automatically converted to nmol/l (÷1000)

## File Organization

### Input Structure
```
your_data_folder/
├── Plate1_data.csv
├── Plate2_data.csv
└── PlateN_data.csv
```

### Output Structure (Auto-created)
```
your_data_folder/
├── Plate1_data.csv
├── Plate2_data.csv
├── PlateN_data.csv
└── output/
    ├── 2025-08-08_102509_sub-pooling.csv
    ├── 2025-08-08_143021_sub-pooling.csv
    └── ...
```

### Workflow Example

1. **Prepare your data folder**:
   ```
   my_experiment/
   ├── PlateA_results.csv
   └── PlateB_results.csv
   ```

2. **Run the script**:
   ```bash
   python pooling.py my_experiment
   ```

3. **Console output shows plate mapping**:
   ```
   Plate Assignment:
     Plate 001: PlateA_results.csv
     Plate 002: PlateB_results.csv
   ```

4. **Results are organized automatically**:
   ```
   my_experiment/
   ├── PlateA_results.csv
   ├── PlateB_results.csv
   └── output/
       └── 2025-08-08_143021_sub-pooling.csv
   ```

### Plate Assignment Logic

- **Single CSV**: Assigned to `SourcePlate[001]`
- **Multiple CSVs**: Sequential assignment (`001`, `002`, etc.) in alphabetical order
- **Console reporting**: Shows which file maps to each plate number
- **Automatic filtering**: Previous output files are ignored during processing

## Algorithm Logic

### Sample Classification

1. **Dimer identification**: Fragments with From ≤ 160bp and To ≤ 200bp
2. **Library identification**: Fragments with From ≥ 160bp
### Volume Constraints

**Strong pools** (>5 nmol/l): **3-7 μl per sample**
**Weak pools** (≤5 nmol/l): **7-20 μl per sample**
**Total pool volume**: **100-150 μl**
**Maximum samples per pool**: Configurable (default 48)

### Pooling Strategy

1. **Sort samples**: By library molarity (strongest to weakest)
2. **Initialize sub-pool**: Start with strongest available sample
3. **Add compatible samples**: Only samples of same pool type (strong/weak)
4. **Volume calculation**: Ensure equimolar contribution based on strongest sample
5. **Constraint checking**:
   - Volume per sample: 3-7μl (strong) or 7-20μl (weak)
   - Total pool volume: 100-150μl
   - Maximum samples per pool: Configurable (default 48)

### Quality Controls

- **Error detection**: Flags samples requiring <3μl or >20μl
- **Volume warnings**: Notes pools below 100μl minimum
- **Data validation**: Ensures no duplicate regions per sample
- **File filtering**: Automatically excludes previous output files from processing

## Output Format

The script generates timestamped CSV files in an `output` subfolder with the naming format: `YYYY-MM-DD_HHMMSS_sub-pooling.csv`

**Example**: `2025-08-08_102509_sub-pooling.csv`

### Output Columns

| Column | Description |
|--------|-------------|
| FileName | Original HSD1000 filename |
| Tape Well | Well identifier |
| Dimer Conc. | Dimer concentration [pg/μl] |
| Dimer Molarity | Dimer molarity [nmol/l] |
| Lib Conc. | Library concentration [pg/μl] |
| Lib Molarity | Library molarity [nmol/l] |
| target ratio | Target/dimer molarity ratio |
| sub-pool number | Assigned sub-pool number |
| volume added | Volume to add to sub-pool [μl] |
| target molarity contribution | Molar contribution to pool |
| sub-pool volume | Total sub-pool volume [μl] |
| sub-pool samples | Number of samples in sub-pool |
| notes | Warnings and flags |

### TECAN Liquid Handling Robot Columns

Additional columns formatted for TECAN liquid handling systems (compatible with Freedom EVOware and Fluent Control):

| Column | Description | TECAN Format |
|--------|-------------|--------------|
| SourcePlateLocation | Source plate labware identifier | "SourcePlate[001]", "SourcePlate[002]", etc. |
| SourceWellPosition | Source well position (1-96 format) | 1-96 (A1=1, A2=2, ..., H12=96) |
| VolSample | Sample volume to aspirate [μl] | Decimal format (e.g., 5.5, 12.3) |
| BufferLocation | Buffer plate labware identifier | "TEBuffer[001]" (constant) |
| BufferWellPosition | Buffer well position | 1 (constant - single buffer well) |
| VolBuffer | Buffer volume to add [μl] | 0.0 (no buffer dilution) |
| DestinationPlate | Destination plate labware identifier | "DestinationPlate[001]" (constant) |
| DestinationWellPosition | Destination well position | 1, 2, 3, ... (equals sub-pool number) |

#### TECAN Integration Notes

- **Labware Names**: Use exact labware names as defined in your TECAN worklist
- **Well Numbering**: Standard 96-well format (A1=1, A2=2, ..., A12=12, B1=13, ..., H12=96)
- **Volume Precision**: Volumes rounded to 0.1μl precision for pipetting accuracy
- **Plate Capacity**: Each destination well represents one sub-pool
- **Buffer Integration**: Ready for optional buffer addition workflows

#### Example TECAN Worklist Entry

```
Sample 1: SourcePlate[001], Well 1 (A1) → 5.5μl → DestinationPlate[001], Well 1
Sample 2: SourcePlate[001], Well 2 (A2) → 7.2μl → DestinationPlate[001], Well 1  
Sample 3: SourcePlate[002], Well 1 (A1) → 3.8μl → DestinationPlate[001], Well 2
```

## Example Output Summary

```
Processing CSV files in /data/results...

Plate Assignment:
  Plate 001: Sample_Batch_A.csv
  Plate 002: Sample_Batch_B.csv

Pooling strategy saved to: /data/results/output/2025-08-08_143021_sub-pooling.csv
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
- `Too strong - requires <3μl`: Sample molarity too high for practical pipetting
- `Too weak - requires >20μl`: Sample molarity too low for efficient pooling
- `Pool below 100μl minimum`: Sub-pool volume insufficient for handling

## Algorithm Parameters

| Parameter | Default | Description | TECAN Consideration |
|-----------|---------|-------------|-------------------|
| Strong pool threshold | 5 nmol/l | Molarity cutoff for strong vs weak classification | Optimized for standard tip volumes |
| Strong volume range | 3-7μl | Volume limits for strong samples | Within TECAN precision range |
| Weak volume range | 7-20μl | Volume limits for weak samples | Single tip capacity limit |
| Pool volume range | 100-150μl | Total volume limits per sub-pool | Optimal for downstream processing |
| Max samples per pool | 48 | Maximum number of samples per sub-pool | Half-plate processing efficiency |
| Minimum starting volume | 3μl | Minimum volume for strongest sample | Above TECAN dead volume |

## TECAN Integration

### Workflow Compatibility

This script generates output that is directly compatible with TECAN Freedom EVOware and Fluent Control software. The CSV output can be imported as a worklist or used to generate TECAN scripts.

### Supported TECAN Features

- **Multi-plate processing**: Handles samples from multiple source plates
- **Variable volume pipetting**: Optimized volumes for equimolar pooling
- **Plate barcode integration**: Ready for labware tracking systems
- **Quality control flagging**: Samples with volume issues clearly marked

### TECAN Worklist Import

1. **Load CSV**: Import the generated CSV file into EVOware/Fluent Control
2. **Map Labware**: Ensure labware names match your deck configuration:
   - Source plates: "SourcePlate[001]", "SourcePlate[002]", etc.
   - Destination plate: "DestinationPlate[001]"  
   - Buffer reservoir: "TEBuffer[001]"
3. **Configure Tips**: Use appropriate tip types for volume ranges (1.5-20μl)
4. **Run Protocol**: Execute automated pooling with liquid level detection

### Error Handling

Samples flagged in the 'notes' column require manual intervention:
- **"Too strong"**: Requires dilution before automated processing
- **"Too weak"**: May need concentration or manual pooling
- **"Pool below 100μl"**: Consider combining with another sub-pool


## Use Cases

- **High-throughput sequencing prep**: Optimize library pooling for Illumina sequencing
- **TECAN automation**: Generate worklists for Freedom EVOware and Fluent Control systems
- **Quality control**: Identify problematic samples before automated processing
- **Resource planning**: Calculate total volumes and pool requirements for liquid handlers
- **Automation integration**: Generate pipetting instructions compatible with robotic systems
- **Laboratory standardization**: Consistent pooling protocols across multiple TECAN instruments
