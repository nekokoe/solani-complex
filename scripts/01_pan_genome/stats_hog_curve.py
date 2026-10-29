#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HOG Curve Statistics Data Calculation

This script reads a TSV file containing SampleID, HOG_ID, and Count columns,
and performs the following statistical analysis:
1. Generate pan-genome saturation curve data
2. Generate HOG count change data for each presence rate interval
3. Save data as TSV files

Usage:
    python stats_hog_curve.py -i input.tsv -o output_dir -r 10
"""

import argparse
import pandas as pd
import random
from collections import defaultdict
import os


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='HOG Curve Statistics Data Calculation')
    parser.add_argument('-i', '--input_file', required=True,
                        help='Input TSV file with columns: SampleID HOG_ID Count (no header)')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for data files')
    parser.add_argument('-r', '--replicates', type=int, default=10,
                        help='Number of replicates for each sample size (default: 10)')
    return parser.parse_args()


def main():
    """Main function to perform HOG curve statistics data calculation."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Read input file (no header)
        print(f"Reading input file: {args.input_file}")
        df = pd.read_csv(args.input_file, sep='\t', header=None, names=['SampleID', 'HOG_ID', 'Count'])
        
        # Skip rows where Count is 0
        initial_rows = len(df)
        df = df[df['Count'] > 0]
        skipped_rows = initial_rows - len(df)
        print(f"Skipped {skipped_rows} rows with Count = 0")
        
        # Get unique samples and HOGs
        total_samples = df['SampleID'].nunique()
        total_hogs = df['HOG_ID'].nunique()
        
        print(f"Total samples: {total_samples}")
        print(f"Total HOGs: {total_hogs}")
        
        # Create interval bins and labels (same as stats_hog_count.py)
        intervals = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
        interval_labels = ['0%-20%', '20%-40%', '40%-60%', '60%-80%', '80%-100%']
        
        # Prepare data for curves
        print("Generating saturation curves data...")
        
        # Get all unique samples
        all_samples = df['SampleID'].unique().tolist()
        
        # Prepare data storage
        curve_data = []
        interval_curve_data = []
        
        # Generate data for different sample sizes
        sample_sizes = range(1, total_samples + 1)
        
        for n in sample_sizes:
            print(f"Processing sample size: {n}")
            
            # Repeat for replicates
            for rep in range(args.replicates):
                # Randomly select n samples
                selected_samples = random.sample(all_samples, n)
                
                # Get HOGs present in these samples
                selected_df = df[df['SampleID'].isin(selected_samples)]
                hogs_in_selected = selected_df['HOG_ID'].unique()
                
                # Calculate pan-genome size (total HOGs)
                pan_genome_size = len(hogs_in_selected)
                
                # Calculate core genome size (HOGs present in all selected samples)
                # For each HOG in selected samples, count how many selected samples it appears in
                hog_sample_count = selected_df.groupby('HOG_ID')['SampleID'].nunique().reset_index()
                core_hogs = hog_sample_count[hog_sample_count['SampleID'] == n]['HOG_ID']
                core_genome_size = len(core_hogs)
                
                # Add to curve data
                curve_data.append({
                    'Sample_Size': n,
                    'Replicate': rep,
                    'Pan_Genome_Size': pan_genome_size,
                    'Core_Genome_Size': core_genome_size
                })
                
                # Calculate HOG counts per interval for selected samples
                interval_counts = defaultdict(int)
                
                # Get presence count for each HOG in selected samples
                hog_sample_count = selected_df.groupby('HOG_ID')['SampleID'].nunique().reset_index()
                
                # Assign each HOG to interval based on presence count in selected samples
                for _, row in hog_sample_count.iterrows():
                    hog_in_samples = row['SampleID']
                    # Directly use presence count / sample size to determine interval
                    # No need to calculate actual percentage
                    presence_ratio = hog_in_samples / n
                    
                    if presence_ratio < 0.2:
                        interval_counts['0%-20%'] += 1
                    elif presence_ratio < 0.4:
                        interval_counts['20%-40%'] += 1
                    elif presence_ratio < 0.6:
                        interval_counts['40%-60%'] += 1
                    elif presence_ratio < 0.8:
                        interval_counts['60%-80%'] += 1
                    else:  # 0.8 <= presence_ratio <= 1.0
                        interval_counts['80%-100%'] += 1
                
                # Add to interval curve data
                interval_data = {
                    'Sample_Size': n,
                    'Replicate': rep
                }
                for label in interval_labels:
                    interval_data[label] = interval_counts.get(label, 0)
                interval_curve_data.append(interval_data)
        
        # Convert data to DataFrames
        curve_df = pd.DataFrame(curve_data)
        interval_curve_df = pd.DataFrame(interval_curve_data)
        
        # Save plotting data to TSV files
        print("Saving plotting data to TSV files...")
        curve_df.to_csv(os.path.join(args.output_dir, 'pan_core_curve_data.tsv'), sep='\t', index=False)
        interval_curve_df.to_csv(os.path.join(args.output_dir, 'interval_curve_data.tsv'), sep='\t', index=False)
        
        print("Data calculation completed successfully!")
        print(f"Results saved in: {args.output_dir}")
        print(f"  - pan_core_curve_data.tsv")
        print(f"  - interval_curve_data.tsv")
        print(f"\nTo generate plots, run: python plot_hog_curve.py -d {args.output_dir}")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == '__main__':
    main()
