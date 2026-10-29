#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HOG Count Statistics Analysis

This script reads a TSV file containing SampleID, HOG_ID, and Count columns,
and performs the following statistical analysis:
1. Count the total number of samples and HOGs
2. Calculate the presence rate of each HOG across samples
3. Group HOGs by their presence rate in 10% intervals
4. For each sample, count HOGs in each presence rate interval
5. Output statistics results and sample-interval matrix to TSV files
"""

import argparse
import pandas as pd
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='HOG Count Statistics Analysis')
    parser.add_argument('-i', '--input_file', required=True,
                        help='Input TSV file with columns: SampleID HOG_ID Count (no header)')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Output TSV file for statistics results')
    parser.add_argument('-m', '--matrix_output_file', required=True,
                        help='Output TSV file for sample-HOG interval matrix')
    parser.add_argument('-l', '--hog_interval_file', required=True,
                        help='Output TSV file for HOG presence rate interval mapping')
    parser.add_argument('-p', '--hog_presence_file', required=True,
                        help='Output TSV file for each HOG specific presence rate')
    return parser.parse_args()


def main():
    """Main function to perform HOG count statistics."""
    args = parse_arguments()
    
    try:
        # Read input file (no header)
        print(f"Reading input file: {args.input_file}")
        df = pd.read_csv(args.input_file, sep='\t', header=None, names=['SampleID', 'HOG_ID', 'Count'])
        
        # Skip rows where Count is 0
        initial_rows = len(df)
        df = df[df['Count'] > 0]
        skipped_rows = initial_rows - len(df)
        print(f"Skipped {skipped_rows} rows with Count = 0")
        
        # 1. Count total number of samples and HOGs
        total_samples = df['SampleID'].nunique()
        total_hogs = df['HOG_ID'].nunique()
        
        print(f"Total samples: {total_samples}")
        print(f"Total HOGs: {total_hogs}")
        
        # 2. Calculate presence of each HOG across samples (Count > 0)
        print("Calculating HOG presence across samples...")
        # 注意：Count>0的筛选已在读取文件后立即执行（行号39-44）
        # 这里直接使用已筛选的数据计算HOG在样本中的存在情况
        hog_presence = df.groupby('HOG_ID')['SampleID'].nunique().reset_index()
        hog_presence.columns = ['HOG_ID', 'Present_Samples']
        
        # 3. Calculate presence rate and group into 10% intervals
        print("Grouping HOGs by presence rate intervals...")
        # Calculate presence rate
        hog_presence['Presence_Rate'] = (hog_presence['Present_Samples'] / total_samples) * 100
        
        # Create interval bins and labels (20% intervals)
        intervals = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
        interval_labels = ['0%-20%', '20%-40%', '40%-60%', '60%-80%', '80%-100%']
        
        # Count HOGs in each interval
        interval_counts = defaultdict(int)
        
        # Create a dictionary to map HOG_ID to its presence interval
        hog_to_interval = {}
        
        for _, row in hog_presence.iterrows():
            hog_id = row['HOG_ID']
            rate = row['Presence_Rate']
            interval_assigned = False
            for i, (lower, upper) in enumerate(intervals):
                if lower <= rate < upper:
                    interval_counts[interval_labels[i]] += 1
                    hog_to_interval[hog_id] = interval_labels[i]
                    interval_assigned = True
                    break
            # Handle the case where rate is exactly 100%
            if not interval_assigned or rate == 100:
                interval_counts['80%-100%'] += 1
                hog_to_interval[hog_id] = '80%-100%'
        
        # 4. For each sample, count HOGs in each presence rate interval
        print("Calculating HOG interval distribution per sample...")
        # Add interval information to the original dataframe
        df['Presence_Interval'] = df['HOG_ID'].map(hog_to_interval)
        
        # Create a pivot table to count HOGs per sample per interval
        sample_interval_counts = df.pivot_table(
            index='SampleID', 
            columns='Presence_Interval', 
            values='HOG_ID', 
            aggfunc='count',
            fill_value=0
        )
        
        # Ensure all interval labels are present as columns
        for label in interval_labels:
            if label not in sample_interval_counts.columns:
                sample_interval_counts[label] = 0
        
        # Sort columns in the correct order
        sample_interval_counts = sample_interval_counts[interval_labels]
        
        # 5. Create statistics result
        print("Generating statistics results...")
        # Create result dataframe
        stats_data = []
        
        # Add total counts
        stats_data.append(['Total_Samples', total_samples])
        stats_data.append(['Total_HOGs', total_hogs])
        stats_data.append(['', ''])
        # Add interval counts in reverse order (highest to lowest)
        for label in reversed(interval_labels):
            stats_data.append([label, interval_counts.get(label, 0)])
        
        # Create dataframe and save to TSV
        stats_df = pd.DataFrame(stats_data, columns=['Category', 'Count'])
        stats_df.to_csv(args.output_file, sep='\t', index=False)
        
        # 6. Save sample-interval matrix
        print(f"Saving sample-HOG interval matrix to: {args.matrix_output_file}")
        sample_interval_counts.to_csv(args.matrix_output_file, sep='\t')
        
        # 7. Save HOG presence rate interval mapping
        print(f"Creating HOG presence rate interval mapping file...")
        # Create a DataFrame with HOG_ID, Presence_Rate, and Presence_Rate_Interval
        hog_interval_df = pd.DataFrame({
            'HOG': hog_presence['HOG_ID'],
            'Presence_Rate': hog_presence['Presence_Rate'],
            'Presence_Rate_Interval': hog_presence['HOG_ID'].map(hog_to_interval)
        })
        # Save to TSV file
        hog_interval_df.to_csv(args.hog_interval_file, sep='\t', index=False)
        print(f"HOG presence rate interval mapping saved to: {args.hog_interval_file}")
        
        # 8. Save each HOG specific presence rate
        print(f"Creating each HOG specific presence rate file...")
        # Create a DataFrame with HOG_ID and Presence_Rate
        hog_presence_df = pd.DataFrame({
            'HOG_ID': hog_presence['HOG_ID'],
            'Presence_Rate': hog_presence['Presence_Rate']
        })
        # Save to TSV file
        hog_presence_df.to_csv(args.hog_presence_file, sep='\t', index=False)
        print(f"HOG specific presence rate saved to: {args.hog_presence_file}")
        
        print(f"Statistics results saved to: {args.output_file}")
        print("Analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == '__main__':
    main()