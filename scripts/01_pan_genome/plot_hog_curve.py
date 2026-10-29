#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HOG Curve Plotting Script

This script reads pre-calculated curve data and generates plots with confidence intervals.

Usage:
    python plot_hog_curve.py -d data_dir -o output_prefix
    
Requirements:
    - pan_core_curve_data.tsv in data directory
    - interval_curve_data.tsv in data directory
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='HOG Curve Plotting Script')
    parser.add_argument('-d', '--data_dir', required=True,
                        help='Directory containing curve data files (pan_core_curve_data.tsv and interval_curve_data.tsv)')
    parser.add_argument('-o', '--output_prefix', default='./hog_curve',
                        help='Prefix for output files (will generate output_prefix.png and output_prefix.pdf) [default: ./hog_curve]')
    parser.add_argument('-w', '--fig_width', type=float, default=14.0,
                        help='Figure width in inches (default: 14.0)')
    parser.add_argument('-ht', '--fig_height', type=float, default=8.0,
                        help='Figure height in inches (default: 8.0)')
    parser.add_argument('--no_pan_core', action='store_true',
                        help='Skip pan-genome and core-genome curve plot')
    parser.add_argument('--no_interval', action='store_true',
                        help='Skip interval curve plot')
    return parser.parse_args()


def plot_pan_core_curve(curve_df, output_file, figsize):
    """Plot pan-genome and core-genome saturation curve with quartile shading."""
    print("Generating pan-genome saturation curve with quartile shading...")
    
    # ggpubr NPG color palette
    npg_colors = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F']
    
    # Calculate quartiles and median for each sample size
    curve_stats = curve_df.groupby('Sample_Size').agg({
        'Pan_Genome_Size': ['median', 'mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), 'count'],
        'Core_Genome_Size': ['median', 'mean', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), 'count']
    }).reset_index()
    
    # Flatten column names
    curve_stats.columns = ['Sample_Size', 
                          'Pan_Median', 'Pan_Mean', 'Pan_Q1', 'Pan_Q3', 'Pan_Count',
                          'Core_Median', 'Core_Mean', 'Core_Q1', 'Core_Q3', 'Core_Count']
    
    # Create plot
    plt.figure(figsize=figsize)
    
    # Plot pan-genome curve with quartile shading
    plt.plot(curve_stats['Sample_Size'], curve_stats['Pan_Median'], 
            color=npg_colors[0], linewidth=1, label='Pan-Genome')
    plt.fill_between(curve_stats['Sample_Size'], 
                    curve_stats['Pan_Q1'], 
                    curve_stats['Pan_Q3'], 
                    alpha=0.25, color=npg_colors[0])
    
    # Plot core-genome curve with quartile shading
    plt.plot(curve_stats['Sample_Size'], curve_stats['Core_Median'], 
            color=npg_colors[1], linewidth=1, label='Core-Genome')
    plt.fill_between(curve_stats['Sample_Size'], 
                    curve_stats['Core_Q1'], 
                    curve_stats['Core_Q3'], 
                    alpha=0.25, color=npg_colors[1])
    
    plt.title('Pan-Genome Saturation Curve', fontsize=10)
    plt.xlabel('Number of Samples', fontsize=8)
    plt.ylabel('HOG Count', fontsize=8)
    plt.legend(loc='best', fontsize=7)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save as PNG
    plt.savefig(f'{output_file}_pan_core.png', dpi=300)
    
    # Save as PDF
    plt.savefig(f'{output_file}_pan_core.pdf')
    
    plt.close()
    
    # Save curve statistics to TSV
    curve_stats.to_csv(f'{output_file}_pan_core_stats.tsv', sep='\t', index=False)
    
    print(f"  Saved: {output_file}_pan_core.png and {output_file}_pan_core.pdf")


def plot_interval_curves(interval_curve_df, output_prefix, figsize):
    """Plot HOG count change by presence rate interval with quartile shading."""
    print("Generating interval curves with quartile shading...")
    
    # ggpubr NPG color palette
    npg_colors = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F']
    
    # Define colors for each interval
    interval_colors = {
        '0%-20%': npg_colors[0],      # Red
        '20%-40%': npg_colors[1],    # Blue
        '40%-60%': npg_colors[2],     # Green
        '60%-80%': npg_colors[3],      # Purple
        '80%-100%': npg_colors[4]     # Orange
    }
    
    # Define interval labels in order
    interval_labels = ['0%-20%', '20%-40%', '40%-60%', '60%-80%', '80%-100%']
    
    plt.figure(figsize=figsize)
    
    # Calculate statistics for each interval
    interval_stats_list = []
    
    for interval in interval_labels:
        if interval not in interval_curve_df.columns:
            continue
            
        # Calculate quartiles, median and mean for each sample size
        interval_stats = interval_curve_df.groupby('Sample_Size')[interval].agg([
            'mean',
            lambda x: x.quantile(0.25),
            lambda x: x.quantile(0.75),
            'median',
            'count'
        ]).reset_index()
        
        interval_stats.columns = ['Sample_Size', 'Mean', 'Q1', 'Q3', 'Median', 'Count']
        interval_stats['Interval'] = interval
        
        interval_stats_list.append(interval_stats)
        
        # Plot the curve with quartile shading
        color = interval_colors.get(interval, '#808080')
        plt.plot(interval_stats['Sample_Size'], interval_stats['Median'], 
                linewidth=1, label=interval, color=color)
        plt.fill_between(interval_stats['Sample_Size'], 
                       interval_stats['Q1'], 
                       interval_stats['Q3'], 
                       alpha=0.25, color=color)
    
    plt.title('HOG Count Change by Presence Rate Interval', fontsize=10)
    plt.xlabel('Number of Samples', fontsize=8)
    plt.ylabel('HOG Count', fontsize=8)
    plt.legend(loc='best', fontsize=6)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save as PNG
    plt.savefig(f'{output_prefix}_interval.png', dpi=300)
    
    # Save as PDF
    plt.savefig(f'{output_prefix}_interval.pdf')
    
    plt.close()
    
    # Save interval statistics to TSV
    if interval_stats_list:
        interval_stats_df = pd.concat(interval_stats_list, ignore_index=True)
        interval_stats_df.to_csv(f'{output_prefix}_interval_stats.tsv', sep='\t', index=False)
    
    print(f"  Saved: {output_prefix}_interval.png and {output_prefix}_interval.pdf")


def main():
    """Main function to generate plots."""
    args = parse_arguments()
    
    try:
        # Check if data directory exists
        if not os.path.isdir(args.data_dir):
            print(f"Error: Data directory not found: {args.data_dir}")
            exit(1)
        
        # Check for required data files
        pan_core_file = os.path.join(args.data_dir, 'pan_core_curve_data.tsv')
        interval_file = os.path.join(args.data_dir, 'interval_curve_data.tsv')
        
        # Generate pan-core curve
        if not args.no_pan_core:
            if os.path.exists(pan_core_file):
                print(f"Reading data from: {pan_core_file}")
                curve_df = pd.read_csv(pan_core_file, sep='\t')
                plot_pan_core_curve(curve_df, args.output_prefix, (args.fig_width, args.fig_height))
            else:
                print(f"Warning: {pan_core_file} not found, skipping pan-core curve")
        
        # Generate interval curves
        if not args.no_interval:
            if os.path.exists(interval_file):
                print(f"Reading data from: {interval_file}")
                interval_df = pd.read_csv(interval_file, sep='\t')
                plot_interval_curves(interval_df, args.output_prefix, (args.fig_width, args.fig_height))
            else:
                print(f"Warning: {interval_file} not found, skipping interval curve")
        
        print("\nPlotting completed successfully!")
        print(f"Results saved with prefix: {args.output_prefix}")
        
    except Exception as e:
        print(f"Error during plotting: {str(e)}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == '__main__':
    main()
