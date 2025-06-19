#!/usr/bin/env python3
import sys

def compute_genetic_map(input_file, output_file, rate_cM_Mb=1.0):
    """
    Computes a simple genetic map assuming a uniform recombination rate.
    
    Parameters:
    - input_file: Path to the input file (must contain 'chrom pos id' columns)
    - output_file: Path where the output file will be written
    - rate_cM_Mb: Recombination rate in cM per Mb (default = 1.0)
    """
    try:
        with open(input_file) as f_in, open(output_file, "w") as f_out, open("ano_eagle_format.txt", "w") as f_eagle:
            # Write headers
            f_out.write("pos id distance(cM) cumulative(cM) rate(cM/Mb)\n")
            f_eagle.write("chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")

            prev_pos = None
            cumulative_cm = 0.0

            for line_num, line in enumerate(f_in, start=1):
                parts = line.strip().split()
                if len(parts) != 3:
                    print(f"Warning: Line {line_num} does not have exactly 3 columns, skipping.")
                    continue

                chrom, pos_str, snp_id = parts
                try:
                    pos = int(pos_str)
                except ValueError:
                    print(f"Warning: Invalid position on line {line_num}, skipping.")
                    continue

                distance_cm = 0.0
                if prev_pos is not None:
                    distance_bp = pos - prev_pos
                    distance_cm = (distance_bp / 1_000_000) * rate_cM_Mb
                    cumulative_cm += distance_cm
                    combined_rate = (distance_cm / distance_bp) * 1_000_000
                else:
                    combined_rate = 0.0  # First SNP has no previous position
                    cumulative_cm = 0.0

                # Write outputs
                f_out.write(f"{pos} {snp_id} {distance_cm:.6f} {cumulative_cm:.6f} {rate_cM_Mb:.2f}\n")
                f_eagle.write(f"{chrom} {pos} {combined_rate:.2f} {cumulative_cm:.6f}\n")

                prev_pos = pos

    except FileNotFoundError:   
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python create_map.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    compute_genetic_map(input_file, output_file)
