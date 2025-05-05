import pandas as pd
import numpy as np

def get_n_e_b_list_from_file(fn):
    '''
    Get num of uniq kmers, eps list, segment count list
    from file
    '''
    with open(fn) as f:
        lines = f.readlines()
    
    x_list = []
    y_list = []
    is_first = True
    for line in lines:
        if is_first:
            is_first = False
            uniq_kmers = int(line)
            continue
        info = line.split()
        x_list.append(int(info[0]))
        y_list.append(int(info[1]))
    return uniq_kmers, x_list, y_list

def get_ribbon(df, n, alpha):
    ratio = n / (df.segments * df.epsilon**alpha)
    l, h = ratio.min(), ratio.max()
    return l, h, h - l

def get_capla(file_path, epsilon_filter=None):
    from scipy.optimize import golden

    with open(file_path, "r") as f:
        n = int(f.readline().strip())
        df = pd.read_csv(f, sep=" ", header=None, names=["epsilon", "segments"])

    if epsilon_filter is not None:
        df = df[df.epsilon.isin(epsilon_filter)]

    if 1 not in df.epsilon.values:
        raise ValueError("Epsilon 1 not found in the file")

    df = df[df.segments > 1]
    df_no_eps1 = df[df.epsilon != 1]
    eps1_segments = df[df.epsilon == 1].segments.values[0]
    log_eps_ratio = np.log(eps1_segments / df_no_eps1.segments) / np.log(df_no_eps1.epsilon)
    alpha_l, alpha_h = np.min(log_eps_ratio), np.max(log_eps_ratio)
    opt_alpha = golden(lambda a: get_ribbon(df, n, a)[2], brack=(alpha_l, alpha_h))
    l, h, w = get_ribbon(df, n, opt_alpha) # beta_min, beta_max, beta_diff

    return df, n, opt_alpha, w, l, h


def write_values(out_fn, data):
    import csv
    import os    
    
    header = ['Genome', 'Uniq_kmers', 'K', 'Alpha', 
                'Beta_min', 'Beta_max']
    file_exists = os.path.isfile(out_fn)
    
    with open(out_fn, mode='a') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(header)
        writer.writerow(data)
        

def main():
    import sys
    segment_file = sys.argv[1]
    genome = sys.argv[2]
    kmer = sys.argv[3]
    out_fn = sys.argv[4]
    
    df, n, opt_alpha, w, l, h = get_capla(segment_file)
    
    value_list = [genome, n, kmer, opt_alpha, l, h]
    write_values(out_fn, value_list)

main()