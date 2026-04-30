import argparse
import subprocess
import numpy as np
import os

env = os.environ.copy()


def run_singer_for(start, end, vcf_header, arg_save, ts_save=None,Ne=2e4,m=1e-8,mut_map=None,recomb_map=None, n=50, thin=50, polar=0.99,seed=None):
    if ts_save is None:
        ts_save = arg_save

    command = [
        "path/to/singer",
        "-Ne", str(Ne),
        "-vcf", vcf_header,
        "-output", str(arg_save),
        "-start", str(int(start)),
        "-end", str(int(end)),
        "-n", str(n),
        "-thin", str(thin),
        "-polar", str(polar),
    ]

    if mut_map is not None:
        command.extend(["-mut_map", mut_map])
    else:
        command.extend(["-m", str(m)])

    if recomb_map is not None:
        command.extend(["-recomb_map", recomb_map])
    
    if seed is None:
        pass
    else:
        command.append("-seed")
        command.append(str(seed))
    print("run start")
    
    
    subprocess.run(command, check=True, env=env, text=True)

    # Second command
    command2 = [
        'path/to/singer/convert_to_tskit', 
        '-input', str(arg_save),
        '-output', str(ts_save),
        '-start', '0',
        '-end', str(n),
        '-step', '1'
    ]

    subprocess.run(command2, check=True, env=env, text=True)

# def split_bins(start,end,length):
#     l=[]
#     n=np.ceil((end-start)/length)
#     if n==1:
#         return [(start,end)]
#     else:
#         for i in range(n-1):
#             l.append((start+i*length,start+(i+1)*length))
#         l.append((start+(n-1)*length,end))
#         return l

def main():
    # Initialize the parser
    parser = argparse.ArgumentParser(description="Run singer and convert to tskit")

    # Define the command-line arguments
    parser.add_argument("-Ne",type=float,default=2e4,help="Ne")
    parser.add_argument("-m",type=float,default=1e-8,help="mutation rate")
    parser.add_argument("-mut_map",type=str,default=None,help="mutation rate map path")
    parser.add_argument("-recomb_map",type=str,default=None,help="recombiantion rate map path")
    parser.add_argument("-start", type=float, required=True, help="Start index for the simulation")
    parser.add_argument("-end", type=float, required=True, help="End index for the simulation")
    parser.add_argument("-vcf_header", type=str, required=True, help="Path to the VCF header file")
    parser.add_argument("-arg_save", type=str, required=True, help="Path to save the VCF file")
    parser.add_argument("-ts_save", type=str, help="Path to save the TS file", default=None)
    parser.add_argument("-n", type=int, default=50, help="The number of iterations (default 50)")
    parser.add_argument("-thin", type=int, default=50, help="Thinning parameter (default 100)")
    parser.add_argument("-polar", type=float, default=0.99, help="Polar value (default 0.99)")
    parser.add_argument("-length",type=float,default=1e6,help="the segment length for running singer")
    parser.add_argument("-seed",type=int,default=None,help="seed")
    # Parse the arguments
    args = parser.parse_args()
    # l=split_bins(args.start,args.end,args.length)
    # # Call the function with parsed arguments
    # print("len of l",len(l))
    # for p in l:
        # scientific_s = "{:.1e}".format(p[0])
    run_singer_for(args.start, args.end, args.vcf_header, args.arg_save, args.ts_save,args.Ne,args.m,args.mut_map, args.recomb_map, args.n, args.thin, args.polar,args.seed)

if __name__ == "__main__":
    main()
