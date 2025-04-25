import random
import sys

def random_sample_large_file(input_file, num_lines, output_file):
    random.seed(42)

    reservoir = []
    
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            if i < num_lines:
                reservoir.append(line)
            else:
                j = random.randint(0, i)
                if j < num_lines:
                    reservoir[j] = line

    with open(output_file, 'w') as out_f:
        out_f.writelines(reservoir)

    print(f"picked {num_lines} lines to {output_file}")


def main():
    if len(sys.argv) != 4:
        print("python random_pick_lines_no_repeat.py file 100000 output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    num_lines = int(sys.argv[2])
    output_file = sys.argv[3]

    random_sample_large_file(input_file, num_lines, output_file)


if __name__ == "__main__":
    main()