import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("num1", help="Number 1 1")
    parser.add_argument("num2", help="Number 2")
    args = parser.parse_args()

    print("Number 1:", args.num1)
    print("Number 2:", args.num2)

    print("Sum: {}".format(int(args.num1) + int(args.num2)))
