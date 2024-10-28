import os
import fypp
import argparse


def apply_command(folder, with_qp, with_xqp, with_hp):
    # Get a list of all files in the folder
    # Filter files based on the provided extension
    lfiles = os.listdir(folder)
    files = [folder+os.sep+file for file in lfiles if file.endswith(".fypp")]

    args = []
    if with_qp:
        args.append("-DWITH_QP=True")
    if with_xqp:
        args.append("-DWITH_XQP=True")
    if with_hp:
        args.append("-DWITH_HP=True")

    optparser = fypp.get_option_parser()
    options, leftover = optparser.parse_args(args=args)
    tool = fypp.Fypp(options)
    # Apply the command line to each file
    for file in files:
        source_file = file
        target_file = file.removesuffix('fypp')+'f90'
        tool.process_file(source_file, target_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Preprocess LightKrylov source files.')
    # parser.add_argument("--with_qp", action="store_true",
    #                     help="Include WITH_QP in the command")
    # parser.add_argument("--with_xqp", action="store_true",
    #                     help="Include WITH_XQP in the command")
    # parser.add_argument("--with_hp", action="store_true",
    #                     help="Include WITH_HP in the command")

    args = parser.parse_args()
    # apply_command(args.with_qp, args.with_xqp, args.with_hp)
    apply_command("./src/", False, False, False)
    apply_command("./test/", False, False, False)
