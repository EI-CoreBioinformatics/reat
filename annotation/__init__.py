import sys

VERSION = '0.4.1'


def report_errors(errors, samples):
    if any([len(error_list) for error_list in errors.values()]):
        print(f"File {samples.name} parsing failed, errors found:\n", file=sys.stderr)
        for line, error_list in errors.items():
            if not error_list:
                continue
            print("Line:", line.strip(), sep='\n\t', file=sys.stderr)
            print("was not parsed successfully, the following errors were found:", file=sys.stderr)
            [print("\t-", e, file=sys.stderr) for e in error_list]
        raise ValueError(f"Could not parse file {samples.name}")
