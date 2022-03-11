import sys
import io


class FastqRecord:
    name: str
    sequence: str
    qualities: str

    def __init__(self, name, sequence, qualities):
        if len(sequence) != len(qualities):
            raise ValueError("Sequence and qualities must be of the same length.")
        self.name = name
        self.sequence = sequence
        self.qualities = qualities

    def to_record_string(self):
        return f"@{self.name}\n{self.sequence}\n+\n{self.qualities}\n"


def fastq_iter(filename: str):
    with open(filename, "rt", encoding="ascii") as file:
        # This saves one attribute lookup per line. Since FASTQ consists of
        # four lines per record, this has a noticeable impact on performance.
        readline = file.readline
        while True:
            name = readline()
            if not name:
                return
            if not name.startswith("@"):
                raise ValueError("FASTQ record should start with @.")
            sequence = readline()
            second_header = readline()
            qualities = readline()
            if not (sequence and second_header and qualities):
                raise EOFError(f"Truncated FASTQ file at {name}")
            if not second_header.startswith("+"):
                raise ValueError("Second header should start with +.")
            name = name[1:-1]  # strip @ and \n
            # Rstrip is faster than slicing for removing the final newline.
            sequence = sequence.rstrip()
            qualities = qualities.rstrip()
            # Length check in FastqRecord init.
            yield FastqRecord(name, sequence, qualities)


def read():
    for record in fastq_iter(sys.argv[1]):
        pass


def read_and_write():
    with open(sys.argv[2], "wt", encoding="ascii") as writer:
        # Given the amount of fastq records this has a notable impact.
        write = writer.write
        for record in fastq_iter(sys.argv[1]):
            write(record.to_record_string())


if __name__ == "__main__":
    read_and_write()

