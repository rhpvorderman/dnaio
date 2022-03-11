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

    def fastq_bytes(self):
        return f"@{self.name}\n{self.sequence}\n+\n{self.qualities}\n".encode("ascii")


def fastq_iter(filename: str):
    with open(filename, "rt", encoding="ascii") as file:
        while True:
            name = file.readline()
            if not name:
                return
            if not name.startswith("@"):
                raise ValueError("FASTQ record should start with @.")
            sequence = file.readline()
            second_header = file.readline()
            qualities = file.readline()
            if not (sequence and second_header and qualities):
                raise EOFError(f"Truncated FASTQ file at {name}")
            if not second_header.startswith("+"):
                raise ValueError("Second header should start with +.")
            name = name[1:-1]  # strip @ and \n
            sequence = sequence.rstrip()  # Rstrip is faster than slicing.
            qualities = qualities.rstrip()
            # Length check in FastqRecord init.
            yield FastqRecord(name, sequence, qualities)


if __name__ == "__main__":
    for record in fastq_iter(sys.argv[1]):
        pass
