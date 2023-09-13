sequences = {
    '1csp': "MLEGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTKEA"
}

def get_sequence(name : str) -> str:

    for key,seq in sequences.items():
        if name.find(key) >= 0:
            return seq

    raise ValueError("Unknown sequece for '%s'" % name)