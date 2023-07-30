import os

os.system('cls')


class NumMatrix:
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.matrix = [[0.0] * cols for _ in range(rows)]
        
    def __str__(self):
        return '\n'.join(['\t'.join(map(str, row)) for row in self.matrix])
        

class PWM:
    def __init__(self, lst_of_words=[], bio_type='DNA', alphabet='ACGT'):
        self.lst_of_words = lst_of_words
        self.bio_type = bio_type
        self.alphabet = alphabet
        self.pwm = NumMatrix(len(alphabet), len(lst_of_words[0])) if lst_of_words else NumMatrix(len(alphabet), 0)
        
        if lst_of_words:
            for i in range(len(lst_of_words[0])):
                col = [word[i] for word in lst_of_words]
                for j, symbol in enumerate(self.alphabet):
                    self.pwm.matrix[j][i] = col.count(symbol) / len(col)
    
    def log_odds(self, bckgrd_freq={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        for i in range(self.pwm.rows):
            bk = bckgrd_freq[self.alphabet[i]]
            for j in range(self.pwm.cols):
                Mkj = self.pwm.matrix[i][j]
                self.pwm.matrix[i][j] = (Mkj / bk) if Mkj > 0 else -100.0
                
    def test(self, target_seq):
        print(f'List of words: {self.lst_of_words}')
        print(f'Biotype: {self.bio_type}')
        print(f'Alphabet: {self.alphabet}')
        print('PWM before log odds calculation:')
        print(self.pwm)
        self.log_odds()
        print('PWM after log odds calculation:')
        print(self.pwm)
        print(f'Target sequence: {target_seq}')
        
        
# Test the implementation
print("\n################# TEST 1 #################\n")
pwm1 = PWM()
pwm1.test('CTCCAGAGTTCCTGCCCTGCCACTCCCTTGACTTTCTCTCTTTTTCTCCCTCCCCTCCGCCGCCGCCCACACCCGAGTGCTAGCTGTGGACCTAGGGGCTCTTGATCTACCACTGAGCCACACCCAAGTGTCCAGGCTTATCATCTTGCCTTAGCTTTTTCTGTAGCACATTATCAATGCCTGGCCAGTGTCTTGTGTATTTATCTTCTTTTTCCTTTCTCTCCACCTGCTTGGCTGTCTCATGATAAC')

print("\n################# TEST 2 #################\n")
pwm2 = PWM(['AACGTG', 'AAGTTC', 'ACGTGC', 'CCGTGG', 'TACGTT'])
pwm2.test('CTCCAGAGTTCCTGCCCTGCCACTCCCTTGACTTTCTCTCTTTTTCTCCCTCCCCTCCGCCGCCGCCCACACCCGAGTGCTAGCTGTGGACCTAGGGGCTCTTGATCTACCACTGAGCCACACCCAAGTGTCCAGGCTTATCATCTTGCCTTAGCTTTTTCTGTAGCACATTATCAATGCCTGGCCAGTGTCTTGTGTATTTATCTTCTTTTTCCTTTCTCTCCACCTGCTTGGCTGTCTCATGATAAC')

print()