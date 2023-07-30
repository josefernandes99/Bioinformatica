class MyBlast:

    def __init__(self, filename = None, w = 3):
        if filename is not None:
            self.readDatabase(filename)
        else:
            self.db = []
        self.w = w
        self.map = None

    def readDatabase(self, filename):
        self.db = []
        with open(filename) as f:
            for line in f:
                self.db.append(line)
    
    def addSequenceDB(self, seq):
        self.db.append(seq)

    def buildMap (self, query):
        self.map = {}
        for i in range(0, len(query)-self.w+1):
            kmer = query[i:i+self.w]
            if kmer in self.map:
                self.map[kmer].append(i)
            else:
                self.map[kmer] = [i]

    def getHits (self, seq, query):
        hits = []
        for i in range(0, len(seq)-self.w+1):
            kmer = seq[i:i+self.w]
            if kmer in self.map:
                for j in self.map[kmer]:
                    hits.append((j,i))
        return hits

    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]:
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk

        k = 0
        matbw = 0
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]:
                matbw+=1
                bestk = k+1
            k+=1
        size += bestk

        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)

    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best


    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res

def test1():
    mb = MyBlast("seqBlast.txt", 11)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.bestAlignment(query)
    print(r)


def test2():
    mb = MyBlast("seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.bestAlignment(query2)
    print(r)
    
test1()
