

# DEFINING A CLASS FOR AMINO ACIDS, p_alpha & p_beta DENOTE THE PROPENSITY VALUES
# defining a class for amino acids
# here p_alpha and p_beta denotes their propensity values

class AminoAcid:
    def __init__(self, p_alpha, p_beta, abbrv):
        # self.name = name
        self.abbrv = abbrv
        self.p_alpha = p_alpha
        self.p_beta = p_beta


# intiallised different amino acids as  object of the amino acid class

am1 = AminoAcid(1.45, 0.97, 'A')
am2 = AminoAcid(0.77, 1.30, 'C')
am3 = AminoAcid(0.98, 0.80, 'D')
am4 = AminoAcid(1.53, 0.26, 'E')
am5 = AminoAcid(1.12, 1.28, 'F')
am6 = AminoAcid(0.53, 0.81, 'G')
am7 = AminoAcid(1.24, 0.71, 'H')
am8 = AminoAcid(1, 1.6, 'I')
am9 = AminoAcid(1.07, 0.74, 'K')
am10 = AminoAcid(1.34, 1.22, 'L')
am11 = AminoAcid(1.2, 1.67, 'M')
am12 = AminoAcid(0.73, 0.65, 'N')
am13 = AminoAcid(0.59, 0.62, 'P')
am14 = AminoAcid(1.17, 1.23, 'Q')
am15 = AminoAcid(0.79, 0.9, 'R')
am16 = AminoAcid(0.79, 0.72, 'S')
am17 = AminoAcid(0.82, 1.2, 'T')
am18 = AminoAcid(1.14, 1.65, 'V')
am19 = AminoAcid(1.14, 1.19, 'W')
am20 = AminoAcid(0.61, 1.29, 'Y')


def get_alpha_score(seq):  # FUNCTION THAT RETURNS ALPHA SCORE OF A SEQUENCE OF AMINO ACIDS
    alpha_score = 0
    for i in range(len(seq)):
        if seq[i] == 'A':
            alpha_score += am1.p_alpha
        elif seq[i] == 'C':
            alpha_score += am2.p_alpha
        elif seq[i] == 'D':
            alpha_score += am3.p_alpha
        elif seq[i] == 'E':
            alpha_score += am4.p_alpha
        elif seq[i] == 'F':
            alpha_score += am5.p_alpha
        elif seq[i] == 'G':
            alpha_score += am6.p_alpha
        elif seq[i] == 'H':
            alpha_score += am7.p_alpha
        elif seq[i] == 'I':
            alpha_score += am8.p_alpha
        elif seq[i] == 'K':
            alpha_score += am9.p_alpha
        elif seq[i] == 'L':
            alpha_score += am10.p_alpha
        elif seq[i] == 'M':
            alpha_score += am11.p_alpha
        elif seq[i] == 'N':
            alpha_score += am12.p_alpha
        elif seq[i] == 'P':
            alpha_score += am13.p_alpha
        elif seq[i] == 'Q':
            alpha_score += am14.p_alpha
        elif seq[i] == 'R':
            alpha_score += am15.p_alpha
        elif seq[i] == 'S':
            alpha_score += am16.p_alpha
        elif seq[i] == 'T':
            alpha_score += am17.p_alpha
        elif seq[i] == 'V':
            alpha_score += am18.p_alpha
        elif seq[i] == 'W':
            alpha_score += am19.p_alpha
        else:
            alpha_score += am20.p_alpha
    return alpha_score


def get_beta_score(seq):  # FUNCTION THAT RETURNS BETA SCORE OF A SEQUENCE OF AMINO ACIDS
    beta_score = 0
    for i in range(len(seq)):
        if seq[i] == 'A':
            beta_score += am1.p_beta
        elif seq[i] == 'C':
            beta_score += am2.p_beta
        elif seq[i] == 'D':
            beta_score += am3.p_beta
        elif seq[i] == 'E':
            beta_score += am4.p_beta
        elif seq[i] == 'F':
            beta_score += am5.p_beta
        elif seq[i] == 'G':
            beta_score += am6.p_beta
        elif seq[i] == 'H':
            beta_score += am7.p_beta
        elif seq[i] == 'I':
            beta_score += am8.p_beta
        elif seq[i] == 'K':
            beta_score += am9.p_beta
        elif seq[i] == 'L':
            beta_score += am10.p_beta
        elif seq[i] == 'M':
            beta_score += am11.p_beta
        elif seq[i] == 'N':
            beta_score += am12.p_beta
        elif seq[i] == 'P':
            beta_score += am13.p_beta
        elif seq[i] == 'Q':
            beta_score += am14.p_beta
        elif seq[i] == 'R':
            beta_score += am15.p_beta
        elif seq[i] == 'S':
            beta_score += am16.p_beta
        elif seq[i] == 'T':
            beta_score += am17.p_beta
        elif seq[i] == 'V':
            beta_score += am18.p_beta
        elif seq[i] == 'W':
            beta_score += am19.p_beta
        else:
            beta_score += am20.p_beta
    return beta_score


def get_alpha_helix(seq):  # FUNCTION THAT RETURNS THE SITES WHERE ALPHA HELIX COULD BE FORMED
    s = []
    for i in range(len(seq)):
        s.append('_')
    if len(seq) < 6:
        print("Sequence length must be atleast 6 for this algorithm to work")
        return
    for i in range(len(seq)-5):
        n = 0
        for j in range(i, i+6):
            if get_alpha_score(seq[j]) >= 1:
                n += 1
        if n >= 4:
            for j in range(i, i + 6):
                if s[j] != 'H':
                    s[j] = 'H'
            p1 = i+6
            p2 = i-1
            while p1 < len(seq):
                if get_alpha_score(seq[p1-3:p1+1]) >= 4:
                    s[p1] = 'H'
                else:
                    break
                p1 += 1
            while p2 >= 0:
                if get_alpha_score(seq[p2:p2+4]) >= 4:
                    s[p2] = 'H'
                else:
                    break
                p2 -= 1
    alpha_helix = ""
    for i in range(len(s)):
        alpha_helix += s[i]
    return alpha_helix


def get_beta_sheet(seq):  # FUNCTION THAT RETURNS THE SITES WHERE BETA SHEETS COULD BE FORMED
    s = []
    for i in range(len(seq)):
        s.append('_')
    if len(seq) < 5:
        print("Sequence length must be atleast 5 for this algorithm to work")
        return
    for i in range(len(seq)-4):
        n = 0
        for j in range(i, i+5):
            if get_beta_score(seq[j]) >= 1:
                n += 1
        if n >= 3:
            for j in range(i, i+5):
                if s[j] != 'S':
                    s[j] = 'S'
            p1 = i+5
            p2 = i-1
            while p1 < len(seq):
                if get_beta_score(seq[p1-3:p1+1]) >= 4:
                    s[p1] = 'S'
                else:
                    break
                p1 += 1
            while p2 >= 0:
                if get_beta_score(seq[p2:p2+4]) >= 4:
                    s[p2] = 'S'
                else:
                    break
                p2 -= 1
    beta_sheet = ""
    for i in range(len(s)):
        beta_sheet += s[i]
    return beta_sheet






seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"  # INPUT SEQUENCE
seq1 = get_alpha_helix(seq)
 # SITES WHERE ALPHA HELIXES ARE FORMED
seq2 = get_beta_sheet(seq)
print("\n \n")
# SITES WHERE BETA SHEETS ARE FORMED
# final_seq = conf_reso(seq1, seq2, seq)  # FINAL SEQUENCE OBTAINED AFTER CONFLICT RESOLUTION
res = []
for i in range(len(seq)):
    res.append("")
    i = 0
while i < len(seq):
    if (seq1[i] == 'H' and seq2[i] == '_') or (seq1[i] == '_' and seq2[i] == 'H'):
        res[i] = 'H'
        i += 1
    elif (seq1[i] == 'S' and seq2[i] == '_') or (seq1[i] == '_' and seq2[i] == 'S'):
        res[i] = 'S'
        i += 1
    elif seq1[i] == '_' and seq2[i] == '_':
        res[i] = '-'
        i += 1
    else:
        n = 0
        while (seq1[i] == 'S' and seq2[i] == 'H') or (seq1[i] == 'H' and seq2[i] == 'S'):
                n += 1
                if i < len(seq)-1:
                    i += 1
                else:
                    break
        if i == len(seq)-1:
                i += 1
        p1 = get_alpha_score(seq[i-n:i])
        p2 = get_beta_score(seq[i-n:i])
        if p1 > p2:
            for k in range(i-n,i):
                    res[k] = 'H'
        else:
            for k in range(i-n,i):
                    res[k] = 'S'
ans = ""
for i in range(len(res)):
        ans += res[i]

print(seq)
print("\n \n")
print("get alpha helix as follows\n")
print(seq1)
print("\n \n")
print("get beta sheet is as follows")
print(seq2)
print("\n \n")

print("final struture as follows ")

print(ans)