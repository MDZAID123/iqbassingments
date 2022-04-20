
# Sequence obtained from my implementation fo chou fastman algorithm

found1 ="_HHHHHHHHHHHSSSSSSSSSSSSSSHHHHHHHHHSSSSHHHHHHHHHHHH__HHHHHHHHHHSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS____HSSSSSSSSSSSSSSSSSSSSS__SSSSSSSSSS___HHHHHH_________"
# Seqeuence obtained from stride web server

found2 = "TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE "
for i, j in enumerate(found1):
  print(f"At index {i+1} -  {j} and {found2[i]}")
s=""
for i in range(len(found1)):
    if found2[i] != found1[i]:
        s += 'n'
    else:
        s += "_"
print(s)
# INDICES MARKED WITH n MEANS THAT PROGRAM OUTPUT AND STRIDE OUTPUT WERE DIFFERENT AT THESE POSITIONS