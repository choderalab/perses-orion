import os



ligand_pairs = [(1, 9), (1, 7), (0, 2), (4, 0), (0, 10), (1, 10), (0, 3), (3, 7), (6, 1), (4, 8), (1, 5), (0, 5), (6, 0), (1, 3), (3, 2), (9, 8)]

directories = [x for x in os.listdir() if x[0:3] == 'lig']


not_done = []
for a in ligand_pairs:
        foldername = f'lig{a[0]}to{a[1]}'
        if foldername not in directories:
            not_done.append(a)
print(not_done)
print(len(not_done))
