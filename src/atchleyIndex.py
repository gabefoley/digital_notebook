import pandas as pd

def getAtchleyIndexForSeq(seq, index):
    df = pd.DataFrame({'pah':[-0.60677211, 1.25436934, 1.20123120, 1.29519207, -0.95950963, -0.52784144, 0.32715696, -1.27633548, 1.82248357, -0.97854190, -0.72854024, 0.89146155, 0.21811094, 0.90061235, 1.52748331, 0.03475833, -0.01042641, -1.32574165, -0.61057537, 0.06016331],
                           'pss':[-1.08940440, 0.46126134, 0.24369180, -1.46537642, -0.43941703, 1.75873553, -0.52087575, -0.59487265, -0.53137304, -1.12494822, -1.41768155, 1.11359534, 1.96849543, -0.53619651, -0.01654753, 0.97046643, 0.65604029, -0.37200927, 0.01937737, 0.91703885],
                           'ms':[-0.92798080, -0.35464046, -1.47196724, -0.72868426, 1.15581211, -2.14269672, -0.42326104, 1.03001883, -0.02567446, 0.87095396, 0.30142692, -0.85232822, -0.15443545, 0.02364059, 0.65573629, -0.84138070, 0.01085433, 1.01019249, 1.50913663, 1.35527720],
                           'cc':[1.4726202, -1.1040471, -0.3243624, 0.1482127, -0.5634368, 0.7848519, -1.4308910, 0.4256597, -0.2432211, 1.4290018, -1.0227830, -0.4088754, 0.3903353, -0.4459318, 0.2563726, 0.9825769, 0.8716081, 1.4746391, -1.8958892, -0.7964405],
                           'ec':[0.251417050, -0.115742254, 2.127034065, 2.078732004, -0.175433772, -0.002451762, -0.494428309, 0.014711047, -1.571322338, -0.183751797, -0.094615066, 0.246465610, 0.472209055, 0.234557443, -2.460536571, 0.108360575, 0.093194098, 0.121129607, -0.441521905, -0.208006781]},
                           index = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q','R','S','T','V','W','Y'])

    seqscore = 0
    for idx in seq:
        seqscore += df[index][idx]

    return seqscore
