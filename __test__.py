exec(open('rpythonplus.py').read())


def test_2() : #r_liste cree avec des py_listes

    import numpy

    robjects.r.source('/home/pierre/_BIOINFORMATIQUE/M2_IPFB/conversion_data.r')

    tableau = base.data_frame(x = numpy.arange(1, 10, 1).tolist(), y = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])

    empilement = robjects.globalenv['empilement']

    destack_dataframe = robjects.globalenv['destack_dataframe']

    gamma = destack_dataframe(empilement(tableau), 1, True)

    sortie = r_vers_python_(gamma)

    print('')

    print('---------------test 2---------------')

    print('r_liste cree avec des py_listes, vers py_dict')

    print('')

    print(sortie)



def test_5() : #r_liste de test statistique

    stats = importr('stats')

    r_df = rpy2.robjects.r['iris']

    temp = stats.kruskal_test(r_df.rx2(1), r_df.rx2(5))

    sortie = r_vers_python_(temp)

    print('')

    print('---------------test 5---------------')

    print('r_liste de test statistique, vers py_dict')

    print('')

    print(sortie)



def test_1() :#r_liste vers py_dict nb r_vecteurs > nb_cpu

    vec1 = python_vers_r_(['0', '1', '2'])

    vec2 = python_vers_r_(['a', 'b', 'c'])

    vec3 = python_vers_r_(['c', 'd', 'e'])

    vec4 = python_vers_r_(['8', '10', '12e'])

    vec5 = python_vers_r_([0, 0, 10])

    liste = base.list(vec1, vec2, vec3, vec4, vec5)

    print('')

    print('---------------test 1---------------')

    print('r_liste (nb r_vecteurs superieur a nb_cpu), vers py_dict')

    print('')

    sortie = r_vers_python_(liste)

    print(sortie)



def test_10() : #py_dataframe une colonne vers r_dataframe

    test_df = pandas.read_csv('/home/pierre/_BIOINFORMATIQUE/_EXPERIMENTATIONS/conversion_r_vers_python/fichiers_test/test10.tsv', sep  ='\t', index_col = 0)

    print('')

    print('---------------test 10---------------')

    print('py_dataframe 1 colonne, vers r_dataframe')

    print('')

    sortie = python_vers_r_(test_df)

    print(sortie)



def test_11() : #py_dataframe autant de colonnes que de cpu vers r_dataframe

    test_df = pandas.read_csv('/home/pierre/_BIOINFORMATIQUE/_EXPERIMENTATIONS/conversion_r_vers_python/fichiers_test/test11.tsv', sep  ='\t')

    print('')

    print('---------------test 11---------------')

    print('py_dataframe avec autant de colonnes que de cpu, vers r_dataframe')

    print('')

    sortie = python_vers_r_(test_df)

    print(sortie)



def test_12() : #py_dataframe autant de colonnes que de cpu vers r_dataframe

    test_df = pandas.read_csv('/home/pierre/_BIOINFORMATIQUE/_EXPERIMENTATIONS/conversion_r_vers_python/fichiers_test/test12.tsv', sep  = '\t')

    print('')

    print('---------------test 12---------------')

    print('py_dataframe plus de colonnes que de cpu, vers r_dataframe')

    print('')

    sortie = python_vers_r_(test_df)

    print(sortie)



def test_13() :

    test_df = utils.read_delim(file = '/home/pierre/_BIOINFORMATIQUE/_EXPERIMENTATIONS/conversion_r_vers_python/fichiers_test/test13.tsv', sep  ='\t', header = True)

    print('')

    print('---------------test 13---------------')

    print('r_dataframe plus de colonnes que de cpu, vers py_dataframe')

    print('')

    sortie = r_vers_python_(test_df)

    print(sortie)

def test_14() :

    import numpy

    mat = numpy.random.randint(10, size=(3, 3))

    print('')

    print('---------------test 14---------------')

    print('matrice numpy vers matrice R')

    print('')

    sortie = python_vers_r_(mat)

    print(sortie)

def test_15() :

    print('')

    print('---------------test 15---------------')

    print('matrice R vers matrice numpy')

    print('')

    mat = base.matrix(python_vers_r_(list(numpy.arange(1, 10))), ncol = 3, nrow = 3)

    sortie = r_vers_python_(mat)

    print(sortie)

def test_16() :

    print('')

    print('---------------test 16---------------')

    print('matrice pandas dataframe 1 ligne vers R')

    print('')

    df1 = pandas.DataFrame({'Prenom' : ['Pierre'], 'Age' : [18], 'Ville' : ['Saint_Brieuc']})

    sortie = python_vers_r_(df1)

    print(sortie)



print('')

print('--------------------------------------------------')

print('---------------r_liste vers py_dictionnaire---------------')

print('--------------------------------------------------')

print('')

test_2()

test_5()

test_1()

print('')

print('--------------------------------------------------')

print('---------------py_dataframe vers r_dataframe---------------')

print('--------------------------------------------------')

print('')

test_10()

test_11()

test_12()

test_16()

print('')

print('-----------------------------------------------------------')

print('---------------r_dataframe vers py_dataframe---------------')

print('-----------------------------------------------------------')

print('')

test_13()

print('')

print('--------------------------------------------------')

print('---------------gestion des matrices---------------')

print('--------------------------------------------------')

print('')

test_14()
