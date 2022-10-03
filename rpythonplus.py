#RPythonPlus

#Author Pierre GUYOMARD
#Contact piguyomard@outlook.com

import multiprocessing

import math

import re

import pandas

import numpy

import rpy2

import gc

import rpy2.robjects as robjects

from rpy2.robjects.packages import importr

utils = importr('utils')

base = importr('base')

robjects.r.source('/home/pierre/miniconda3/envs/julia_env2/lib/python3.7/site-packages/' + 'rpythonplus_toolbox.r')

is_numeric = (int, float, numpy.float16, numpy.float32, numpy.float64, numpy.float128, numpy.int8, numpy.int16, numpy.int32, numpy.int64)


###########################################################################################################################
###########################################################################################################################
#							R -> PYTHON
###########################################################################################################################
###########################################################################################################################



###########################################################################################################################
#FONCTION PRINCIPALE
###########################################################################################################################

def r_vers_python_(objet) :

	vecteur_type_r_liste = (rpy2.robjects.vectors.IntVector, rpy2.robjects.vectors.FloatVector, rpy2.robjects.vectors.StrVector, rpy2.robjects.vectors.IntVector)

	if isinstance(objet, rpy2.robjects.vectors.DataFrame) == True :

		sortie = r_dataframe_vers_python_dataframe(objet)

	elif len(re.findall('Matrix', str(type(objet)))) > 0 :

		sortie = base.as_data_frame(objet)

		sortie = r_dataframe_vers_python_dataframe(sortie)

		sortie = sortie.to_numpy()

	elif isinstance(objet, vecteur_type_r_liste) == True :

		sortie = r_vecteur_vers_python_liste(objet)

	elif isinstance(objet, rpy2.rinterface_lib.sexp.NULLType) == True : #les r_vecteurs de longueur nul ont ce typage

		sortie = []

	elif isinstance(objet, rpy2.robjects.vectors.ListVector) == True and isinstance(objet.names, rpy2.rinterface_lib.sexp.NULLType) != True :

		sortie = r_liste_vers_python_dictionnaire_names(objet)

	elif isinstance(objet, rpy2.robjects.vectors.ListVector) == True and isinstance(objet.names, rpy2.rinterface_lib.sexp.NULLType) == True : #les listes vides sont gérées ici

		sortie = r_liste_vers_python_dictionnaire_nonames(objet)

	else :

		sortie = objet

		print('\n')

		print(type(objet), ': pas de conversion disponible')

		print('\n')



	return sortie

########################################################################################################################
#R-LISTE VERS PYTHON-DICTIONNAIRE
########################################################################################################################


def r_liste_vers_python_dictionnaire_nonames(objet) :

	if len(objet) == 0 : #si liste vide

		sortie = {}

	else :

		alpha = r_liste_vers_p_dictionnaire_noname(objet, 4)

		#multiprocessing.set_start_method('spawn', force=True)

		pool = multiprocessing.Pool(len(alpha.index))

		beta = pool.map(alpha.process_k_mers, [k for k in range(1, len(alpha.k_mers_elements) + 1)])

		pool.close()

		gc.collect()

		gamma = {}

		c_fin = 0

		while c_fin <= len(beta) - 1 :

			temp_keys = list(beta[c_fin].keys())

			c_fin2 = 0

			while c_fin2 <= len(temp_keys) - 1 :

				gamma[temp_keys[c_fin2]] = beta[c_fin][temp_keys[c_fin2]]

				gc.collect()

				c_fin2 = c_fin2 + 1

			gc.collect()

			c_fin = c_fin + 1

		gamma.pop('*!*!*!*')

	return gamma





class r_liste_vers_p_dictionnaire_noname :

	def __init__(self, liste, nb_cpu) :

		self.r_liste = liste

		self.names = list(range(1, len(self.r_liste) + 1))



		self.index = {}

		self.k_mers_elements = {}

		self.k_mers_names = {}





		k = 1

		minimum = 0

		maximum = math.ceil(len(self.r_liste) / nb_cpu)

		etendue = math.ceil(len(self.r_liste) / nb_cpu)

		while k <= nb_cpu :

			if maximum > len(self.r_liste) :

				maximum = len(self.r_liste)

			self.index[k] = list(range(minimum, maximum, 1))

			self.k_mers_elements[k] = self.r_liste[minimum:maximum]

			self.k_mers_names[k] = self.names[minimum:maximum]



			minimum = maximum

			maximum = minimum + etendue

			gc.collect()

			k = k + 1

		if nb_cpu > 1 :

			self.ajout = - (len(self.k_mers_elements[len(self.k_mers_elements)]) - len(self.k_mers_elements[len(self.k_mers_elements) - 1]))

			r_add_vecteur = robjects.globalenv['r_add_vecteur']

			self.k_mers_elements[len(self.k_mers_elements)] = r_add_vecteur(self.k_mers_elements[len(self.k_mers_elements)], self.ajout)

			self.k_mers_names[len(self.k_mers_names)] = self.k_mers_names[len(self.k_mers_names)] + ['*!*!*!*'] * self.ajout

		else :

			pass

	def process_k_mers(self, k) :

		sortie_dict = {}

		fehu = 0

		while fehu <= len(self.k_mers_elements[k]) - 1 :

			sortie_dict[self.k_mers_names[k][fehu]] = self.k_mers_elements[k][fehu]

			gc.collect()

			fehu = fehu + 1

		return sortie_dict











def r_liste_vers_python_dictionnaire_names(objet) : #convertir un vecteur r en liste python

	nb_cpu = multiprocessing.cpu_count()



	if nb_cpu >= len(objet) or math.floor((len(objet)/(2*nb_cpu))) == 1 :

		nb_cpu = math.ceil(min(len(objet), nb_cpu) / 10)

	pre_sortie = r_liste_vers_p_dictionnaire_names(objet, nb_cpu)



	if nb_cpu > 1 :

		#multiprocessing.set_start_method('spawn', force=True)

		pool = multiprocessing.Pool(len(pre_sortie.index))

		dictionnaires_k_mers = pool.map(pre_sortie.process_k_mers, [k for k in range(1, len(pre_sortie.k_mers_elements) + 1)]) #le premier k-mer a pour indice 1 et le dernier k-mer a pour indice la longueur de l'index + 1

		pool.close()

	else :

		return pre_sortie.process_monomere()



	compteur = 0 #cette boucle rassemble toutes les clefs de tous les k-mers

	pre_liste_clefs = []

	while compteur <= len(dictionnaires_k_mers) - 1:

		pre_liste_clefs.extend(list(dictionnaires_k_mers[compteur].keys()))

		gc.collect()

		compteur = compteur + 1



	omicron = 0 #cette boucle supprime les doublons de clefs dans pre_liste_clefs

	liste_clefs = []

	dictionnaire_merge_k_mers = {}

	while omicron <= len(pre_liste_clefs) - 1 :

		if pre_liste_clefs[omicron] not in liste_clefs :

			liste_clefs.append(pre_liste_clefs[omicron])

			dictionnaire_merge_k_mers[pre_liste_clefs[omicron]] = []

		gc.collect()

		omicron = omicron + 1




	#fusionne les dictionnaires issus des k-mers

	compteur_ouf = 0

	while compteur_ouf <= len(dictionnaires_k_mers) - 1 : #boucle sur toutes les clefs d'un k-mer

		omega = 0

		while omega <= len(liste_clefs) - 1 : #boucle sur les clefs du dictionnaire global (dictionnaire_merge_k_mers)

			prec = dictionnaire_merge_k_mers[liste_clefs[omega]] #capture les sortie deja presentes pour la clef

			temp_ouf = []

			if liste_clefs[omega] in dictionnaires_k_mers[compteur_ouf] :



				manant = dictionnaires_k_mers[compteur_ouf][liste_clefs[omega]]



				liste_typage = (str, float, int)

				if isinstance(manant, liste_typage) == True :

					temp_ouf.append(manant)

				else :

					temp_ouf.extend(manant)

				if len(temp_ouf) == 1 :

					dictionnaire_merge_k_mers[liste_clefs[omega]] = temp_ouf[0]

				elif prec == [] :

					dictionnaire_merge_k_mers[liste_clefs[omega]] = temp_ouf

				else :

					try :

						dictionnaire_merge_k_mers[liste_clefs[omega]] = prec + temp_ouf

					except TypeError :

						dictionnaire_merge_k_mers[liste_clefs[omega]] = list(prec) + temp_ouf

			else :

				pass

			gc.collect()

			omega = omega + 1

		gc.collect()

		compteur_ouf = compteur_ouf + 1

	try :

		dictionnaire_merge_k_mers.pop('*!*!*!*')

	except KeyError : #si aucune 'fausse' clef n'a pas ete ajoutee

		pass

	gc.collect()

	return dictionnaire_merge_k_mers








class r_liste_vers_p_dictionnaire_names :

	def __init__(self, liste, nb_cpu) :

		self.r_liste = liste

		self.index = {}

		self.k_mers_elements = {}

		self.k_mers_names = {}

		if nb_cpu > 1 :

			k = 1

			minimum = 0

			maximum = math.ceil(len(self.r_liste) / nb_cpu)

			etendue = math.ceil(len(self.r_liste) / nb_cpu)

			while k <= nb_cpu :

				if maximum > len(self.r_liste) :

					maximum = len(self.r_liste)

				self.index[k] = list(range(minimum, maximum, 1))

				self.k_mers_elements[k] = self.r_liste[minimum:maximum]

				self.k_mers_names[k] = self.r_liste.names[minimum:maximum]



				minimum = maximum

				maximum = minimum + etendue


				k = k + 1

			self.ajout = len(self.k_mers_elements[len(self.k_mers_elements) - 1]) - len(self.k_mers_elements[len(self.k_mers_elements)])

		if nb_cpu == 1 :

			k = 1

			self.index[k] = list(range(1, len(self.r_liste), 1))

			self.k_mers_elements[k] = liste

			self.k_mers_names[k] = self.r_liste.names

			self.ajout = 0





		r_add_vecteur = robjects.globalenv['r_add_vecteur']

		self.k_mers_elements[len(self.k_mers_elements)] = r_add_vecteur(self.k_mers_elements[len(self.k_mers_elements)], self.ajout)

		self.k_mers_names[len(self.k_mers_names)] = r_add_vecteur(self.k_mers_names[len(self.k_mers_names)], self.ajout)

		gc.collect()



	def process_k_mers(self, k) :

		liste_noms = []

		intervalles = []

		liste_elements = []



		sigma = 0

		while sigma <= len(self.k_mers_names[k]) - 1 :

			liste_elements.append(self.k_mers_elements[k][sigma][0]) # chaque element est un vecteur de longeur un

			if self.k_mers_names[k][sigma].split('.')[0] not in liste_noms : #traitement des noms sous format rpy2-liste

				liste_noms.append(self.k_mers_names[k][sigma].split('.')[0])

				intervalles.append(sigma) #note la position ou le nom change

			sigma = sigma + 1





		sortie_dict = dict.fromkeys(liste_noms)

		if len(intervalles) > 1 :

			sowil = 0

			while sowil <= len(intervalles) - 1 :

				min = intervalles[sowil]

				if sowil < len(intervalles) - 1 :

					max = intervalles[sowil + 1]

					sortie_dict[liste_noms[sowil]] = liste_elements[min:max]

					if len(sortie_dict[liste_noms[sowil]]) == 1 : #si l'element est de longueur 1, recuperer ce qu'il y a l'interieur de l'element

						sortie_dict[liste_noms[sowil]] = sortie_dict[liste_noms[sowil]][0]

					else : #si l'element est de longueur superieure a 1, rien

						pass

				if sowil == len(intervalles) - 1 :

					max = intervalles[sowil]

					if (max + 1) < len(liste_elements) :

						sortie_dict[liste_noms[sowil]] = liste_elements[max:(max + 2)]

					else :

						sortie_dict[liste_noms[sowil]] = liste_elements[max]

				sowil = sowil + 1

		else :

			sortie_dict[liste_noms[0]] = liste_elements



		return sortie_dict

	def process_monomere(self) :

		sortie_dict = {}

		uruz = 1

		while uruz <= len(self.r_liste) :

			sortie_dict[self.k_mers_names[1][uruz - 1]] = self.r_liste.rx2(uruz)

			uruz = uruz + 1

		return sortie_dict

########################################################################################################################
#R-VECTEUR VERS PYTHON-LISTE
########################################################################################################################



def r_vecteur_vers_python_liste(objet) : #convertir un vecteur r en liste python

	nb_cpu = multiprocessing.cpu_count()



	if nb_cpu >= len(objet) or math.floor((len(objet)/(2*nb_cpu))) == 1 :

		nb_cpu = math.ceil(min(len(objet), nb_cpu) / 10)

	pre_sortie = r_vers_p_liste(objet, nb_cpu)

	gc.collect()

	#multiprocessing.set_start_method('spawn', force=True)

	pool = multiprocessing.Pool(len(pre_sortie.index))

	result_tranche = pool.map(pre_sortie.process_colonne_dataframe_virtuel, [k for k in range(1, len(pre_sortie.k_mers) + 1)]) #le premier k-mer a pour indice 1 et le dernier k-mer a pour indice la longueur de l'index + 1

	pool.close()



	sortie = []

	sowil = 0

	while sowil <= len(result_tranche) - 1 :

		sortie.extend(result_tranche[sowil])

		sowil = sowil + 1

	sortie = sortie[0:len(sortie) - pre_sortie.ajout]

	return sortie







class r_vers_p_liste :

	def __init__(self, liste, nb_cpu) :

		self.r_liste = liste

		self.index = {}

		self.k_mers = {}

		if nb_cpu > 1 :

			k = 1

			minimum = 0

			maximum = math.ceil(len(self.r_liste) / nb_cpu)

			etendue = math.ceil(len(self.r_liste) / nb_cpu)

			while k <= nb_cpu :

				if maximum > len(self.r_liste) :

					maximum = len(self.r_liste)

				self.index[k] = list(range(minimum, maximum, 1))

				self.k_mers[k] = liste[minimum:maximum]



				minimum = maximum

				maximum = minimum + etendue


				k = k + 1

			self.ajout = len(self.k_mers[len(self.k_mers) - 1]) - len(self.k_mers[len(self.k_mers)])

		else :

			k = 1

			self.index[k] = list(range(1, len(self.r_liste), 1))

			self.k_mers[k] = liste

			self.ajout = 0



		r_add_vecteur = robjects.globalenv['r_add_vecteur']

		self.k_mers[len(self.k_mers)] = r_add_vecteur(self.k_mers[len(self.k_mers)], self.ajout)

		gc.collect()



	def process_colonne_dataframe_virtuel(self, k) :

		liste = []

		nombre_lignes = len(self.k_mers[k])

		theta = 0

		while theta <= nombre_lignes - 1 : #BOUCLE 3

			liste.append(self.k_mers[k][theta])

			theta = theta + 1

			gc.collect()

		return liste





########################################################################################################################
#R-DATAFRAME VERS PYTHON-DATAFRAME
########################################################################################################################


def r_dataframe_vers_python_dataframe(r_dataframe) :

	if r_dataframe.nrow * r_dataframe.ncol == 0 : #si le dataframe est vide (aucune colonne ou aucune ligne)

		sortie = pandas.DataFrame()

	else :

		r_replace_na = robjects.globalenv['r_replace_na']

		r_dataframe = r_replace_na(r_dataframe)

		objet = r_vers_p_dataframe(r_dataframe)

		sortie = []

		k = 1

		while k <= max(list(objet.index.keys())) : #parcours de toutes les tranches de l'objet #BOUCLE 1 qui parcourt les tranches

			liste = objet.parcours_tranche(k)

			if isinstance(liste[0], list) == True :

				sortie.extend(liste)

			else :

				sortie.append(liste)

			gc.collect()

			k = k + 1

		sortie = pandas.DataFrame(sortie).transpose()

		noms_colonnes = []

		for i in range(0, len(r_dataframe.colnames)) :

			noms_colonnes.append(r_dataframe.colnames[i])

		noms_lignes = []

		for j in range(0, len(r_dataframe.rownames)) :

			noms_lignes.append(r_dataframe.rownames[j])

		sortie.set_axis(noms_colonnes, axis = 1, inplace = True)

		sortie.set_axis(noms_lignes, axis = 'index', inplace = True)

		sortie.replace('*-NA-*', numpy.nan, inplace = True)

	return sortie



class r_vers_p_dataframe :

	def __init__(self, dataframe) :

		nb_cpu = multiprocessing.cpu_count()

		self.r_dataframe = dataframe

		self.index = {}

		if nb_cpu > 1 :

			self.nb_tranches = math.ceil(self.r_dataframe.ncol / nb_cpu)

			self.sortie = pandas.DataFrame(index=range(self.r_dataframe.nrow),columns=range(self.r_dataframe.ncol))

			maximum = min(self.r_dataframe.ncol, multiprocessing.cpu_count())

			minimum = 1

			k = 1

			while k <= self.nb_tranches :

				self.index[k] = list(range(minimum, maximum + 1, 1))

				minimum = maximum + 1

				maximum = maximum + multiprocessing.cpu_count()

				if maximum > self.r_dataframe.ncol :

                                    maximum = self.r_dataframe.ncol

				k = k + 1

		if nb_cpu == 1 :

			self.index[1] = list(range(1, self.r_dataframe.ncol, 1))

			self.nb_tranches = 1

	def parcours_tranche(self, k) :

		df = self

		temp = tranche(self.r_dataframe, k, self.index[k]) #preciser les colonnes à cibler

		pool = multiprocessing.Pool(len(self.index[k]))

		borne_inf = 1

		if len(self.index[k])%2 == 0 :

			borne_sup = len(self.index[k])

		if len(self.index[k])%2 != 0 :

			borne_sup = len(self.index[k])

		if borne_inf != borne_sup :

			result_tranche = pool.starmap(temp.process_colonne, [(numero_colonne, self.r_dataframe.nrow) for numero_colonne in range(borne_inf, borne_sup + 1)])

		else :

			result_tranche = temp.process_colonne(borne_inf, self.r_dataframe.nrow)

		return result_tranche

		pool.close()

		gc.collect()


class tranche :

	def __init__(self, dataframe, k, index_k) :

		self.numero_tranche = k

		self.plage_colonnes = index_k

		r_sousdataframe = robjects.globalenv['r_sousdataframe']

		self.temp_r_dataframe = r_sousdataframe(dataframe, start_colonne = self.plage_colonnes[0], fin_colonne = self.plage_colonnes[len(self.plage_colonnes) - 1], force_dataframe = True) #numero des colonnes en R correspond au numero des colonnes dans l'index

		gc.collect()



	def process_colonne(self, numero_colonne, nombre_lignes) :

		temp_column_un = base.as_character(self.temp_r_dataframe.rx2(numero_colonne)) #ce bloc permet de gerer le cas ou le dataframe ne contient qu'une colonne et est donc un vecteur

		if len(temp_column_un) == 1 :

			temp_column_un = base.as_character(self.temp_r_dataframe)

		else :

			pass



		for_evaluation = temp_column_un[0]

		try : #essai de la conversion en entier

			int(for_evaluation)

			temp_column_un = base.as_numeric(temp_column_un)

		except ValueError : #si conversion en entier impossible

			try : #alors essai de la conversion en flotrant

				float(for_evaluation)

				temp_column_un = base.as_double(temp_column_un)

			except ValueError : #si conversion en flottant impossible

				try : #alors essai de la conversion en complexe

					complex(for_evaluation)

					temp_column_un = base.as_complex(temp_column_un)

				except ValueError : #si aucune de ces conversions n'est possible, ce r-vecteur contient des chaines de caracteres

					pass

				pass

			pass

		liste = list(temp_column_un)

		return liste




###########################################################################################################################
###########################################################################################################################
#							PYTHON -> R
###########################################################################################################################
###########################################################################################################################






###########################################################################################################################
#FONCTION PRINCIPALE
###########################################################################################################################



def python_vers_r_(objet) :

    vecteur_type_liste = (dict, list)

    if isinstance(objet, dict) == True :

        sortie = _2_python_dict_vers_r_liste(objet)

    elif isinstance(objet, pandas.core.frame.DataFrame) == True :

        sortie = _2_python_dataframe_vers_r_dataframe(objet)

    elif isinstance(objet, numpy.ndarray) == True :

        sortie = pandas.DataFrame(objet)

        sortie = _2_python_dataframe_vers_r_dataframe(sortie)

        sortie = base.as_matrix(sortie)

    elif isinstance(objet, list) == True :

        sortie = _2_python_liste_vers_r_vecteur(objet)

    elif isinstance(objet, tuple) == True :

        sortie = list(objet)

        sortie = _2_python_liste_vers_r_vecteur(sortie)

    else :

        sortie = objet

        print('\n')

        print(type(objet), ': pas de conversion disponible')

        print('\n')

    return sortie

###########################################################################################################################
#FPYTHON DICTIONNAIRE -> R LISTE
###########################################################################################################################


def _2_python_dict_vers_r_liste(objet) :

	if objet == {} :

		sortie = base.list()

	else :

		creation_r_liste = robjects.globalenv['creation_r_liste']

		pre_sortie = _2_p_dict_vers_r_liste(objet, nb_cpu = 2)

		gc.collect()

		#multiprocessing.set_start_method('spawn', force=True)

		pool = multiprocessing.Pool(pre_sortie.etendue)

		result_parallele = pool.map(pre_sortie.process_k_mers, [k for k in range(1, pre_sortie.nb_kmers + 1)])



		sortie = result_parallele[0]

		ophala = 1

		while ophala <= len(result_parallele) - 1 :

			sortie = base.c(sortie, result_parallele[ophala])

			ophala = ophala + 1

			gc.collect()

	return sortie




class _2_p_dict_vers_r_liste :

	def __init__(self, objet, nb_cpu=4) :

		liste_clefs = list(objet.keys())

		liste_valeurs = list(objet.values())

		if nb_cpu > 1 :

			self.longueur = len(liste_clefs)

			self.etendue = nb_cpu

			self.index = {}

			self.valeurs = {}



			if nb_cpu > 1 :

				self.nb_kmers = math.ceil(self.longueur / self.etendue)

				maximum = min(self.longueur, self.etendue)

				minimum = 1

				k = 1

				while k <= self.nb_kmers :

					self.index[k] = list(range(minimum, maximum + 1, 1))

					self.valeurs[k] = liste_valeurs[(minimum - 1):maximum]



					minimum = maximum + 1

					maximum = maximum + self.etendue

					if maximum > self.longueur :

						maximum = self.longueur

					k = k + 1

			else :

				self.nb_kmers = 1

				self.index[1] = list(range(1, len(objet) + 1, 1))

				self.valeurs[1] = liste_valeurs

	def process_k_mers(self, k) :

		creation_r_liste = robjects.globalenv['creation_r_liste']

		matrice_r_liste = creation_r_liste(self.index[k][0], self.valeurs[k][0])

		raido = 1

		while raido <= len(self.index[k]) - 1 :

			temp_r_liste = creation_r_liste(self.index[k][raido], self.valeurs[k][raido])

			matrice_r_liste = base.c(matrice_r_liste, temp_r_liste)

			raido = raido + 1

			gc.collect()

		return matrice_r_liste




###########################################################################################################################
#PYTHON-DATAFRAME VERS R-DATAFRAME
###########################################################################################################################



def _2_python_dataframe_vers_r_dataframe(python_dataframe) : #fonction d'appel pour la conversion des py-dataframe

	#print('coucou')

	if python_dataframe.size == len(python_dataframe.index) : #si c'est un dataframe d'une colonne

		sortie = list(python_dataframe.iloc[:,0])

		sortie = python_vers_r_(sortie)

		sortie = base.as_data_frame(sortie)

		sortie.colnames = python_dataframe.columns[0]

		sortie.rownames = list(python_dataframe.index)

	elif python_dataframe.size == len(python_dataframe.columns) : #si c'est un dataframe d'une ligne

		sortie = list(python_dataframe.iloc[0,:])

		dict_type_position = {'character' : [], 'numeric' : [], 'bool' : [], 'complex' : []}

		for i, j in enumerate(sortie): #https://stackoverflow.com/questions/30843103/how-to-get-the-index-of-an-integer-from-a-list-if-the-list-contains-a-boolean 09.05.2022 - 15.06.2017 Haresh Shyara

			if isinstance(j, is_numeric) :

				dict_type_position['numeric'].append(i)

			elif isinstance(j, bool) :

				dict_type_position['bool'].append(i)

			elif isinstance(j, complex) :

				dict_type_position['complex'].append(i)

			else :

				dict_type_position['character'].append(i)



		sortie = python_vers_r_(sortie)

		sortie = base.as_data_frame(sortie)

		sortie = base.as_data_frame(base.t(sortie))

		sortie.colnames = list(python_dataframe.columns)

		sortie.rownames = [1]

		r_changetype_column = robjects.globalenv['r_changetype_column']



		alpha = 0

		liste_clefs = list(dict_type_position.keys())



		while alpha <= len(liste_clefs) - 1 :

			beta = 0

			while beta <= len(dict_type_position[liste_clefs[alpha]]) - 1 :

				sortie = r_changetype_column(sortie, type = liste_clefs[alpha], numero_colonne = dict_type_position[liste_clefs[alpha]][beta] + 1)

				beta = beta + 1

				gc.collect()

			alpha = alpha + 1

	else :

		#print('ici')

		objet = _2_p_dataframe_vers_r_dataframe(python_dataframe) #instanciation de l'objet kmerise

		print('ici2')

		sortie = base.data_frame(base.rep('', python_dataframe.shape[0])) #cree un r-dataframe vide avec autant de lignes que le py-dataframe initial

		sortie.colnames = '"*"'

		liste_df = []

		k = 1

		while k <= max(list(objet.index.keys())) : #parcours de toutes les tranches de l'objet #BOUCLE 1 qui parcourt les tranches

			#print(' k :', k)

			temp_liste = objet._2_p_r_parcours_tranche(k) #parcours des tranches (sous-dataframe avec autant de colonne que de coeurs) dans le py-dataframe kmerise. On obtient en sortie les tranches (r-tranches) mais sous forme de r-dataframe

			liste_df.extend(temp_liste)

			k = k + 1

		theta = 0

		while theta <= len(liste_df) - 1 :

			print('theta :', theta)

			sortie = base.cbind(sortie, liste_df[theta]) #concatenation des r-tranches

			gc.collect()

			theta = theta + 1



		del_colonne_vide = robjects.globalenv['del_colonne_vide']

		sortie = del_colonne_vide(sortie)

		sortie.colnames = list(python_dataframe.columns)[0:base.dim(sortie)[1]]

		sortie.rownames = list(python_dataframe.index)

		print('-BYE-')

	return sortie


class _2_p_dataframe_vers_r_dataframe :

    def __init__(self, py_dataframe) :

        nb_cpu = multiprocessing.cpu_count()

        self.py_dataframe = py_dataframe

        self.nb_colonnes = len(self.py_dataframe.columns)

        self.index = {}

        if nb_cpu > 1 :

            self.nb_tranches = math.ceil(self.nb_colonnes / nb_cpu)

            maximum = min(self.nb_colonnes, multiprocessing.cpu_count())

            minimum = 1

            k = 1

            while k <= self.nb_tranches :

                self.index[k] = list(range(minimum, maximum + 1, 1))

                minimum = maximum + 1

                maximum = maximum + multiprocessing.cpu_count()

                if maximum > self.nb_colonnes :

                    maximum = self.nb_colonnes

                else :

                    pass

                k = k + 1

        if nb_cpu == 1 :

            self.index[1] = list(range(1, self.nb_colonnes, 1))

            self.nb_tranches = 1

    def _2_p_r_parcours_tranche(self, k) :

        #multiprocessing.set_start_method('spawn', force=True)

        pool = multiprocessing.Pool(len(self.index[k]))

        print('longueur :', len(self.index[k]))

        print('k1 :', k)

        temp_tranche = _2_tranche(self.py_dataframe, k, self.index[k])

        result_tranche = pool.map(temp_tranche._2_p_r_process_colonne, [j for j in range(0, len(self.index[k]) - 0)])

        pool.close()

        gc.collect()



        return result_tranche




class _2_tranche :

    def __init__(self, dataframe, k, index_k) :

        self.numero_tranche = k

        self.plage_colonnes = index_k

        self.temp_py_dataframe = dataframe.iloc[:,(self.plage_colonnes[0] - 1):self.plage_colonnes[len(self.plage_colonnes) - 1]] #numero des colonnes donnes grace a l'index



    def _2_p_r_process_colonne(self, numero_colonne) :

        temp_colonne = list(self.temp_py_dataframe.iloc[:,numero_colonne])

        try :

            if isinstance(temp_colonne[0], str) == True :

                temp_colonne = rpy2.robjects.StrVector(temp_colonne)

                #print(utils.tail(temp_colonne))

            elif isinstance(temp_colonne[0], is_numeric) == True :

                temp_colonne = rpy2.robjects.FloatVector(temp_colonne)

            #elif isinstance(temp_colonne[0], int) == True :

                #temp_colonne = rpy2.robjects.IntVector(temp_colonne)

            elif isinstance(temp_colonne[0], bool) == True :

                temp_colonne = rpy2.robjects.BoolVector(temp_colonne)

            else :

                temp_colonne = rpy2.robjects.ComplexVector(temp_colonne)

        except (AttributeError, ValueError) :

            liste = []

            sowil = 0

            while sowil <= len(temp_colonne) - 1 :

                liste.append(str(temp_colonne[sowil]))

                sowil = sowil + 1

            temp_colonne = rpy2.robjects.StrVector(liste)

            gc.collect()

        return temp_colonne


###########################################################################################################################
#PYTHON-LISTE VERS R-VECTEUR
###########################################################################################################################



def _2_fonction_typage(a_typer) :

    if isinstance(a_typer[0], str) == True :

        typage = 'chaine'

    elif isinstance(a_typer[0], is_numeric) == True :

        typage = 'numerique'

    elif isinstance(a_typer[0], bool) == True :

        typage = 'booleen'

    elif isinstance(a_typer[0], complex) == True :

        typage = 'complexe'

    else :

        typage = 'inconnu'

    return typage


def _2_python_liste_vers_r_vecteur(test_liste) :

    test = _2_p_liste_vers_r_vecteur(test_liste)

    gc.collect()

    #multiprocessing.set_start_method('spawn', force=True)

    pool = multiprocessing.Pool(len(test.index))

    result_tranche = pool.map(test._2_process_k_mers, [k for k in range(1, len(test.k_mers) + 1)])

    if test.typage == 'chaine' :

        sortie = rpy2.robjects.StrVector([])

    elif test.typage == 'numerique' :

        sortie = rpy2.robjects.FloatVector([])

    elif test.typage == 'booleen' :

        sortie = rpy2.robjects.BoolVector([])

    elif test.typage == 'complexe' :

        sortie = rpy2.robjects.ComplexVector([])

    else :

        return test





    fin = 0

    while fin <= len(result_tranche) - 1 :

        sortie = base.c(sortie, result_tranche[fin])

        fin = fin + 1

        gc.collect()

    return sortie




class _2_p_liste_vers_r_vecteur :

    def __init__(self, p_liste) :

        self.p_liste = p_liste

        self.typage = _2_fonction_typage(p_liste)

        nb_cpu = multiprocessing.cpu_count()

        self.index = {}

        self.k_mers = {}

        if nb_cpu > 1 :

            minimum = 0

            maximum = math.ceil(len(self.p_liste) / nb_cpu)

            uruz = 1

            while uruz <= nb_cpu :

                if maximum > len(self.p_liste) :

                    maximum = len(self.p_liste)

                else :

                    pass

                self.index[uruz] = list(range(minimum, maximum + 1, 1))

                self.k_mers[uruz] = self.p_liste[minimum:maximum + 1]

                if self.k_mers[uruz] == [] : #supprime les k-mers vide

                    del self.k_mers[uruz]

                    del self.index[uruz]

                else :

                    pass

                minimum = maximum + 1

                maximum = minimum + maximum

                uruz = uruz + 1

                gc.collect()

        else :

            self.index[1] = list(range(0, len(self.p_liste) - 1, 1))

            self.k_mers[1] = self.p_liste



    def _2_process_k_mers(self, k) :

        temp = self.k_mers[k]

        try :

            if self.typage == 'chaine' :

                sortie = rpy2.robjects.StrVector(temp)

            elif self.typage == 'numerique' :

            	sortie = rpy2.robjects.FloatVector(temp)

            elif self.typage == 'booleen' :

            	sortie = rpy2.robjects.BoolVector(temp)

            elif self.typage == 'complexe' :

            	sortie = rpy2.robjects.ComplexVector(temp)

            else :

            	sortie = []

            	j = 0

            	while j <= len(self.k_mers[k]) - 1 :

                	sortie.append(str(self.k_mers[k][j]))

                	gc.collect()

                	j = j + 1 ;

            	sortie = rpy2.robjects.StrVector(sortie)

        except (AttributeError, ValueError) :

        	sortie = []

        	j = 0

        	while j <= len(self.k_mers[k]) - 1 :

        		sortie.append(str(self.k_mers[k][j]))

        		j = j + 1 ;

        	sortie = rpy2.robjects.StrVector(sortie)

        	gc.collect()

        return sortie


###########################################################################################################################
###########################################################################################################################
#							AUTOMTATIC CALCUL OF REQUIRED THREAD NUMBER
###########################################################################################################################
###########################################################################################################################


def nombre_coeurs(objet, nb_cpu = multiprocessing.cpu_count()) :

	dict_coeurs = {'cpu' : nb_cpu, 'etendue' : 0}

	dict_coeurs['etendue'] = math.ceil(len(objet) / dict_coeurs['cpu'])

	while dict_coeurs['cpu'] > dict_coeurs['etendue'] :

		dict_coeurs = nombre_coeurs(objet, nb_cpu - 1)

		if dict_coeurs['cpu'] == 1 :

			dict_coeurs['cpu'] = 2

			break

		else :

			pass

	if dict_coeurs['cpu'] == dict_coeurs['etendue'] and dict_coeurs['cpu'] - 1 != 1 :

		dict_coeurs['cpu'] = dict_coeurs['cpu'] - 1

	dict_coeurs['etendue'] = math.ceil(len(objet) / dict_coeurs['cpu'])

	print(dict_coeurs)

	return dict_coeurs
