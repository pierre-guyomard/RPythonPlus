stacking_liste_ggplot <- function(liste) { #prend en entree une liste contenant des dataframe stackes ; retourne un dataframe unique; stacking d'une liste contenant des dataframes stackes

	##print('--log-fonction : stacking_liste_ggplot')
	
	df_sortie <- liste[[1]]

	sigma <- 2

	while (sigma <= length(liste)) {
                
		df_sortie <- rbind(df_sortie, liste[[sigma]])
                
		sigma <- sigma + 1
                
	}
        
	return(df_sortie)
	
	##print('--log-fonction : fait')
        
}

colbind_liste <- function(liste) { #cree un dataframe a partir de vecteurs de taille differente  
        
	##print('--log-fonction : colbind_liste')
	
	max_dim <- 0
        
	i <- 1
        
	while (i <= length(liste)) {
                
		if (length(liste[[i]]) > max_dim) {

			max_dim <- length(liste[[i]])
                        
		}
                
		i <- i + 1
                
	}
        
	sortie <- data.frame(x = rep(NA, max_dim + 1))
        
	i <- 1
        
	while (i <= length(liste)) {
                
		ad <- max_dim - length(liste[[i]])
                
		temp_v <- rep(NA, ad)
                
		sortie <- cbind(sortie, c(median(liste[[i]]), liste[[i]], temp_v))
                
		i <- i + 1
                
	}
        
	sortie <- sortie[, -1]
        
	colnames(sortie) <- names(liste)
        
	rownames(sortie)[1] <- 'mediane'
        
	sortie <- as.data.frame(t(sortie))
        
	sortie <- sortie[order(sortie$mediane), ]
        
	sortie <- as.data.frame(t(sortie))
        
	sortie <- sortie[-1, ]
        
	return(sortie)
	
	##print('--log-fonction : fait')
        
}              



colbind_df <- function(dataframe, vecteur, nom = '*ajout*', position='apres') { #cree un dataframe a partir de vecteurs de taille differente  
	
	
	
	ad <- dim(dataframe)[1] - length(vecteur)
	
	
	
	theta_nom <- colnames(dataframe)

	if (ad > 0) {
                
		vecteur <- c(vecteur, rep(NA, ad))
		
		cat(paste('\n  <NA> completed as arguments imply differing number of rows: ', dim(dataframe)[1], ', ', length(vecteur) - ad, sep = ''))
		
		cat('\n')

	}

	else if (ad == 0) {

		dataframe <- dataframe

	}

	else { #ad < 0
	
		temp <- rep(NA, dim(dataframe)[2])
		
		i <- 1
		
		while (i <= -ad ) {
		
			dataframe <- rbind(dataframe, temp)
			
			i <- i + 1
			
		}

		#vecteur <- vecteur[1:dim(dataframe)[1]]
		
		cat(paste('\n  <NA> completed as arguments imply differing number of rows: ', dim(dataframe)[1], ', ', length(vecteur) - ad, sep = ''))
		
		cat('\n')

	}
	
	if (position == 'apres') {	
	
		dataframe <- cbind(dataframe, vecteur)
		
		colnames(dataframe) <- c(theta_nom, nom)	
	
	}
	
	else {
	
		dataframe <- cbind(vecteur, dataframe)
		
		colnames(dataframe) <- c(nom, theta_nom)
		
	}

	cat('\n')
	
	return(dataframe)
        
}



rowbind_df <- function(dataframe, vecteur, nom = '*ajout*', position='dessous') { #cree un dataframe a partir de vecteurs de taille differente  
	
	
	
	ad <- dim(dataframe)[2] - length(vecteur)
	
	
	
	theta_nom <- rownames(dataframe)

	if (ad > 0) {
                
		vecteur <- c(vecteur, rep(NA, ad))
		
		cat(paste('\n  <NA> completed as arguments imply differing number of columns: ', dim(dataframe)[2], ', ', length(vecteur) - ad, sep = ''))
		
		cat('\n')

	}

	else if (ad == 0) {

		dataframe <- dataframe

	}

	else { #ad < 0
	
		assign('NO_U_V_E_AU', rep(NA, dim(dataframe)[1]))
		
		i <- 1
		
		while (i <= -ad ) {
		
			dataframe <- cbind(dataframe, NO_U_V_E_AU)
			
			i <- i + 1
			
		}

		#vecteur <- vecteur[1:dim(dataframe)[1]]
		
		cat(paste('\n  <NA> completed as arguments imply differing number of columns: ', dim(dataframe)[2], ', ', length(vecteur) - ad, '\n  - new column named : "NO_U_V_E_AU"', sep = ''))
		
		cat('\n')

	}
	
	if (position == 'dessous') {	
	
		dataframe <- rbind(dataframe, vecteur)
		
		rownames(dataframe) <- c(theta_nom, nom)	
	
	}
	
	else {
	
		dataframe <- rbind(vecteur, dataframe)
		
		rownames(dataframe) <- c(nom, theta_nom)
		
	}

	cat('\n')
	
	return(dataframe)
        
}        



destack_dataframe <- function(stacked_dataframe, colonne_facteur, sortie_liste = FALSE) { #prend en entree un dataframe stacke et renvoie ou une liste ou un dataframe
        
	#print('--log-fonction : destack_dataframe')
	
	facteur <- levels(as.factor(stacked_dataframe[, colonne_facteur]))
        
	liste <- list()
        
	theta <- 1
        
	while (theta <= length(facteur)) {
                
		temp <- subset(stacked_dataframe, stacked_dataframe[, colonne_facteur] == facteur[theta])
                
		liste[[theta]] <- temp[, 2]
                
		names(liste)[theta] <- facteur[theta]
                
		theta <- theta + 1
                
	}
        
	if (sortie_liste == TRUE) {
                
		return(liste)
                
	}
        
	if (sortie_liste == FALSE) {
                
		return(colbind_liste(liste))
                
	}
	
	#print('--log-fonction : fait')
        
} 
    
stacking <- function(vecteur, name) { #permet de transformer un vecteur en element stacke
        
	#print('--log-fonction : stacking')
	
	temp <- data.frame(x = rep(name, length(vecteur)), y = vecteur)
        
	colnames(temp) <- c('my', 'value')
        
	return(temp)
	
	#print('--log-fonction : fait')
        
}


stacking_subset <- function(facteur, ligne) { #stacke un vecteur selon un facteur

	#print('--log-fonction : stacking_subset')
	
	facteur <- as.vector(as.factor(facteur))

	ligne <- as.numeric(ligne)
        
	ad <- length(facteur) - length(ligne)

	if (ad < 0) {

		ad <- (-1)*ad

	}

	if (ad == 0 ) {

		sortie <- data.frame(x = facteur, y = ligne)

	}

	if (ad > 0 ) {

		sortie <- data.frame(x = facteur, y = c(ligne, rep(NA, ad)))

	}

	colnames(sortie) <- c('my', 'value')

	return(sortie)
	
	#print('--log-fonction : fait')

}


              
comptage <- function(dataframe, debut, fin, numero_colonne) { #retourne une table de comptage à partir d'un dataframe bidimensionnel
        
	#print('--log-fonction : comptage')
	
	sortie <- aggregate(formula = dataframe[, debut] ~ dataframe[, numero_colonne], data = dataframe, FUN = sum)
        
	temp_nom <- c(colnames(dataframe)[numero_colonne], colnames(dataframe)[debut])
	
	
        
	compteur <- debut + 1
        
	while (compteur <= fin ) {
                
		temp_aggregate <- aggregate(formula = dataframe[, compteur] ~ dataframe[, numero_colonne], data = dataframe, FUN = sum)

		#print(compteur)
                
		temp_nom <- c(temp_nom, colnames(dataframe)[compteur])

		

		
                
		sortie <- colbind_df(sortie, as.numeric(temp_aggregate[, 2]))
                
		colnames(sortie) <- temp_nom

		compteur <- compteur + 1
                
	}
	
	
	
	#print('comptage - fait')
        
	return(sortie)
        
}

empilement <- function(dataframe) { #stacke un dataframe dans lequel les donnees sont organisées par colonne

	#print('--log-fonction : empilement')
	
	df_sortie <- stacking(dataframe[, 1], colnames(dataframe)[1])
	
	colnames(df_sortie) <- c('my', 'value')
	
	compteur <- 2
	
	while (compteur <= dim(dataframe)[2] ) {

		df_sortie <- rbind(df_sortie, stacking(dataframe[, compteur], colnames(dataframe)[compteur]))
		
		compteur <- compteur + 1
		
	}
	
	return(df_sortie)
	
	#print('--log-fonction : fait')
	
}
		
liste_vers_dataframe <- function(liste) { #prend en entree une liste de vecteurs et renvoie un dataframe

	df_sortie <- as.data.frame(liste[[1]])
	
	i <- 2
	
	while (i <= length(liste)) {
	
		df_sortie <- colbind_df(df_sortie, liste[[i]])
		
		i <- i + 1
		
	}
	
	colnames(df_sortie) <- names(liste)
	
	return(df_sortie)
	
	#print('--log-fonction : fait')
	
}

r_subset <- function(dataframe, numero_colonne, critere, facteur) {#fonction 'subset' dans R

	numero_colonne <- as.character(numero_colonne)
	
	critere <- as.character(critere)
	
	facteur <- as.character(facteur)
	
	
	
	chaine <- paste("subset(dataframe, dataframe[ ,", numero_colonne, "] ", critere, facteur, ")" , sep= "" )
	
	sortie <- eval(parse(text=chaine))
	
	return(sortie)
	
}

r_subvecteur <- function(dataframe, critere, facteur) {

	chaine <- paste("dataframe[dataframe ", critere, ' ', facteur, "]", sep= "" )
	
	sortie <- eval(parse(text=chaine))
	
	return(sortie)
	
}

r_position_vecteur <- function(vecteur, start, stop) {

	sortie <- vecteur[start:stop]
	
	return(sortie)
	
}

r_dataframe_order <- function(dataframe, var, numero, decroissant=FALSE) {

	if (var == 'colonne') {

		chaine <- paste('dataframe[order(dataframe[, ', numero, '], decreasing = ', decroissant, '), ]', sep= '')
		
	}
	
	else {
	
		chaine <- paste('dataframe[, order(dataframe[', numero, ', ], decreasing = ', decroissant, ')]', sep= "")
	
	}
	
	sortie <- eval(parse(text=chaine))
	
	return(sortie)
	
}

r_changetype_column <- function(dataframe, type, numero_colonne) {

	chaine <- paste("dataframe[, numero_colonne] <- as.", type, "(dataframe[, numero_colonne])" , sep= "" )
	
	eval(parse(text=chaine))
	
	return(dataframe)
	
}

r_add_vecteur <- function(liste, nombre) {

	liste <- c(liste, rep("*!*!*!*", nombre))

	return(liste)

}

r_sousdataframe <- function(dataframe, start_ligne='', fin_ligne='', start_colonne='', fin_colonne='', force_dataframe=FALSE) {

	if (fin_ligne == '') {
	
		if (start_ligne == '') {
		
			start_ligne <- 1
			
			fin_ligne <- dim(dataframe)[1]
			
		}
		
		else {#si start_ligne != NULL et fin_ligne == NULL  df[3, ]
		
			fin_ligne <- start_ligne
		
		}
		
		
	}
	
	else { #si fin_colonne != ''
	
		if (start_ligne == '' && fin_ligne > 0) {
		
			start_ligne <- 1
			
		}
		
		if (start_ligne == '' && fin_ligne < 0) {
		
			start_ligne <- -1
			
		}
		
	} 
	
	if (fin_colonne == '') {
	
		if (start_colonne == '') {
		
			start_colonne <- 1
			
			fin_colonne <- dim(dataframe)[2]
			
		}
		
		else {
		
			fin_colonne <- start_colonne
			
		}
		
	}
	
	else { #si fin_colonne != ''
	
		if (start_colonne == '' && fin_colonne > 0) {
		
			start_colonne <- 1			
			
		}
		
		if (start_colonne == '' && fin_colonne < 0) {
		
			start_colonne <- -1			
			
		}
	
	}
	
	
		
	
	sortie <- dataframe[start_ligne:fin_ligne, start_colonne:fin_colonne]
	
	
	
	if (force_dataframe == TRUE) {
	
		sortie <- as.data.frame(sortie)
		
		colnames(sortie) <- colnames(dataframe)[start_colonne:fin_colonne]
	
	}
	
	

	return(sortie)

}

creation_r_liste <- function(clef, valeur) {

	liste <- list(valeur)

	names(liste) <- c(clef)

	return(liste)

}

del_colonne_vide <- function(dataframe) {

	dataframe <- dataframe[, (2:(dim(dataframe)[2]))]

	return(dataframe)

}

concatenation_r_vecteur <- function(matrice, elongation) {

	matrice <- c(matrice, elongation)

	return(matrice)

}

r_liste_vers_r_vecteur <- function(r_liste) {

	r_vecteur_sortie <- c()

	uruz <- 1
	
	while (uruz <= length(r_liste)) {
	
		r_vecteur_sortie <- c(r_vecteur_sortie, r_liste[[uruz]])
	
		uruz <- uruz + 1
	
	}
	
	return(r_vecteur_sortie)
	
}

r_replace_na <- function(dataframe_vecteur) { #remplace les NA dans un r_dataframe ou un r_vecteur

	if ( class(dataframe_vecteur) == 'data.frame' ) {
	
		library(dplyr, warn.conflicts = FALSE)

		dataframe_vecteur <- mutate_all(dataframe_vecteur, ~replace(., is.na(.), '*-NA-*'))
		
	}
	
	else {
	
		dataframe_vecteur[is.na(dataframe_vecteur)] <- '*-NA-*'
		
	}
	
	return(dataframe_vecteur)
	
}

r_changetype_matrix_elements <- function(matrix, type) { #permet de changer le type des elements dans une matrice. Depuis https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix 05.07.2022 answered Sinan U, Jan 2, 2021 at 22:02

	if (type == 'str') {
	
		type <- 'character'
		
	}
	
	else if (type == 'int' || type == 'float') {
	
		type <- 'numeric'
		
	}
	
	else {
	
		type <- type
		
	}	
	
	chaine <- paste("apply(matrix, 2, as.", type, ")", sep = "")
	
	sortie <- eval(parse(text=chaine))
	
	return(sortie)
	
}

r_outliers_dataframe <- function(dataframe, numero_colonne, action) { #https://statsandr.com/blog/outliers-detection-in-r/ 03.08.2022



	outliers <- boxplot.stats(dataframe[ , numero_colonne])$out
	
	
	
	if (length(outliers) > 0) {
	
		outliers_index <- which(dataframe[ , numero_colonne] %in% c(outliers))
	
		#dataframe[outliers_index, numero_colonne]
	
		if (action == 'remove') {
	
			dataframe <- dataframe[-outliers_index, ]
		
		}
	
		else if (action == 'return') {
	
			dataframe <- outliers_index
		
		}
	
		else {
	
			omega <- 1
		
			while (omega <= length(outliers_index)) {
		
				temp <- dataframe[outliers_index[omega], numero_colonne]
	
				chaine <- paste('temp ', action, sep = "")
			
				dataframe[outliers_index[omega], numero_colonne] <- eval(parse(text=chaine))
			
				omega <- omega + 1
			
			}	
		
		}
		
	}
		
	else {
	
		temp <- dataframe
		
		dataframe <- temp 
	
	}
	
	return(dataframe)
	
}



r_outliers_vecteur <- function(vecteur, action) { ##https://statsandr.com/blog/outliers-detection-in-r/ 03.08.2022



	outliers <- boxplot.stats(vecteur)$out
	
	
	
	if (length(outliers) > 0) {
	
		outliers_index <- which(vecteur %in% c(outliers))
	
		#dataframe[outliers_index, numero_colonne]
	
		if (action == 'remove') {
	
			vecteur <- vecteur[-outliers_index]
		
		}
	
		else if (action == 'return') {
	
			vecteur <- outliers_index
		
		}
	
		else {
	
			omega <- 1
		
			while (omega <= length(outliers_index)) {
		
				temp <- vecteur[outliers_index[omega]]
	
				chaine <- paste('temp ', action, sep = "")
			
				vecteur[outliers_index[omega]] <- eval(parse(text=chaine))
			
				omega <- omega + 1
			
			}	
		
		}
		
	}
		
	else {
	
		temp <- vecteur
		
		vecteur <- temp 
	
	}
	
	return(vecteur)
	
}

r_dataframe_vers_r_liste <- function(dataframe) { #prend en entree un dataframe et renvoie une liste dont les entree sont les colonnes

	sortie <- list()

	hagalaz = 1
	
	while (hagalaz <= dim(dataframe)[2]) {
	
		sortie[[hagalaz]] <- dataframe[ , hagalaz]
		
		hagalaz <- hagalaz + 1
		
	}
	
	names(sortie) <- colnames(dataframe)
	
	return(sortie)
	
}

r_retrieve_pval <- function(test_stat) {

	sortie <- test_stat$p.value
	
	return(sortie)
	
}
	
		
r_chaine_vers_graphe <- function(chaine, x_min=-10, x_max=10, y_min = -10, y_max = 10, add=FALSE, log_invert=FALSE, is_log=FALSE, axe_abcisses=FALSE, axe_ordonnees=FALSE, is_identity=FALSE) { #https://stackoverflow.com/questions/23648059/defining-a-curve-by-a-string





        fonction_name <- chaine

        if (is_identity == TRUE) {

        	fonction_name = 'x'

        }



        if (is_log == TRUE) { #remplace les 'log' dans le nom de la fonction par des 'ln'

                fonction_name <- stringr::str_replace_all(fonction_name, "log", "ln")

        }

        if (log_invert == TRUE) { #si l'argument du logarithme est négatif, il faut inverser l'axe

                temp <- -1*x_max

                x_max <- x_min

                x_min <- temp

        }



        x_maximum <- x_max

        x_minimum <- x_min

        if (x_max == 0 || abs(x_max) == Inf) {

                x_maximum <- 1

        }

        if (x_min == 0 || abs(x_min) == Inf ) {

                x_minimum <- 1

        }



        slides <- round(abs(101*x_maximum*x_minimum), digits = 0) #nombre de points à tracer

        evaluer <- eval(substitute(curve(qq, xlim = c(x_min, x_max), ylim = c(y_min, y_max), n = slides), list(qq = parse(text = chaine)[[1]])))



        #plot(evaluer, type = 'l', xlab = 'x', ylab = 'y', main = paste('y = ', fonction_name, sep =''), col = 'blue')



        if ( x_min < 0 & 0 < x_max || axe_abcisses == TRUE) {

                if ( abs(0 - x_min) > 5 || abs(0 - x_max) > 5 || axe_abcisses ==  TRUE) {

                        abline(h = 0, lwd=3, lty=2, col = '#216045')

                }

        }

        if ( min(evaluer$x) < 0 & 0 < max(evaluer$x) || axe_ordonnees == TRUE ) {

                if ( abs(0 - min(evaluer$x)) > 5 || abs(0 - max(evaluer$x)) > 5 || axe_ordonnees == TRUE ) {

                        abline(v = 0, lwd=3, lty=2, col = '#216045')

                }

        }

        return(evaluer)

}










	
