set.seed(10)
DATE = Sys.Date()

PATH_DATA = 'C:/Users/Felix/DropBox/Felix_Kevin_Han-seniorthesis2015-16/data/'
PATH_SAVE = PATH_DATA

lapply(grep('^source_[:alnum:]+', dir(), value = T), source)
