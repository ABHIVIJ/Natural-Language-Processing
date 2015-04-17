def getClass(entity):
	if(entity == "O"):
		return "other"
	else :
		return entity[2:]

def getData(text):
	data = text.split("\n")
	train_data = []
	element = []
	for d in data:
		if d :						#some lines are empty in input file
			element = d.split("\t")
			element[1] = getClass(element[1])
		train_data.append(element)
	return train_data

def featureCount(featureFunction, data, *arg):
	count = {"DNA":0, "RNA":0, "protein":0, "cell_type":0, "cell_line":0, "other":0}
	for e in data:
		if(arg and featureFunction(arg[0], e[0])):	#for isPresentInDict function
			count[e[1]] += 1
		elif((not arg) and featureFunction(e[0])):	#for other feature functions
			count[e[1]] += 1
	return count			

def ifPresentInDict(dictionary, word):
	if word in dictionary:
		return True
	else:
		return False

def containsHyphen(word):
	if(word.find("-") == -1):
		return False
	else:
		return True

def getDictionary(filename):
	dictionary = (file(filename).read()).split("\r")
	for i in range(len(dictionary)):
		dictionary[i] = dictionary[i][1:]
	#print(dictionary[:5])
	return dictionary	

train_data = getData(file("Genia4ER_train.txt").read())

#print(train_data[:5])

dictionary = getDictionary("comWord.txt")

print(featureCount(ifPresentInDict, train_data, dictionary))
print(featureCount(containsHyphen, train_data))



