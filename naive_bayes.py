import math

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
		if(d and ("###MEDLINE:" not in d)) :		#some lines are empty in input file and some lines in test file contains ###MEDLINE:
			element = d.split("\t")
			element[1] = getClass(element[1])

		train_data.append(element)
	return train_data

def featureCount(featureFunction, data, *arg):
	count = {"DNA":1, "RNA":1, "protein":1, "cell_type":1, "cell_line":1, "other":1}	#Laplace Smoothing
	for e in data:
		if(arg and featureFunction(arg[0], e[0])):	#for isPresentInDict function
			count[e[1]] += 1
		elif((not arg) and featureFunction(e[0])):	#for other feature functions
			count[e[1]] += 1
	return count			

def ifPresentInDict(dictionary, word):
	return word in dictionary

def containsHyphen(word):
	return "-" in word

def ifAllCaps(word):
	return word.isupper()

def ifAllSmall(word):
	return word.islower()

def ifNumber(word):
	return word.isdigit()

def getAddData(filename):					#manually put a space in first line of both additional data files to use this function
	text = (file(filename).read()).split("\r")
	for i in range(len(text)):
		text[i] = text[i][1:]
	#print(dictionary[:5])
	return text	

def computeTotalCount(count, entity):
	total = 0
	for f in count:
		total += f[entity]
	return total

def getFeatureProb(train_data, ffunc_list, dict_func, dictionary):
	count = []
	for f in ffunc_list:					#add a new row to count_table for each feature
		if f in dict_func:
			count.append(featureCount(f, train_data, dictionary))
		else:
			count.append(featureCount(f, train_data))

	total_count = {} 
	for e in count[0]:
		total_count[e] = computeTotalCount(count, e)	

	for f in count:						#finding probability
		for e in total_count:
			f[e] = math.log(f[e]*1.0/total_count[e])#Storing the log, so that addition can be done and underflow prevented

	return count

def findProb(word, prob, featureFunc, *arg):
	if(arg and featureFunc(arg[0], word)):
		return prob	
	elif((not arg) and featureFunc(word)):
		return prob
	else:
		return {}

def dictSum(prob, featureProb):
	for e in featureProb:
		prob[e] += featureProb[e]
	return prob
	
def dictMaxKey(prob):
	k = list(prob.keys())
	v = list(prob.values())
	return k[v.index(max(v))]

def predict(word, prob_table, featureFuncList, dict_func, dictionary):
	prob = {"DNA":0, "RNA":0, "protein":0, "cell_type":0, "cell_line":0, "other":0}
	for i in range(len(prob_table)):
		if(featureFuncList[i] in dict_func):
			if(findProb(word, prob_table[i], featureFuncList[i], dictionary)):
				prob = dictSum(prob, prob_table[i])
		else:
			if(findProb(word, prob_table[i], featureFuncList[i])):
				prob = dictSum(prob, prob_table[i])
	#print(word, prob)	
	return dictMaxKey(prob)
	
def getRow(keys):
	row = {}
	for k in keys:
		row[k] = 0
	return row

def evaluate(data, prob_table, ffunc_list, dict_func, dictionary):
	contigency_matrix = {}
	cont_matrix_row = ["DNA", "RNA", "protein", "cell_type", "cell_line", "other"]
	contigency_matrix["predicted"] = getRow(cont_matrix_row)
	contigency_matrix["actual"] = getRow(cont_matrix_row)
	contigency_matrix["got_right"] = getRow(cont_matrix_row)

	for e in data:
		pred = predict(e[0], prob_table, ffunc_list, dict_func,  dictionary)
		contigency_matrix["actual"][e[1]] += 1 
		contigency_matrix["predicted"][pred] += 1
		if(e[1] == pred):	
			contigency_matrix["got_right"][e[1]] += 1

	print(contigency_matrix)

	print("Results\n")
	correct = 0
	total = 0	
	for e in cont_matrix_row:
		correct += contigency_matrix["got_right"][e]
		total += contigency_matrix["actual"][e]
		print("\n")
		print(e)
		#print(contigency_matrix["got_right"][e])
		#print(contigency_matrix["predicted"][e])
		#print(contigency_matrix["actual"][e])
		print("\nPrecision")
		if(contigency_matrix["predicted"][e] != 0):
			p = (contigency_matrix["got_right"][e]*1.0)/contigency_matrix["predicted"][e]
			print(p)
		else:
			print("Undefined")
		print("\nRecall")
		if(contigency_matrix["actual"][e] != 0):
			r = (contigency_matrix["got_right"][e]*1.0)/contigency_matrix["actual"][e]
			print(r)
		else:
			print("Undefined")
		print("\nF-Measure")
		if(contigency_matrix["actual"][e] != 0 or contigency_matrix["predicted"][e] != 0):
			f = (2.0*contigency_matrix["got_right"][e])/(contigency_matrix["actual"][e]+contigency_matrix["predicted"][e])
			print(f)
		else:
			print("Undefined")		
	
	acc = (1.0 * correct)/total
	print("\nAccuracy")
	print(acc)

def getsubstrFeatures(text):
	substr = text.split('\n')
	return substr
	
train_data = getData(file("Genia4ER_train.txt").read())

#print(train_data[:5])

dictionary = getAddData("comWord.txt")
#print(dictionary[:10])

substrFeatures = getAddData("prefix_suffix list.txt")
#print(substrFeatures[:10])	

#print(featureCount(ifPresentInDict, train_data, dictionary))
#print(featureCount(ifAllCaps, train_data))
#print(featureCount(ifAllSmall, train_data))
#print(featureCount(ifNumber, train_data))
#print(featureCount(containsHyphen, train_data))

ffunc_list = [ifPresentInDict, ifAllCaps, ifAllSmall, ifNumber, containsHyphen]
dict_func = [ifPresentInDict]

prob_table = getFeatureProb(train_data, ffunc_list, dict_func, dictionary)

#print(predict("Number", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("GR", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("lymphocytes", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("interleukin-2", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("mRNA", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("tightly", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("atomic", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("and", prob_table, ffunc_list, dict_func,  dictionary))

print("\nEvaluation on training data\n")
evaluate(train_data, prob_table, ffunc_list, dict_func, dictionary)

test_data = getData(file("Genia4ER_test.txt").read())

print("\nEvaluation on test data\n")
evaluate(test_data, prob_table, ffunc_list, dict_func, dictionary)

