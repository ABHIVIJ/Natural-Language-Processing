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
		if d :						#some lines are empty in input file
			element = d.split("\t")
			element[1] = getClass(element[1])
			print(element)
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
	if word in dictionary:
		return True
	else:
		return False

def containsHyphen(word):
	if(word.find("-") == -1):
		return False
	else:
		return True

def ifAllCaps(word):
	if(word.upper()):
		return True
	else:
		return False

def ifAllSmall(word):
	if(word.lower()):
		return True
	else:
		return False

def ifNumber(word):
	if(word.lower()):
		return True
	else:
		return False

def getDictionary(filename):
	dictionary = (file(filename).read()).split("\r")
	for i in range(len(dictionary)):
		dictionary[i] = dictionary[i][1:]
	#print(dictionary[:5])
	return dictionary	

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

def evaluate(data, prob_table, ffunc_list, dict_func, dictionary):
	contigency_matrix = {}
	cont_matrix_row = {"DNA":0, "RNA":0, "protein":0, "cell_type":0, "cell_line":0, "other":0}
	contigency_matrix["predicted"] = cont_matrix_row
	contigency_matrix["actual"] = cont_matrix_row
	contigency_matrix["got_right"] = cont_matrix_row
	
	for e in data:
		pred = predict(e[0], prob_table, ffunc_list, dict_func,  dictionary)
		contigency_matrix["actual"][e[1]] += 1 
		contigency_matrix["predicted"][pred] += 1
		if(e[1] == pred):	
			contigency_matrix["got_right"][e[1]] += 1

	print("Results")
	correct = 0
	total = 0	
	for e in cont_matrix_row:
		correct += contigency_matrix["got_right"][e]
		total += contigency_matrix["actual"][e]
		p = (contigency_matrix["got_right"][e]*1.0)/contigency_matrix["predicted"][e]
		r = (contigency_matrix["got_right"][e]*1.0)/contigency_matrix["actual"][e]	
		f = (2*p*r)/(p+r)
		print(e)
		print(contigency_matrix["got_right"][e])
		print(contigency_matrix["predicted"][e])
		print(contigency_matrix["actual"][e])
		print("Precision")
		print(p)
		print("Recall")
		print(r)
		print("F-Measure")
		print(f)		
	
	acc = (1.0 * correct)/total
	print("Accuracy")
	print(acc)
	
#train_data = getData(file("Genia4ER_train.txt").read())

#print(train_data[:5])

dictionary = getDictionary("comWord.txt")	

#print(featureCount(ifPresentInDict, train_data, dictionary))
#print(featureCount(ifAllCaps, train_data))
#print(featureCount(ifAllSmall, train_data))
#print(featureCount(ifNumber, train_data))
#print(featureCount(containsHyphen, train_data))

ffunc_list = [ifPresentInDict, ifAllCaps, ifAllSmall, ifNumber, containsHyphen]
dict_func = [ifPresentInDict]

#prob_table = getFeatureProb(train_data, ffunc_list, dict_func, dictionary)

#print(predict("Number", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("GR", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("lymphocytes", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("interleukin-2", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("mRNA", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("tightly", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("atomic", prob_table, ffunc_list, dict_func,  dictionary))
#print(predict("and", prob_table, ffunc_list, dict_func,  dictionary))

print("Evaluation on training data")
#evaluate(train_data, prob_table, ffunc_list, dict_func, dictionary)

test_data = getData(file("Genia4ER_test.txt").read())

print("Evaluation on test data")
evaluate(test_data, prob_table, ffunc_list, dict_func, dictionary)

